//! A tool for interacting with the POST endpoints of Ensembl REST API.
//!  * Spawns an async task that repeatedly polls requests made to its [Client] objects.
//!  * Bundles those requests and posts them to the Ensembl endpoints, asyncronously returning [serde::Deserialize] objects representing the result.
//! ```
//! # tokio::runtime::Builder::new_current_thread()
//! #       .enable_all()
//! #       .build()
//! #       .unwrap()
//! #       .block_on(async {
//! use rs_embl::{Getter, VEPAnalysis};
//!
//! //Wrap in an async tokio runtime
//!
//! let v = Getter::<VEPAnalysis>::new();
//! let handles: Vec<_> = ["3:g.46373453_46373484del", "10:g.72346580_72346583dup"]
//!     .iter()
//!     .map(|id| {
//!         let v = v.client();
//!         tokio::spawn(async move { v.get(id.to_string()).await })
//!     })
//!     .collect();
//! for h in handles.into_iter() {
//!     let vep: Option<VEPAnalysis> = h.await.unwrap();
//!     println!("{:#?}", vep.unwrap());
//! }
//! # });
//! ```

use serde::{de::DeserializeOwned, Deserialize, Serialize};
use std::{
    collections::HashMap,
    sync::{
        atomic::{AtomicBool, Ordering},
        Arc,
    },
};
use tokio::sync::mpsc;
use tokio::time::{sleep, Duration};

/// The minimum time between post operations.
pub const WAIT_DELAY: Duration = Duration::from_millis(500);
const ENSEMBL_SERVER: &str = r#"https://rest.ensembl.org"#;
/// Encapsulates Ensembl REST API calls to allow multiple entries to be condensed into a single POST request.
///  * This will spawn a new asyncronous task that will periodically poll for new requests and handle them.
///  * The task will abort when the Getter object is dropped.
/// ```
/// # tokio::runtime::Builder::new_current_thread()
/// #       .enable_all()
/// #       .build()
/// #       .unwrap()
/// #       .block_on(async {
/// use rs_embl::{Getter, VEPAnalysis};
///
/// // Wrap in an async tokio runtime
///
/// // Instaniate a new Getter object, VEPAnalysis implements serde::Deserialize from the
/// // json returned from Ensembl and implements the EnsemblPostEndpoint trait.
/// let v = Getter::<VEPAnalysis>::new();
/// let handles: Vec<_> = ["3:g.46373453_46373484del", "10:g.72346580_72346583dup"]
///     .iter()
///     .map(|id| {
///         let v = v.client();
///         tokio::spawn(async move { v.get(id.to_string()).await })
///     })
///     .collect();
/// for h in handles.into_iter() {
///     let vep: Option<VEPAnalysis> = h.await.unwrap();
///     println!("{:#?}", vep.unwrap());
/// }
/// # });
/// ```
#[derive(Debug)]
pub struct Getter<T: EnsemblPostEndpoint + Send + DeserializeOwned> {
    is_alive: Arc<AtomicBool>,
    tx: mpsc::Sender<(String, tokio::sync::oneshot::Sender<T>)>,
}
impl<T: EnsemblPostEndpoint + Send + DeserializeOwned> Drop for Getter<T> {
    fn drop(&mut self) {
        self.is_alive.store(false, Ordering::Relaxed);
    }
}
impl<T: 'static + EnsemblPostEndpoint + Send + DeserializeOwned> Default for Getter<T> {
    fn default() -> Self {
        Self::new()
    }
}
impl<T: 'static + EnsemblPostEndpoint + Send + DeserializeOwned> Getter<T> {
    /// Create a new Getter object to return T from Enseble REST endpoint.
    pub fn new() -> Self {
        let (tx, mut rx) = mpsc::channel::<(String, tokio::sync::oneshot::Sender<T>)>(500);
        let is_alive = Arc::new(AtomicBool::new(true));
        let client = reqwest::Client::new();
        {
            let is_alive = Arc::clone(&is_alive);
            tokio::spawn(async move {
                while is_alive.load(Ordering::Relaxed) {
                    sleep(WAIT_DELAY).await;
                    let mut gets = HashMap::new();
                    while let Ok((k, v)) = rx.try_recv() {
                        gets.insert(k, v);
                    }
                    Self::process(gets, &client).await;
                }
                rx.close();
                let mut gets = HashMap::new();
                while let Some((k, v)) = rx.recv().await {
                    gets.insert(k, v);
                }
                Self::process(gets, &client).await;
            });
        }
        Self { is_alive, tx }
    }
    async fn process(
        mut input: HashMap<String, tokio::sync::oneshot::Sender<T>>,
        client: &reqwest::Client,
    ) {
        if input.is_empty() {
            return;
        }
        let ids: Vec<&str> = input.keys().map(|s| s.as_str()).collect();
        let payload = T::payload_template().replace(r"{ids}", &json::stringify(ids));
        let values = client
            .post(String::from(ENSEMBL_SERVER) + T::extension())
            .header("Content-Type", "application/json")
            .header("Accept", "application/json")
            .body(payload)
            .send()
            .await
            .unwrap()
            .text()
            .await
            .unwrap();
        let outputs: Vec<T> = if let Ok(outputs) = serde_json::from_str(&values) {
            outputs
        } else {
            panic!("Failed to parse the following response: {}", values)
        };
        for output in outputs.into_iter() {
            let target = input.remove(output.input()).unwrap();
            let _ = target.send(output); //if the sender's not listening that's it's problem
        }
    }
}
impl<'a, T: 'a + EnsemblPostEndpoint + Send + DeserializeOwned> Getter<T> {
    ///Create a trivially clonable Client that can be sent across async tasks.
    pub fn client(&self) -> Client<'a, T> {
        Client::<T> {
            tx: self.tx.clone(),
            getter: std::marker::PhantomData::<&'a Getter<T>>,
        }
    }
}

/// A Client that can be cloned and sent across async tasks or threads to allow access to the underlying [Getter].
/// * Created by the [Getter::client()] method. [Client::clone()] is equivalent.
/// * Unlike the [Getter], [Client] implements [Send]. Thus, it is usually created in the parent task then passed to workers.
#[derive(Debug, Clone)]
pub struct Client<'a, T: EnsemblPostEndpoint + Send + DeserializeOwned> {
    tx: mpsc::Sender<(String, tokio::sync::oneshot::Sender<T>)>,
    getter: std::marker::PhantomData<&'a Getter<T>>,
}
impl<'a, T: EnsemblPostEndpoint + Send + DeserializeOwned> Client<'a, T> {
    /// Get the Ensembl response for the given identifier.
    /// Under the hood, this request will be bundled with other requests then returned asyncronously.
    /// # Panics
    ///
    /// Panics if the [Getter] has dropped, or the undelying channel has closed.
    pub async fn get(&self, id: String) -> Option<T> {
        let (tx, rx) = tokio::sync::oneshot::channel();
        if let Err(err) = self.tx.send((id, tx)).await {
            panic!(
                "Getter was closed or dropped recieving request: {}",
                err.0 .0
            )
        } // If the channel has closed, we can ignore the result
        rx.await.ok()
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct VEPAnalysis {
    pub input: String,
    pub strand: i8,
    pub assembly_name: String,
    pub seq_region_name: String,
    pub most_severe_consequence: String,
    pub start: u32,
    pub allele_string: String,
    pub transcript_consequences: Vec<TranscriptConsequence>,
}

#[derive(Debug, Serialize, Deserialize, Hash, PartialEq, Eq)]
pub struct TranscriptConsequence {
    pub transcript_id: String,
    pub impact: String,
    pub gene_id: String,
    pub gene_symbol: String,
    pub biotype: String,
    pub consequence_terms: Vec<String>,
    #[serde(default)]
    pub canonical: u8,
    pub nmd: Option<String>,
    #[serde(flatten)]
    pub protein_consequences: Option<ProteinConsequence>,
}
#[derive(Debug, Serialize, Deserialize, Hash, PartialEq, Eq)]
pub struct ProteinConsequence {
    pub hgvsp: String,
    pub hgvsc: String,
    pub protein_start: u32,
    pub protein_end: u32,
    pub codons: String,
    pub exon: String,
    pub amino_acids: String,
    pub cdna_start: u32,
    pub cdna_end: u32,
}

#[derive(Debug, Serialize, Deserialize, Hash, PartialEq, Eq)]
pub struct CodingSequence {
    pub query: String,
    pub id: String,
    pub desc: Option<String>,
    pub seq: String,
}

/// Data required to poll an Ensembl endpoint to create the given output.
pub trait EnsemblPostEndpoint {
    /// Return the URL extension for the Enseml endpoint.
    /// eg ```"/vep/human/hgvs"```
    fn extension() -> &'static str;
    /// Return a teplate for composing the json to be posted to the Ensembl endpoint.
    /// Should contain one insertion site for the list of identifiers requested.
    fn payload_template() -> &'static str;
    /// Get the input string from the Ensembl response. Will usually be &self.input
    fn input(&self) -> &str;
}
impl EnsemblPostEndpoint for VEPAnalysis {
    fn extension() -> &'static str {
        "/vep/human/hgvs"
    }
    fn payload_template() -> &'static str {
        r#"{"hgvs": 1, "numbers": 1, "canonical" : 1, "NMD" : 1, "hgvs_notations" : {ids}}"#
    }
    fn input(&self) -> &str {
        &self.input
    }
}

impl EnsemblPostEndpoint for CodingSequence {
    fn extension() -> &'static str {
        "/sequence/id"
    }
    fn payload_template() -> &'static str {
        r#"{"type": "cdna", "mask_feature" : 1, "ids" : {ids}}"#
    }
    fn input(&self) -> &str {
        &self.query
    }
}
