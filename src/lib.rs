//! A tool for interacting with the POST endpoints of Ensembl REST API.
//!  * Spawns an async task that repeatedly polls requests made to its [Client] objects.
//!  * Bundles those requests and posts them to the Ensembl endpoints, asyncronously returning [Deserialize] objects representing the result.
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
//! let mut handles: Vec<_> = ["18:g.31592898delC", "18:g.31592998_31592999insC"].iter()
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
const VEP_EXT: &str = "/vep/human/hgvs";
/** Encapsulates Ensembl REST API calls to allow multiple entries to be condensed into a single POST request.
 * This will spawn a new asyncronous task that will periodically poll for new requests and handle them.
 * The task will abort when the Getter object is dropped.
```
# tokio::runtime::Builder::new_current_thread()
#       .enable_all()
#       .build()
#       .unwrap()
#       .block_on(async {
use rs_embl::{Getter, VEPAnalysis};

//Wrap in an async tokio runtime

let v = Getter::<VEPAnalysis>::new();
let mut handles: Vec<_> = ["18:g.31592898delC", "18:g.31592998_31592999insC"].iter()
    .map(|id| {
        let v = v.client();
        tokio::spawn(async move { v.get(id.to_string()).await })
    })
    .collect();
for h in handles.into_iter() {
    let vep: Option<VEPAnalysis> = h.await.unwrap();
    println!("{:#?}", vep.unwrap());
}
# });
```
*/
#[derive(Debug)]
pub struct Getter<T: EnsemblEndpoint + Send + DeserializeOwned> {
    is_alive: Arc<AtomicBool>,
    tx: mpsc::Sender<(String, tokio::sync::oneshot::Sender<T>)>,
}
impl<T: EnsemblEndpoint + Send + DeserializeOwned> Drop for Getter<T> {
    fn drop(&mut self) {
        self.is_alive.store(false, Ordering::Relaxed);
    }
}
impl<T: 'static + EnsemblEndpoint + Send + DeserializeOwned> Getter<T> {
    ///Create a new Getter object to return T from Enseble REST endpoint.
    pub fn new() -> Self {
        let (tx, mut rx) = mpsc::channel::<(String, tokio::sync::oneshot::Sender<T>)>(500);
        let is_alive = Arc::new(AtomicBool::new(true));
        let client = reqwest::Client::new();
        {
            let is_alive = Arc::clone(&is_alive);
            let client = client.clone();
            tokio::spawn(async move {
                while is_alive.load(Ordering::Relaxed) {
                    sleep(WAIT_DELAY).await;
                    let mut gets = HashMap::new();
                    while let Ok((k, v)) = rx.try_recv() {
                        gets.insert(k, v);
                    }
                    if gets.is_empty() {
                        break;
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
        let ids: Vec<&str> = input.keys().map(|s| s.as_str()).collect();
        let payload = format!(
            r#"{{"hgvs": 1, "canonical" : 1, "NMD" : 1, "hgvs_notations" : {}}}"#,
            json::stringify(ids)
        );
        let values = client
            .post(String::from(ENSEMBL_SERVER) + VEP_EXT)
            .header("Content-Type", "application/json")
            .header("Accept", "application/json")
            .body(payload)
            .send()
            .await
            .unwrap()
            .text()
            .await
            .unwrap();
        println!("{}", values);
        let veps: Vec<T> = serde_json::from_str(&values).unwrap();
        for vep in veps.into_iter() {
            let target = input.remove(vep.input()).unwrap();
            target.send(vep);
        }
    }
    ///Create a trivially clonable Client that can be sent across async tasks.
    pub fn client(&self) -> Client<T> {
        Client::<T> {
            tx: self.tx.clone(),
        }
    }
}

/// A Client that can be cloned and sent across async tasks or threads to allow access to the underlying [Getter].
/// * Created by the [Getter::client()] method. [Client::clone()] is equivalent.
/// * Unlike the [Getter], Client implements [Send]. Thus, it is usually created in the parent task then passed to workers.
#[derive(Debug, Clone)]
pub struct Client<T> {
    tx: mpsc::Sender<(String, tokio::sync::oneshot::Sender<T>)>,
}
impl<T> Client<T> {
    /// Get the Ensembl response for the given identifier.
    /// Under the hood, this request will be bundled with other requests then returned asyncronously.
    pub async fn get(&self, id: String) -> Option<T> {
        let (tx, rx) = tokio::sync::oneshot::channel();
        self.tx.send((id, tx)).await;
        Some(rx.await.ok()?)
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct VEPAnalysis {
    pub input: String,
    strand: i8,
    assembly_name: String,
    seq_region_name: String,
    start: u32,
    allele_string: String,
    transcript_consequences: Vec<TranscriptConsequence>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct TranscriptConsequence {
    transcript_id: String,
    impact: String,
    gene_id: String,
    biotype: String,
    consequence_terms: Vec<String>,
    #[serde(default)]
    canonical: u8,
    nmd: Option<String>,
    hgvsp: Option<String>,
    hgvsc: Option<String>,
}
/// Data required to poll an Ensembl endpoint to create the given output.
pub trait EnsemblEndpoint {
    /// Return the URL extension for the Enseml endpoint.
    /// eg ```"/vep/human/hgvs"```
    fn extension() -> &'static str;
    /// Return a teplate for composing the json to be posted to the Ensembl endpoint.
    /// Should contain one insertion site for the list of identifiers requested.
    fn payload_template() -> &'static str;
    /// Get the input string from the Ensembl response. Will usually be &self.input
    fn input(&self) -> &str;
}
impl EnsemblEndpoint for VEPAnalysis {
    fn extension() -> &'static str {
        VEP_EXT
    }
    fn payload_template() -> &'static str {
        r#"{{"hgvs": 1, "canonical" : 1, "NMD" : 1, "hgvs_notations" : {}}}"#
    }
    fn input(&self) -> &str {
        &self.input
    }
}
