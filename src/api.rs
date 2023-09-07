use serde::de::DeserializeOwned;
use std::collections::HashMap;
use tokio::sync::mpsc;
use tokio::time::{sleep, Duration};

/// The minimum time between post operations.
pub const WAIT_DELAY: Duration = Duration::from_millis(500);
const ENSEMBL_SERVER: &str = r#"https://rest.ensembl.org"#;
/// Encapsulates Ensembl REST API calls to allow multiple entries to be condensed into a single POST request.
///  * This will spawn a new asyncronous task that will periodically poll for new requests and handle them.
///  * The task will abort when the [Getter] object is dropped.
/// ```
/// # tokio::runtime::Builder::new_current_thread()
/// #       .enable_all()
/// #       .build()
/// #       .unwrap()
/// #       .block_on(async {
/// use rs_embl::{Getter, vep::VEPAnalysis};
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
    //is_alive: Arc<AtomicBool>,
    tx: mpsc::Sender<(String, tokio::sync::oneshot::Sender<T>)>,
}
// impl<T: EnsemblPostEndpoint + Send + DeserializeOwned> Drop for Getter<T> {
//     fn drop(&mut self) {
//         self.is_alive.store(false, Ordering::Relaxed);
//     }
// }
impl<T: 'static + EnsemblPostEndpoint + Send + DeserializeOwned> Default for Getter<T> {
    fn default() -> Self {
        Self::new()
    }
}
impl<T: 'static + EnsemblPostEndpoint + Send + DeserializeOwned> Getter<T> {
    /// Create a new Getter object to return T from Enseble REST endpoint.
    pub fn new() -> Self {
        let (tx, mut rx) = mpsc::channel::<(String, tokio::sync::oneshot::Sender<T>)>(500);
        //let is_alive = Arc::new(AtomicBool::new(true));
        let client = reqwest::Client::new();
        {
            //let is_alive = Arc::clone(&is_alive);
            tokio::spawn(async move {
                loop {
                    sleep(WAIT_DELAY).await;
                    let mut gets = HashMap::new();
                    let Some((key, value)) = rx.recv().await else { break };
                    gets.insert(key, value);
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
        Self { tx }
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