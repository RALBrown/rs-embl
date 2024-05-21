use serde::de::DeserializeOwned;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

use tokio::sync::mpsc;
use tokio::time::{sleep, Duration};

use tokio::spawn;

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
#[cfg(not(target_arch = "wasm32"))]
#[derive(Debug)]
pub struct Getter<T: EnsemblPostEndpoint + Send + DeserializeOwned> {
    //is_alive: Arc<AtomicBool>,
    tx: mpsc::Sender<(String, tokio::sync::oneshot::Sender<T>)>,
}

impl<T: 'static + EnsemblPostEndpoint + Send + DeserializeOwned> Default for Getter<T> {
    fn default() -> Self {
        Self::new()
    }
}
#[cfg(not(target_arch = "wasm32"))]
impl<T: 'static + EnsemblPostEndpoint + Send + DeserializeOwned> Getter<T> {
    /// Create a new Getter object to return T from Enseble REST endpoint.
    pub fn new() -> Self {
        let (tx, mut rx) = mpsc::channel::<(String, tokio::sync::oneshot::Sender<T>)>(500);
        {
            #[cfg(not(target_arch = "wasm32"))]
            let client = reqwest::Client::new();
            {
                spawn(async move {
                    loop {
                        sleep(WAIT_DELAY).await;
                        let mut gets = HashMap::new();
                        let Some((key, value)) = rx.recv().await else {
                            break;
                        };
                        gets.insert(key, value);
                        while let Ok((k, v)) = rx.try_recv() {
                            gets.insert(k, v);
                        }
                        #[cfg(not(target_arch = "wasm32"))]
                        Self::process(gets, &client).await;
                        #[cfg(target_arch = "wasm32")]
                        Self::process(gets).await;
                    }
                    rx.close();
                    let mut gets = HashMap::new();
                    while let Some((k, v)) = rx.recv().await {
                        gets.insert(k, v);
                    }
                    #[cfg(not(target_arch = "wasm32"))]
                    Self::process(gets, &client).await;
                    #[cfg(target_arch = "wasm32")]
                    Self::process(gets).await;
                });
            }
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
        let mut vals = input.drain();
        while vals.len() > 0 {
            let ids = (&mut vals).take(T::max_post_size()).collect::<Vec<_>>();
            let keys = ids.iter().map(|i| i.0.clone()).collect::<Vec<_>>();
            let mut map = HashMap::new();
            ids.into_iter().for_each(|(k, v)| {
                map.insert(k, v);
            });
            let payload = T::payload_template().replace(r"{ids}", &json::stringify(keys));
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
            let outputs: Vec<T> = match serde_json::from_str(&values) {
                Ok(outputs) => outputs,
                Err(err) => {
                    if format!("{err}").as_str()
                        != "invalid type: map, expected a sequence at line 1 column 0"
                    {
                        eprintln!("{err}");
                    }
                    if let Ok(e) = serde_json::from_str::<EnsemblTopLevelError>(&values) {
                        eprintln!("Ensembl Error: {}", e.error);
                        return;
                    }
                    match serde_json::from_str::<HashMap<String, T>>(&values) {
                        Ok(outputs) => outputs.into_values().collect(),
                        Err(e) => {
                            panic!("Failed to parse the following response: {}\n{e:?}", values,);
                        }
                    }
                }
            };
            for output in outputs.into_iter() {
                let target = map.remove(output.input()).unwrap();
                let _ = target.send(output); //if the sender's not listening that's its problem
            }
        }
    }

    #[cfg(target_arch = "wasm32")]
    async fn process(mut input: HashMap<String, tokio::sync::oneshot::Sender<T>>) {
        if input.is_empty() {
            return;
        }
        let ids: Vec<&str> = input.keys().map(|s| s.as_str()).collect();
        let payload = T::payload_template().replace(r"{ids}", &json::stringify(ids));
        let request = ehttp::Request {
            headers: ehttp::headers(&[
                ("Content-Type", "application/json"),
                ("Accept", "application/json"),
            ]),
            ..ehttp::Request::post(
                String::from(ENSEMBL_SERVER) + T::extension(),
                payload.into(),
            )
        };
        let (tx, resp) = tokio::sync::oneshot::channel();
        ehttp::fetch(request, move |result| {
            tx.send(result.unwrap().text().unwrap().to_owned());
        });
        let values = resp.await.unwrap();
        let outputs: Vec<T> = if let Ok(outputs) = serde_json::from_str(&values) {
            outputs
        } else {
            if let Ok(outputs) = serde_json::from_str::<HashMap<String, T>>(&values) {
                outputs.into_values().collect()
            } else {
                panic!("Failed to parse the following response: {}", values);
            }
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
impl<'a, T: 'static + EnsemblPostEndpoint + Send + DeserializeOwned> Client<'a, T> {
    /// Get the Ensembl response for the given identifier.
    /// Under the hood, this request will be bundled with other requests then returned asyncronously.
    /// # Panics
    ///
    /// Panics if the [Getter] has dropped, or the undelying channel has closed.
    pub async fn get(self, id: String) -> Option<T> {
        let (tx, rx) = tokio::sync::oneshot::channel();
        if let Err(err) = self.tx.send((id, tx)).await {
            panic!(
                "Getter was closed or dropped recieving request: {}",
                err.0 .0
            )
        }; // If the channel has closed, we can ignore the result
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
    // Get the maximum number of identifiers that can be sent in a single request.
    fn max_post_size() -> usize {
        50
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Hash, PartialEq, Eq)]
pub struct EnsemblTopLevelError {
    pub error: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Hash, PartialEq, Eq)]
pub struct EnsemblError {
    pub input: String,
    pub error: String,
}

#[cfg(target_arch = "wasm32")]
pub struct Getter<T: EnsemblPostEndpoint + DeserializeOwned> {
    tx: mpsc::Sender<(String, tokio::sync::oneshot::Sender<T>)>,
    rx: mpsc::Receiver<(String, tokio::sync::oneshot::Sender<T>)>,
    last_fetch: std::time::Instant,
    //to_fetch: HashMap<String, Sender<T>>,
}
#[cfg(target_arch = "wasm32")]
impl<T: 'static + EnsemblPostEndpoint + DeserializeOwned> Getter<T> {
    pub fn new() -> Self {
        let (tx, rx) = mpsc::channel(500);
        let last_fetch = std::time::Instant::now();
        Self { tx, rx, last_fetch }
    }

    pub async fn process(&mut self) {
        if (self.last_fetch.elapsed()) < WAIT_DELAY {
            return;
        }
        self.last_fetch = std::time::Instant::now();
        let mut input = HashMap::new();
        let Some((key, value)) = self.rx.recv().await else {
            return;
        };
        input.insert(key, value);
        while let Ok((k, v)) = self.rx.try_recv() {
            input.insert(k, v);
        }
        if input.is_empty() {
            return;
        }
        let ids: Vec<&str> = input.keys().map(|s| s.as_str()).collect();
        let payload = T::payload_template().replace(r"{ids}", &json::stringify(ids));
        let request = ehttp::Request {
            headers: ehttp::headers(&[
                ("Content-Type", "application/json"),
                ("Accept", "application/json"),
            ]),
            ..ehttp::Request::post(
                String::from(ENSEMBL_SERVER) + T::extension(),
                payload.into(),
            )
        };
        let (tx, resp) = tokio::sync::oneshot::channel();
        ehttp::fetch(request, move |result| {
            tx.send(result.unwrap().text().unwrap().to_owned());
        });
        let values = resp.await.unwrap();
        let outputs: Vec<T> = if let Ok(outputs) = serde_json::from_str(&values) {
            outputs
        } else {
            if let Ok(outputs) = serde_json::from_str::<HashMap<String, T>>(&values) {
                outputs.into_values().collect()
            } else {
                panic!("Failed to parse the following response: {}", values);
            }
        };
        for output in outputs.into_iter() {
            let target = input.remove(output.input()).unwrap();
            let _ = target.send(output); //if the sender's not listening that's it's problem
        }
    }
}
