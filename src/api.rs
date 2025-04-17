use serde::de::DeserializeOwned;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use thiserror::Error;
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
///     let vep = h.await.unwrap();
///     println!("{:#?}", vep.unwrap());
/// }
/// # });
/// ```
#[cfg(not(target_arch = "wasm32"))]
#[derive(Debug)]
pub struct Getter<T: EnsemblPostEndpoint + Send + DeserializeOwned> {
    //is_alive: Arc<AtomicBool>,
    tx: mpsc::Sender<(
        String,
        tokio::sync::oneshot::Sender<Result<T, EnsemblError>>,
    )>,
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
        let (tx, mut rx) = mpsc::channel::<(
            String,
            tokio::sync::oneshot::Sender<Result<T, EnsemblError>>,
        )>(500);
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
        mut input: HashMap<String, tokio::sync::oneshot::Sender<Result<T, EnsemblError>>>,
        client: &reqwest::Client,
    ) {
        if input.is_empty() {
            return;
        }
        let mut map = HashMap::new();
        input.drain().into_iter().for_each(|(k, v)| {
            map.insert(k, v);
        });
        let mut ids_vec = map.keys().map(|s| s.clone()).collect::<Vec<_>>();

        while ids_vec.len() > 0 {
            let keys = ids_vec
                .drain(..usize::min(ids_vec.len(), T::max_post_size()))
                .collect::<Vec<_>>();
            let payload = T::payload_template().replace(r"{ids}", &json::stringify(keys.clone()));
            let response = client
                .post(String::from(ENSEMBL_SERVER) + T::extension())
                .header("Content-Type", "application/json")
                .header("Accept", "application/json")
                .body(payload)
                .send()
                .await
                .unwrap();
            match response.status().as_u16() {
                200 => {
                    let values = response.text().await.unwrap();
                    let outputs: Vec<T> = match serde_json::from_str(&values) {
                        Ok(outputs) => outputs,
                        Err(err) => {
                            if format!("{err}").as_str()
                                != "invalid type: map, expected a sequence at line 1 column 0"
                            {
                                eprintln!("{err}");
                            }
                            match serde_json::from_str::<HashMap<String, T>>(&values) {
                                Ok(outputs) => outputs.into_values().collect(),
                                Err(e) => {
                                    panic!(
                                        "Failed to parse the following response: {}\n{e:?}",
                                        values,
                                    );
                                }
                            }
                        }
                    };
                    for output in outputs.into_iter() {
                        let target = map.remove(output.input()).unwrap();
                        let _ = target.send(Ok(output));
                    }
                }
                400 => {
                    let error_message = response
                        .text()
                        .await
                        .unwrap_or("No detail included".to_string());
                    eprintln!("Bad Request: {error_message}");
                    for id in keys.into_iter() {
                        let _ = map.remove(&id).unwrap().send(Err(EnsemblError {
                            status_code: 400,
                            error: format!("Bad Request: {error_message}"),
                            input: id,
                        }));
                    }
                }
                403 => {
                    eprintln!(
                        "403 Forbidden: Too many requests. Waiting for 5 mins before trying again"
                    );
                    for id in keys.into_iter() {
                        let _ = map.remove(&id).unwrap().send(Err(EnsemblError {
                            status_code: 403,
                            error: format!("403 Forbidden: Too many requests."),
                            input: id,
                        }));
                    }
                    sleep(Duration::from_secs(300)).await;
                }
                404 => {
                    eprintln!("Not Found: Check your URL or request format.");
                    for id in keys.into_iter() {
                        let _ = map.remove(&id).unwrap().send(Err(EnsemblError {
                            status_code: 404,
                            error: "Not Found: Badly formatted request.".to_string(),
                            input: id,
                        }));
                    }
                }
                408 => {
                    eprintln!("Request Timeout. Pausing requests for 1 minute");
                    for id in keys.into_iter() {
                        let _ = map.remove(&id).unwrap().send(Err(EnsemblError {
                            status_code: 408,
                            error: "Request Timeout. Pausing requests for 1 minute".to_string(),
                            input: id,
                        }));
                    }
                    sleep(Duration::from_secs(60)).await;
                }
                429 => {
                    let reset_time = response
                        .headers()
                        .get("X-RateLimit-Reset")
                        .and_then(|v| v.to_str().ok())
                        .and_then(|v| v.parse::<u64>().ok())
                        .unwrap_or(60); // Default to 60 seconds if header is missing
                    eprintln!(
                        "Too Many Requests: Rate limit resets in {reset_time} seconds. Waiting..."
                    );
                    for id in keys.into_iter() {
                        let _ = map.remove(&id).unwrap().send(Err(EnsemblError {
                            status_code: 429,
                            error: format!(
                                "Too Many Requests: Rate limit resets in {reset_time} seconds."
                            ),
                            input: id,
                        }));
                    }
                    sleep(Duration::from_secs(reset_time)).await;
                }
                502 => {
                    eprintln!("Bad Gateway: Retrying after a pause...");
                    for id in keys.into_iter() {
                        let _ = map.remove(&id).unwrap().send(Err(EnsemblError {
                            status_code: 502,
                            error: "Bad gateway.".to_string(),
                            input: id,
                        }));
                    }
                    sleep(Duration::from_secs(10)).await;
                }
                503 => {
                    eprintln!("Service Unavailable: Retrying after a pause...");
                    for id in keys.into_iter() {
                        let _ = map.remove(&id).unwrap().send(Err(EnsemblError {
                            status_code: 503,
                            error: "Service Unavailable.".to_string(),
                            input: id,
                        }));
                    }
                    sleep(Duration::from_secs(10)).await;
                }
                status => {
                    let error_message = response
                        .text()
                        .await
                        .unwrap_or("No detail included".to_string());
                    eprintln!("Unexpected status code {status}: {}", &error_message);
                    for id in keys.into_iter() {
                        let _ = map.remove(&id).unwrap().send(Err(EnsemblError {
                            status_code: status as i16,
                            error: error_message.clone(),
                            input: id,
                        }));
                    }
                }
            }
        }
        for (id, sender) in &mut map.drain() {
            let _ = sender.send(Err(EnsemblError{ status_code: 404, error: format!("The input {id} did not give results, usually this means that it is not formated correctly.") , input: id,  }));
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
    tx: mpsc::Sender<(
        String,
        tokio::sync::oneshot::Sender<Result<T, EnsemblError>>,
    )>,
    getter: std::marker::PhantomData<&'a Getter<T>>,
}
impl<'a, T: 'static + EnsemblPostEndpoint + Send + DeserializeOwned + Clone> Client<'a, T> {
    /// Get the Ensembl response for the given identifier.
    /// Under the hood, this request will be bundled with other requests then returned asyncronously.
    /// # Panics
    ///
    /// Panics if the [Getter] has dropped, or the undelying channel has closed.
    pub async fn get(self, id: String) -> Result<T, EnsemblError> {
        let mut retries = 0;
        let (mut tx, mut rx) = tokio::sync::oneshot::channel();
        if let Err(err) = self.tx.send((id.clone(), tx)).await {
            panic!(
                "Getter was closed or dropped recieving request: {}",
                err.0 .0
            )
        }; // If the channel has closed, we can ignore the result
        loop {
            match rx.await {
                Ok(Ok(t)) => return Ok(t),
                Ok(Err(e)) => {
                    retries += 1;
                    if retries < 3 && [403, 408, 429, 502, 503].contains(&e.status_code) {
                        eprintln!("Error getting Ensembl data for {id}. Retry({retries})...\n{e}");
                        (tx, rx) = tokio::sync::oneshot::channel();
                        if let Err(err) = self.tx.send((id.clone(), tx)).await {
                            panic!(
                                "Getter was closed or dropped recieving request after {} retries: {}",
                                retries,
                                err.0 .0
                            )
                        };
                        continue;
                    } else {
                        return Err(e);
                    }
                }
                Err(e) => {
                    return Err(EnsemblError {
                        status_code: 0,
                        input: id,
                        error: e.to_string(),
                    })
                }
            }
        }
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

#[derive(Debug, Error, Clone, Serialize, Deserialize, Hash, PartialEq, Eq)]
#[error("Error accessing Ensembl data. Ensembl gave this response:\n{error}")]
pub struct EnsemblError {
    pub status_code: i16,
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
