<!-- cargo-sync-readme start -->

A tool for interacting with the POST endpoints of Ensembl REST API.
 * Spawns an async task that repeatedly polls requests made to its [Client] objects.
 * Bundles those requests and posts them to the Ensembl endpoints, asyncronously returning [Deserialize] objects representing the result.
```rust
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
```

<!-- cargo-sync-readme end -->
