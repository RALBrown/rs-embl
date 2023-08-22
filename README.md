<!-- cargo-sync-readme start -->

A tool for interacting with the POST endpoints of Ensembl REST API.
 * Spawns an async task that repeatedly polls requests made to its [Client] objects.
 * Bundles those requests and posts them to the Ensembl endpoints, asyncronously returning [serde::Deserialize] objects representing the result.
```rust
use rs_embl::{Getter, VEPAnalysis};

//Wrap in an async tokio runtime

let v = Getter::<VEPAnalysis>::new();
let handles: Vec<_> = ["3:g.46373453_46373484del", "10:g.72346580_72346583dup"]
    .iter()
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
