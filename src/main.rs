use rs_embl::{sequence::CdnaSequence, vep::VEPAnalysis, Getter};
#[cfg_attr(not(target_arch = "wasm32"), tokio::main)]
#[cfg_attr(target_arch = "wasm32", tokio::main(flavor = "current_thread"))]
async fn main() {
    let v = Getter::<VEPAnalysis>::new();
    let handles: Vec<_> = ["3:g.46373453_46373484del", "10:g.72346580_72346583dup"]
        .iter()
        .map(|id| {
            let v = v.client();
            tokio::spawn(async move { v.get(id.to_string()).await })
        })
        .collect();
    drop(v);
    let v2 = Getter::<CdnaSequence>::new();
    let handles2: Vec<_> = ["ENST00000237014", "ENSE00003556666"]
        .iter()
        .map(|id| {
            let v = v2.client();
            tokio::spawn(async move { v.get(id.to_string()).await })
        })
        .collect();
    drop(v2);
    for h in handles.into_iter() {
        let vep: Option<VEPAnalysis> = h.await.unwrap();
        println!("{:#?}", vep.unwrap());
    }
    for h in handles2.into_iter() {
        let vep = h.await.unwrap();
        println!("{:#?}", vep.unwrap());
    }
}
