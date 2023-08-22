use rs_embl::{CodingSequence, Getter, VEPAnalysis};
#[tokio::main]
async fn main() {
    let v = Getter::<VEPAnalysis>::new();
    let handles: Vec<_> = ["3:g.46373453_46373484del", "10:g.72346580_72346583dup"]
        .iter()
        .map(|id| {
            let v = v.client();
            tokio::spawn(async move { v.get(id.to_string()).await })
        })
        .collect();

    let v2 = Getter::<CodingSequence>::new();
    let handles2: Vec<_> = ["ENST00000237014", "ENSE00003556666"]
        .iter()
        .map(|id| {
            let v = v2.client();
            tokio::spawn(async move { v.get(id.to_string()).await })
        })
        .collect();
    for h in handles.into_iter() {
        let vep: Option<VEPAnalysis> = h.await.unwrap();
        println!("{:#?}", vep.unwrap());
    }
    for h in handles2.into_iter() {
        let vep = h.await.unwrap();
        println!("{:#?}", vep.unwrap());
    }
}
