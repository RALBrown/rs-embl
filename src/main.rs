use rs_embl::{Getter, VEPAnalysis};
#[tokio::main]
async fn main() {
    let v = Getter::<VEPAnalysis>::new();
    let mut handles = Vec::new();
    for i in ["18:g.31592898delC", "18:g.31592998_31592999insC"] {
        let v = v.client();
        handles.push(tokio::spawn(async move { v.get(i.to_string()).await }));
    }
    for h in handles.into_iter() {
        let hgvs = h.await.unwrap();
        println!("{:#?}", hgvs.unwrap());
    }
}
