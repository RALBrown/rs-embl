use anyhow::Result;
use rs_embl::{transcript::Transcript, vep::VEPAnalysis, Getter};
#[tokio::main]
async fn main() -> Result<()> {
    let vep_getter = Getter::<VEPAnalysis>::new();
    let transcript_getter = Getter::<Transcript>::new();
    let handles = ["18:g.31592974G>A"]
        .iter()
        .map(|id| {
            let v = vep_getter.client();
            let t = transcript_getter.client();
            tokio::spawn(async move {
                let Some(vep) = v.get((*id).to_owned()).await else {return None};
                let consequences = vep.transcript_consequences.clone();
                let handles = consequences
                    .into_iter()
                    .map(|transcript| {
                        let client = t.clone();
                        tokio::spawn(
                            async move { client.get(transcript.transcript_id.clone()).await },
                        )
                    })
                    .collect::<Vec<_>>();
                Some((vep, handles))
            })
        })
        .collect::<Vec<_>>();
    for v in handles.into_iter() {
        let Some((vep, tcs)) = v.await.unwrap() else {continue};
        println!("{:#?}", vep);
        for tc in tcs {
            let Some(tc) = tc.await.unwrap() else {continue};
            println!("{:#?}", tc)
        }
    }

    Ok(())
}
