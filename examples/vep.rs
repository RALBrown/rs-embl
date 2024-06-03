use anyhow::Result;
use rs_embl::{sequence::GenomicSequence, transcript::Transcript, vep::VEPAnalysis, Getter};
#[tokio::main]
async fn main() -> Result<()> {
    let vep_getter = Getter::<VEPAnalysis>::new();
    let transcript_getter = Getter::<Transcript>::new();
    let sequence_getter = Getter::<GenomicSequence>::new();
    let handles = ["18:g.31592974G>A"]
        .iter()
        .map(|id| {
            let v = vep_getter.client();
            let t = transcript_getter.client();
            tokio::spawn(async move {
                let Ok(vep) = v.get((*id).to_owned()).await else {
                    return None;
                };
                let consequences = vep.transcript_consequences.clone();
                let handles = consequences
                    .into_iter()
                    .map(|transcript| {
                        let client = t.clone();

                        match transcript {
                            rs_embl::vep::TranscriptConsequenceResponse::Parseable(t) => {
                                tokio::spawn(
                                    async move { client.get(t.transcript_id.clone()).await },
                                )
                            }
                            rs_embl::vep::TranscriptConsequenceResponse::Unparseable(_) => todo!(),
                        }
                    })
                    .collect::<Vec<_>>();
                Some((vep, handles))
            })
        })
        .collect::<Vec<_>>();
    for v in handles.into_iter() {
        let Some((vep, tcs)) = v.await.unwrap() else {
            continue;
        };
        println!("{:#?}", vep);
        for tc in tcs {
            let Ok(tc) = tc.await.unwrap() else {
                continue;
            };
            println!("{:#?}", tc);
            let genomic_seq = tc.genomic_sequence(sequence_getter.client()).await;
            let exons = genomic_seq.exons();
            println!("{:#?}", exons);
            let exon1 = exons.into_iter().next();
            let output = match exon1 {
                Some(s) => {
                    let Some(t) = tc.translation else { continue };
                    &(s.get((t.start - tc.start) as usize..).unwrap_or(""))
                }
                None => "",
            };
            println!("{}", output);
        }
    }

    Ok(())
}
