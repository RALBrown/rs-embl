use serde::{Deserialize, Serialize};

use crate::{
    sequence::{CdnaSequence, GenomicSequence},
    Client,
};

/**

*/
#[derive(Debug, Clone, Serialize, Deserialize, Hash, PartialEq, Eq)]
pub struct Transcript {
    pub id: String,
    pub display_name: String,
    pub start: u32,
    pub end: u32,
    pub strand: i8,
    #[serde(rename = "Translation")]
    pub translation: Option<Translation>,
    #[serde(rename = "UTR")]
    pub utrs: Vec<Utr>,
    #[serde(rename = "Exon")]
    pub exons: Vec<Exon>,
    #[serde(default, rename = "is_canonical")]
    pub canonical: crate::Canonical,
    pub species: String,
    #[serde(default)]
    biotype: crate::Biotype,
}
impl Transcript {
    pub async fn cdna_sequence(&self, client: Client<'static, CdnaSequence>) -> CdnaSequence {
        let id = self.id.clone();
        tokio::spawn(async move { client.get(id).await })
            .await
            .unwrap()
            .unwrap()
    }

    pub async fn genomic_sequence(
        &self,
        client: Client<'static, GenomicSequence>,
    ) -> GenomicSequence {
        let id = self.id.clone();
        tokio::spawn(async move { client.get(id).await })
            .await
            .unwrap()
            .unwrap()
    }
}

impl crate::EnsemblPostEndpoint for Transcript {
    fn extension() -> &'static str {
        "/lookup/id"
    }
    fn payload_template() -> &'static str {
        r#"{"expand": 1, "utr" : 1, "ids" : {ids}}"#
    }
    fn input(&self) -> &str {
        &self.id
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Hash, PartialEq, Eq)]
pub struct Utr {
    pub id: String,
    #[serde(rename = "Parent")]
    pub parent: String,
    pub start: u32,
    pub end: u32,
    pub strand: i8,
    #[serde(rename = "type")]
    pub utr_type: UtrType,
}

#[derive(Debug, Clone, Serialize, Deserialize, Hash, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum UtrType {
    FivePrimeUtr,
    ThreePrimeUtr,
}

#[derive(Debug, Clone, Serialize, Deserialize, Hash, PartialEq, Eq)]
pub struct Exon {
    pub id: String,
    pub start: u32,
    pub end: u32,
    pub strand: i8,
}

#[derive(Debug, Clone, Serialize, Deserialize, Hash, PartialEq, Eq)]
pub struct Translation {
    pub id: String,
    pub start: u32,
    pub end: u32,
    pub length: u32,
}
/// Provide the reverse complementary sequence of a nucleotide sequence
pub fn reverse_complement(seq: &str) -> String {
    let mut output = String::new();
    for b in seq.chars().rev() {
        output.push(match b {
            'a' => 't',
            'A' => 'T',
            'c' => 'g',
            'C' => 'G',
            'g' => 'c',
            'G' => 'C',
            't' => 'a',
            'T' => 'A',
            _ => panic!("{b} is not a recognized nucletide base"),
        });
    }
    output
}

pub fn translate(seq: &str) -> String {
    let mut output = String::new();
    for codon in seq
        .chars()
        .map(|c| {
            if c == 'T' {
                return 'U';
            } else {
                return c;
            }
        })
        .filter(|c| c.is_uppercase())
        .collect::<Vec<char>>()
        .chunks(3)
    {
        let aa = match codon {
            ['G', 'C', _] => 'A',
            ['U', 'G', 'U'] | ['U', 'G', 'C'] => 'C',
            ['G', 'A', 'U'] | ['G', 'A', 'C'] => 'D',
            ['G', 'A', 'A'] | ['G', 'A', 'G'] => 'E',
            ['U', 'U', 'U'] | ['U', 'U', 'C'] => 'F',
            ['G', 'G', _] => 'G',
            ['C', 'A', 'U'] | ['C', 'A', 'C'] => 'H',
            ['A', 'U', 'U'] | ['A', 'U', 'C'] | ['A', 'U', 'A'] => 'I',
            ['A', 'A', 'A'] | ['A', 'A', 'G'] => 'K',
            ['C', 'U', _] | ['U', 'U', 'A'] | ['U', 'U', 'G'] => 'L',
            ['A', 'U', 'G'] => 'M',
            ['A', 'A', 'U'] | ['A', 'A', 'C'] => 'N',
            ['C', 'C', _] => 'P',
            ['C', 'A', 'A'] | ['C', 'A', 'G'] => 'Q',
            ['C', 'G', _] | ['A', 'G', 'A'] | ['A', 'G', 'G'] => 'R',
            ['U', 'C', _] | ['A', 'G', 'U'] | ['A', 'G', 'C'] => 'S',
            ['A', 'C', _] => 'T',
            ['G', 'U', _] => 'V',
            ['U', 'G', 'G'] => 'W',
            ['U', 'A', 'U'] | ['U', 'A', 'C'] => 'Y',
            ['U', 'A', 'A'] | ['U', 'A', 'G'] | ['U', 'G', 'A'] => '*',
            _ => panic!("{codon:?} is not a recognized codon"),
        };
        output.push(aa);
        if aa == '*' {
            break;
        }
    }
    output
}

pub fn make_consequences(
    seq: &GenomicSequence,
    transcript: &Transcript,
    start: u32,
    end: u32,
    variant_allele: &str,
) -> Consequences {
    let mut edited_sequence: String = String::default();
    let upstream = &seq.seq[..(start - transcript.start) as usize];
    let downstream = &seq.seq[(end - transcript.start) as usize..];
    match (
        downstream.chars().next().unwrap().is_lowercase(),
        upstream.chars().last().unwrap().is_lowercase(),
    ) {
        (true, true) => return Consequences::Intron,
        (true, false) | (false, true) => {
            return Consequences::DisruptedSpliceSite;
        }
        (false, false) => {}
    }
    edited_sequence.push_str(upstream);
    edited_sequence.push_str(variant_allele);
    edited_sequence.push_str(downstream);

    let mut protein_sequence = String::default();
    if let Some(translation) = &transcript.translation {
        protein_sequence = if transcript.strand == 1 {
            translate(&edited_sequence[(translation.start - transcript.start) as usize..])
        } else {
            translate(&reverse_complement(
                &edited_sequence[(transcript.end - translation.end) as usize..],
            ))
        }
    }
    Consequences::Coding {
        genomic_sequence: edited_sequence,
        protein_sequence,
        nmd: true,
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Hash, PartialEq, Eq)]
pub enum Consequences {
    DisruptedSpliceSite,
    Coding {
        genomic_sequence: String,
        protein_sequence: String,
        nmd: bool,
    },
    Intron,
}
