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
