//! Structures for the Variant Effect Predictor (VEP) endpoint of the Ensembl API.

use std::str::FromStr;

use serde::{Deserialize, Serialize};
use thiserror::Error;

#[derive(Debug, Clone, Serialize, Deserialize, Hash, PartialEq, Eq)]
#[serde(untagged)]
pub enum VEPResult {
    Success(VEPAnalysis),
    EnsemblError(crate::api::EnsemblError),
    Error,
}
impl VEPResult {
    pub fn input(&self) -> &str {
        match self {
            VEPResult::Success(analysis) => &analysis.input,
            VEPResult::EnsemblError(error) => &error.input,
            VEPResult::Error => "ERROR",
        }
    }
}

/// A successful VEP analysis result.
#[derive(Debug, Clone, Serialize, Deserialize, Hash, PartialEq, Eq)]
pub struct VEPAnalysis {
    pub input: String,
    #[serde(default)]
    pub id: String,
    pub strand: i8,
    pub assembly_name: String,
    pub seq_region_name: String,
    pub most_severe_consequence: String,
    pub start: u32,
    pub end: u32,
    #[serde(rename = "allele_string")]
    pub allele: Allele,
    #[serde(default)]
    pub transcript_consequences: Vec<TranscriptConsequence>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Hash, PartialEq, Eq)]
pub struct TranscriptConsequence {
    pub transcript_id: String,
    pub impact: Option<String>,
    #[serde(default)]
    pub gene_id: String,
    #[serde(default)]
    pub gene_symbol: String,
    pub biotype: Option<String>,
    #[serde(default)]
    pub consequence_terms: Vec<String>,
    #[serde(default)]
    pub canonical: crate::Canonical,
    pub nmd: Option<String>,
    #[serde(flatten)]
    pub protein_consequences: Option<ProteinConsequence>,
    pub cdna_start: Option<u32>,
    pub cdna_end: Option<u32>,
    pub exon: Option<String>,
    pub intron: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Hash, PartialEq, Eq)]
pub struct ProteinConsequence {
    pub hgvsp: String,
    pub hgvsc: String,
    pub cds_start: u32,
    pub cds_end: u32,
    pub protein_start: u32,
    pub protein_end: u32,
    pub codons: String,
    pub amino_acids: String,
}

impl crate::EnsemblPostEndpoint for VEPAnalysis {
    fn extension() -> &'static str {
        "/vep/human/hgvs"
    }
    fn payload_template() -> &'static str {
        r#"{"hgvs": 1, "numbers": 1, "canonical" : 1, "NMD" : 1, "hgvs_notations" : {ids}}"#
    }
    fn input(&self) -> &str {
        &self.input
    }
    fn max_post_size() -> usize {
        200
    }
}

impl crate::EnsemblPostEndpoint for VEPResult {
    fn extension() -> &'static str {
        "/vep/human/hgvs"
    }
    fn payload_template() -> &'static str {
        r#"{"hgvs": 1, "numbers": 1, "canonical" : 1, "NMD" : 1, "hgvs_notations" : {ids}}"#
    }
    fn input(&self) -> &str {
        self.input()
    }
    fn max_post_size() -> usize {
        200
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Hash, PartialEq, Eq)]
#[serde(try_from = "String")]
pub struct Allele {
    pub normal: String,
    pub variant: String,
}
impl TryFrom<String> for Allele {
    type Error = AlleleParseError;
    fn try_from(value: String) -> Result<Self, Self::Error> {
        Self::from_str(&value)
    }
}
impl FromStr for Allele {
    type Err = AlleleParseError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let Some(mid) = s.find("/") else {
            return Err(AlleleParseError::NoSlash);
        };
        let mut normal = s[..mid].to_owned();
        let mut variant = s[mid + 1..].to_owned();
        if normal == "-" {
            normal = "".to_owned();
        };
        if variant == "-" {
            variant = "".to_owned();
        };
        Ok(Self { normal, variant })
    }
}

#[derive(Error, Debug)]
pub enum AlleleParseError {
    #[error("Allele strings need to conatain a /")]
    NoSlash,
}
