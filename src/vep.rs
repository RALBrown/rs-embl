//! Structures for the Variant Effect Predictor (VEP) endpoint of the Ensembl API.

use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VEPAnalysis {
    pub input: String,
    pub strand: i8,
    pub assembly_name: String,
    pub seq_region_name: String,
    pub most_severe_consequence: String,
    pub start: u32,
    pub allele_string: String,
    pub transcript_consequences: Vec<TranscriptConsequence>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Hash, PartialEq, Eq)]
pub struct TranscriptConsequence {
    pub transcript_id: String,
    pub impact: String,
    pub gene_id: String,
    pub gene_symbol: String,
    pub biotype: String,
    pub consequence_terms: Vec<String>,
    #[serde(default)]
    pub canonical: crate::Canonical,
    pub nmd: Option<String>,
    #[serde(flatten)]
    pub protein_consequences: Option<ProteinConsequence>,
}
#[derive(Debug, Clone, Serialize, Deserialize, Hash, PartialEq, Eq)]
pub struct ProteinConsequence {
    pub hgvsp: String,
    pub hgvsc: String,
    pub protein_start: u32,
    pub protein_end: u32,
    pub codons: String,
    pub exon: String,
    pub amino_acids: String,
    pub cdna_start: u32,
    pub cdna_end: u32,
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
}
