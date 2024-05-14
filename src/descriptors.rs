use serde::de::{self, Deserializer, Visitor};
use serde::{Deserialize, Serialize};
use std::fmt;
use std::str::FromStr;
use thiserror::Error;

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Eq, Ord, Hash, Serialize, Deserialize)]
#[serde(from = "i32", into = "i32")]
pub enum Canonical {
    CANONICAL,
    NONCANONICAL,
}
impl Default for Canonical {
    fn default() -> Self {
        Self::NONCANONICAL
    }
}
impl Into<i32> for Canonical {
    fn into(self) -> i32 {
        match &self {
            Canonical::CANONICAL => 1,
            Canonical::NONCANONICAL => 0,
        }
    }
}
impl From<i32> for Canonical {
    fn from(value: i32) -> Self {
        match value {
            0 => Self::NONCANONICAL,
            _ => Self::CANONICAL,
        }
    }
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize)]
pub enum Strand {
    PLUS,
    MINUS,
}
impl TryFrom<i32> for Strand {
    type Error = StrandError;
    fn try_from(value: i32) -> Result<Self, Self::Error> {
        match value {
            1 => Ok(Strand::PLUS),
            -1 => Ok(Strand::MINUS),
            _ => Err(StrandError::InvalidStrandNumber(value)),
        }
    }
}
impl FromStr for Strand {
    type Err = StrandError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "+" | "1" => Ok(Strand::PLUS),
            "-" | "-1" => Ok(Strand::MINUS),
            _ => Err(StrandError::InvalidStrandString(s.to_owned())),
        }
    }
}
impl TryFrom<&str> for Strand {
    type Error = StrandError;
    fn try_from(value: &str) -> Result<Self, Self::Error> {
        Self::from_str(value)
    }
}
impl std::fmt::Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Strand::PLUS => write!(f, "+"),
            Strand::MINUS => write!(f, "-"),
        }
    }
}
impl std::fmt::Debug for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self)
    }
}

#[derive(Error, Debug)]
pub enum StrandError {
    #[error("`{0}` is not a valid strand designator. Use '+' or '-'.")]
    InvalidStrandString(String),
    #[error("`{0}` is not a valid strand designator. Use '1' or '-1'.")]
    InvalidStrandNumber(i32),
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[allow(non_camel_case_types)]
pub enum Biotype {
    TR_V_gene,
    LRG_gene,
    miRNA,
    rRNA,
    tRNA,
    unprocessed_pseudogene,
    transcribed_pseudogene,
    transcribed_processed_pseudogene,
    IG_D_gene,
    TR_D_gene,
    lncRNA,
    IG_V_pseudogene,
    TR_J_gene,
    nonsense_mediated_decay,
    transcribed_unitary_pseudogene,
    sRNA,
    TR_J_pseudogene,
    RNase_MRP_RNA,
    telomerase_RNA,
    unitary_pseudogene,
    snRNA,
    Y_RNA,
    IG_C_pseudogene,
    cdna_update,
    IG_pseudogene,
    ribozyme,
    ncRNA_pseudogene,
    artifact,
    pseudogene,
    rRNA_pseudogene,
    Mt_rRNA,
    IG_J_pseudogene,
    non_stop_decay,
    aligned_transcript,
    TR_C_gene,
    mRNA,
    processed_pseudogene,
    TEC,
    scRNA,
    ccds_gene,
    scaRNA,
    vault_RNA,
    Mt_tRNA,
    other,
    misc_RNA,
    ncRNA,
    retained_intron,
    TR_V_pseudogene,
    protein_coding,
    IG_J_gene,
    snoRNA,
    RNase_P_RNA,
    translated_processed_pseudogene,
    antisense_RNA,
    transcribed_unprocessed_pseudogene,
    protein_coding_LoF,
    protein_coding_CDS_not_defined,
    processed_transcript,
    IG_V_gene,
    IG_C_gene,
    Unknown,
}
impl Default for Biotype {
    fn default() -> Self {
        Self::Unknown
    }
}

const FIELDS: &'static [&'static str] = &["+/1", "-/-1"];
impl<'de> Deserialize<'de> for Strand {
    fn deserialize<D>(deserializer: D) -> Result<Strand, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct StrandVisitor;

        impl<'de> Visitor<'de> for StrandVisitor {
            type Value = Strand;

            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                formatter.write_str("+ or - / 1 or -1")
            }

            fn visit_str<E>(self, value: &str) -> Result<Strand, E>
            where
                E: de::Error,
            {
                match value {
                    "+" | "1" => Ok(Strand::PLUS),
                    "-" | "-1" => Ok(Strand::MINUS),
                    _ => Err(de::Error::unknown_variant(value, FIELDS)),
                }
            }
            fn visit_i64<E>(self, value: i64) -> Result<Strand, E>
            where
                E: de::Error,
            {
                match value {
                    1 => Ok(Strand::PLUS),
                    -1 => Ok(Strand::MINUS),
                    _ => Err(de::Error::unknown_variant(
                        format!("{value}").as_str(),
                        FIELDS,
                    )),
                }
            }
        }

        deserializer.deserialize_any(StrandVisitor)
    }
}
