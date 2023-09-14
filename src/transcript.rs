use serde::{Deserialize, Serialize};

/**

*/
#[derive(Debug, Clone, Serialize, Deserialize, Hash, PartialEq, Eq)]
pub struct Transcript {
    pub id: String,
    pub display_name: String,
    pub start: u32,
    pub end: u32,
    pub strand: i8,
    #[serde(default, rename = "is_canonical")]
    pub canonical: crate::Canonical,
    pub species: String,
    #[serde(default)]
    biotype: crate::Biotype,
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
