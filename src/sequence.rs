//! Structures for the Sequence endpoint of the Ensembl API.
use serde::{Deserialize, Serialize};

#[derive(Debug, Serialize, Deserialize, Hash, PartialEq, Eq)]
pub struct CdnaSequence {
    pub query: String,
    pub id: String,
    pub desc: Option<String>,
    pub seq: String,
}

#[derive(Debug, Serialize, Deserialize, Hash, PartialEq, Eq)]
pub struct GenomicSequence {
    pub query: String,
    pub id: String,
    pub desc: Option<String>,
    pub seq: String,
}

#[derive(Debug, Serialize, Deserialize, Hash, PartialEq, Eq)]
pub struct CodingSequence {
    pub query: String,
    pub id: String,
    pub desc: Option<String>,
    pub seq: String,
}

impl GenomicSequence {
    /// ```
    /// use rs_embl::sequence::*;
    /// let test_seq = GenomicSequence{
    /// query: "".to_owned(),
    /// id: "".to_owned(),
    /// desc: None,
    /// seq: "acACGTacgtACGTacgt".to_owned(),
    /// };
    /// assert_eq!(test_seq.exons(), vec!["ACGT","ACGT"]);
    /// ```
    /// ```
    /// # tokio::runtime::Builder::new_current_thread()
    /// #       .enable_all()
    /// #       .build()
    /// #       .unwrap()
    /// #       .block_on(async {
    /// use rs_embl::Getter;
    /// use rs_embl::sequence::*;
    /// // Wrap in an async tokio runtime
    ///
    /// // Instaniate a new Getter object,
    /// let v = Getter::<GenomicSequence>::new();
    /// let client = v.client();
    /// let handle = tokio::spawn(async move { client.get("ENST00000237014".to_owned()).await.expect("Is your internet connected?") });
    /// let sequence = handle.await.unwrap();
    /// let exons = sequence.exons();
    /// println!("{:#?}", exons);
    /// # });
    ///```
    pub fn exons(&self) -> Vec<&str> {
        let mut output = Vec::new();
        let mut start = 0;
        let mut end = 1;
        let mut is_upper = char::from(self.seq.as_bytes()[start]).is_uppercase();
        while end < self.seq.len() {
            if is_upper == char::from(self.seq.as_bytes()[end]).is_uppercase() {
                end += 1;
                continue;
            } else {
                if is_upper {
                    output.push(&self.seq[start..end]);
                }
                start = end;
                is_upper = char::from(self.seq.as_bytes()[start]).is_uppercase();
                end += 1;
            }
        }
        if end > 1 {
            if is_upper {
                output.push(&self.seq[start..]);
            }
        }
        output
    }
}

impl crate::EnsemblPostEndpoint for CodingSequence {
    fn extension() -> &'static str {
        "/sequence/id"
    }
    fn payload_template() -> &'static str {
        r#"{"type": "cds", "mask_feature" : 1, "ids" : {ids}}"#
    }
    fn input(&self) -> &str {
        &self.query
    }
}
impl crate::EnsemblPostEndpoint for CdnaSequence {
    fn extension() -> &'static str {
        "/sequence/id"
    }
    fn payload_template() -> &'static str {
        r#"{"type": "cdna", "mask_feature" : 1, "ids" : {ids}}"#
    }
    fn input(&self) -> &str {
        &self.query
    }
}
impl crate::EnsemblPostEndpoint for GenomicSequence {
    fn extension() -> &'static str {
        "/sequence/id"
    }
    fn payload_template() -> &'static str {
        r#"{"type": "genomic", "mask_feature" : 1, "ids" : {ids}}"#
    }
    fn input(&self) -> &str {
        &self.query
    }
}
