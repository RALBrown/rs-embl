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
pub struct CodingSequence {
    pub query: String,
    pub id: String,
    pub desc: Option<String>,
    pub seq: String,
}

impl CodingSequence {
    /// ```
    /// use rs_embl::sequence::*;
    /// let test_seq = CodingSequence{
    /// query: "".to_owned(),
    /// id: "".to_owned(),
    /// desc: None,
    /// seq: "ACGTacgtACGTacgt".to_owned(),
    /// };
    /// assert_eq!(test_seq.exons(), vec!["ACGT".to_owned(),"ACGT".to_owned(),"ACGT".to_owned(), "ACGT".to_owned()]);
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
    /// let v = Getter::<CodingSequence>::new();
    /// let client = v.client();
    /// let handle = tokio::spawn(async move { client.get("ENST00000237014".to_owned()).await.expect("Is your internet connected?") });
    /// let cds = handle.await.unwrap();
    /// let exons = cds.exons();
    /// println!("{:#?}", exons);
    /// # });
    ///```
    pub fn exons(&self) -> Vec<String> {
        let mut output = Vec::new();
        let mut start = 0;
        let mut end = 1;
        while end < self.seq.len() {
            if char::from(self.seq.as_bytes()[start]).is_uppercase()
                == char::from(self.seq.as_bytes()[end]).is_uppercase()
            {
                end += 1;
                continue;
            } else {
                output.push(self.seq[start..end].to_uppercase());
                start = end;
                end += 1;
            }
        }
        if end > 1 {
            output.push(self.seq[start..].to_uppercase());
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
