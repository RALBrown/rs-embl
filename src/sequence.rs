//! Structures for theSequence endpoint of the Ensembl API.
use serde::{Deserialize, Serialize};

#[derive(Debug, Serialize, Deserialize, Hash, PartialEq, Eq)]
pub struct CodingSequence {
    pub query: String,
    pub id: String,
    pub desc: Option<String>,
    pub seq: String,
}
impl CodingSequence {
    /**
    ```
    use rs_embl::sequence::*;
    let test_seq = CodingSequence{
        query: "".to_owned(),
        id: "".to_owned(),
        desc: None,
        seq: "ACGTacgtACGTacgt".to_owned(),
    };
    assert_eq!(test_seq.exons(), vec!["ACGT","acgt","ACGT", "acgt"]);
    ```
     */
    pub fn exons(&self) -> Vec<&str> {
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
                output.push(&self.seq[start..end]);
                start = end;
                end += 1;
            }
        }
        if end > 1 {
            output.push(&self.seq[start..]);
        }
        output
    }
}

impl crate::EnsemblPostEndpoint for CodingSequence {
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
