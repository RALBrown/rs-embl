use itertools::Itertools;
use regex::Regex;
use serde::{Deserialize, Serialize};

use crate::{
    sequence::{CdnaSequence, GenomicSequence},
    Client,
};

const LAST_EJC_REGEX: &str = r".+([A-Z][a-z]+[A-Z]+)$";
const EMPTY_STR: &str = "";
/**

*/
#[derive(Debug, Clone, Serialize, Deserialize, Hash, PartialEq, Eq)]
pub struct Transcript {
    #[serde(default)]
    pub id: String,
    #[serde(default)]
    pub display_name: String,
    #[serde(default)]
    pub start: u32,
    #[serde(default)]
    pub end: u32,
    #[serde(default)]
    pub strand: i8,
    #[serde(rename = "Translation")]
    pub translation: Option<Translation>,
    #[serde(rename = "UTR", default)]
    pub utrs: Vec<Utr>,
    #[serde(rename = "Exon", default)]
    pub exons: Vec<Exon>,
    #[serde(default, rename = "is_canonical")]
    pub canonical: crate::Canonical,
    pub species: String,
    #[serde(default)]
    pub biotype: crate::Biotype,
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

    fn max_post_size() -> usize {
        1000
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

#[derive(Debug, Default, Clone, Serialize, Deserialize, Hash, PartialEq, Eq)]
pub struct TranslationConsequence {
    pub protein_sequence: String,
    pub stop_index: Option<usize>,
    pub last_ejc_index: Option<usize>,
    pub translation_type: TranslationType,
}

#[derive(Debug, Default, Clone, Serialize, Deserialize, Hash, PartialEq, Eq)]
pub enum TranslationType {
    NORMAL,
    NMD,
    NONSTOP,
    #[default]
    ERROR,
}

pub fn translate(seq: &str) -> TranslationConsequence {
    let last_ejc_capture = Regex::new(LAST_EJC_REGEX).unwrap().captures(seq);
    let last_ejc_index = match last_ejc_capture {
        Some(capture) => Some(capture.get(1).unwrap().start()),
        None => None,
    };
    let mut output = String::new();
    let mut counter: usize = 0;
    for codon in seq
        .chars()
        .map(|c| {
            counter += 1;
            if c == 'T' {
                return 'U';
            } else if c == 't' {
                return 'u';
            } else {
                return c;
            }
        })
        .filter(|c| c.is_uppercase())
        .tuples()
    {
        let aa = match codon {
            ('G', 'C', _) => 'A',
            ('U', 'G', 'U') | ('U', 'G', 'C') => 'C',
            ('G', 'A', 'U') | ('G', 'A', 'C') => 'D',
            ('G', 'A', 'A') | ('G', 'A', 'G') => 'E',
            ('U', 'U', 'U') | ('U', 'U', 'C') => 'F',
            ('G', 'G', _) => 'G',
            ('C', 'A', 'U') | ('C', 'A', 'C') => 'H',
            ('A', 'U', 'U') | ('A', 'U', 'C') | ('A', 'U', 'A') => 'I',
            ('A', 'A', 'A') | ('A', 'A', 'G') => 'K',
            ('C', 'U', _) | ('U', 'U', 'A') | ('U', 'U', 'G') => 'L',
            ('A', 'U', 'G') => 'M',
            ('A', 'A', 'U') | ('A', 'A', 'C') => 'N',
            ('C', 'C', _) => 'P',
            ('C', 'A', 'A') | ('C', 'A', 'G') => 'Q',
            ('C', 'G', _) | ('A', 'G', 'A') | ('A', 'G', 'G') => 'R',
            ('U', 'C', _) | ('A', 'G', 'U') | ('A', 'G', 'C') => 'S',
            ('A', 'C', _) => 'T',
            ('G', 'U', _) => 'V',
            ('U', 'G', 'G') => 'W',
            ('U', 'A', 'U') | ('U', 'A', 'C') => 'Y',
            ('U', 'A', 'A') | ('U', 'A', 'G') | ('U', 'G', 'A') => '*',
            _ => panic!("{codon:?} is not a recognized codon"),
        };
        output.push(aa);
        if aa == '*' {
            return TranslationConsequence {
                protein_sequence: output,
                stop_index: Some(counter),
                last_ejc_index,
                translation_type: match last_ejc_index {
                    None => TranslationType::NORMAL,
                    Some(last_ejc_index) => {
                        if counter + 50 < last_ejc_index {
                            TranslationType::NMD
                        } else {
                            TranslationType::NORMAL
                        }
                    }
                },
            };
        }
    }
    TranslationConsequence {
        protein_sequence: output,
        stop_index: None,
        last_ejc_index,
        translation_type: TranslationType::NONSTOP,
    }
}

pub fn make_consequences(
    seq: &GenomicSequence,
    transcript: &Transcript,
    start: u32,
    end: u32,
    variant_allele: &str,
) -> Consequences {
    let mut edited_sequence: String = String::default();
    let upstream;
    let downstream;
    if transcript.strand == 1 {
        upstream = &seq.seq[..usize::min((start - transcript.start) as usize, seq.seq.len())];
        downstream = seq
            .seq
            .get((end - transcript.start + 1) as usize..)
            .unwrap_or(EMPTY_STR);
    } else {
        upstream = &seq.seq[..usize::min((transcript.end - end) as usize, seq.seq.len())];
        downstream = seq
            .seq
            .get((transcript.end - start + 1) as usize..)
            .unwrap_or(EMPTY_STR);
    }
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
    if let Some(translation) = &transcript.translation {
        if transcript.strand == 1 {
            if start < translation.start {
                if end > translation.start {
                    return Consequences::LostStart;
                } else {
                    return Consequences::FivePrimeUTR;
                }
            }
        } else {
            if end > translation.end {
                if start < translation.end {
                    return Consequences::LostStart;
                } else {
                    return Consequences::FivePrimeUTR;
                }
            }
        }
    }
    edited_sequence.push_str(upstream);
    if transcript.strand == 1 {
        edited_sequence.push_str(variant_allele);
    } else {
        edited_sequence.push_str(&reverse_complement(variant_allele));
    }
    edited_sequence.push_str(downstream);

    let mut edited_protein_sequence = TranslationConsequence::default();
    let mut unedited_protein_sequence = TranslationConsequence::default();
    if let Some(translation) = &transcript.translation {
        edited_protein_sequence = translate(
            &edited_sequence[if transcript.strand == 1 {
                (translation.start - transcript.start) as usize
            } else {
                (transcript.end - translation.end) as usize
            }..],
        );
        unedited_protein_sequence = translate(
            &seq.seq[if transcript.strand == 1 {
                (translation.start - transcript.start) as usize
            } else {
                (transcript.end - translation.end) as usize
            }..],
        );
    }
    Consequences::Coding {
        edited_genomic_sequence: edited_sequence,
        edited_protein_sequence,
        unedited_protein_sequence,
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Hash, PartialEq, Eq)]
pub enum Consequences {
    DisruptedSpliceSite,
    LostStart,
    FivePrimeUTR,
    Coding {
        edited_genomic_sequence: String,
        edited_protein_sequence: TranslationConsequence,
        unedited_protein_sequence: TranslationConsequence,
    },
    Intron,
}

#[cfg(test)]
mod tests {
    use crate::sequence::GenomicSequence;

    const TTR_GENOME_SEQ: &str = "ACAGAAGTCCACTCATTCTTGGCAGGATGGCTTCTCATCGTCTGCTCCTCCTCTGCCTTGCTGGACTGGTATTTGTGTCTGAGGCTGGCCCTACGgtgagtgtttctgtgacatcccattcctacatttaagattcacgctaaatgaagtagaagtgactccttccagctttgccaaccagcttttattactagggcaagggtacccagcatctatttttaatataattaattcaaacttcaaaaagaatgaagttccactgagcttactgagctgggacttgaactctgagcattctacctcattgctttggtgcattaggtttgtaatatctggtacctctgtttcctcagatagatgatagaaataaagatatgatattaaggaagctgttaatactgaattttcagaaaagtatccctccataaaatgtatttgggggacaaactgcaggagattatattctggccctatagttattcaaaacgtatttattgattaatctttaaaaggcttagtgaacaatattctagtcagatatctaattcttaaatcctctagaagaattaactaatactataaaatgggtctggatgtagttctgacattattttataacaactggtaagagggagtgactatagcaacaactaaaatgatctcaggaaaacctgtttggccctatgtatggtacattacatcttttcagtaattccactcaaatggagacttttaacaaagcaactgttctcaggggacctattttctcccttaaaattcattatacacatccctggttgatagcagtgtgtctggaggcagaaaccattcttgctttggaaacaattacgtctgtgttatactgagtagggaagctcattaattgtcgacacttacgttcctgataatgggatcagtgtgtaattcttgtttcgctccagatttctaataccacaaagaataaatcctttcactctgatcaattttgttaacttctcacgtgtcttctctacacccagGGCACCGGTGAATCCAAGTGTCCTCTGATGGTCAAAGTTCTAGATGCTGTCCGAGGCAGTCCTGCCATCAATGTGGCCGTGCATGTGTTCAGAAAGGCTGCTGATGACACCTGGGAGCCATTTGCCTCTGGgtaagttgccaaagaaccctcccacaggacttggttttatcttcccgtttgcccctcacttggtagagagaggctcacatcatctgctaaagaatttacaagtagattgaaaaacgtaggcagaggtcaagtatgccctctgaaggatgccctctttttgttttgcttagctaggaagtgaccaggaacctgagcatcatttaggggcagacagtagagaaaagaaggaatcagaactcctctcctctagctgtggtttgcaacccttttgggtcacagaacactttatgtaggtgatgaaaagtaaacattctatgcccagaaaaaatgcacagatacacacacatacaaaatcatatatgtgattttaggagtttcacagattccctggtgtccctgggtaacaccaaagctaagtgtccttgtcttagaattttaggaaaaggtataatgtgtattaacccattaacaaaaggaaaggaattcagaaatattattaaccaggcatctgtctgtagttaatatggatcacccaaaacccaaggcttttgcctaatgaacactttggggcacctactgtgtgcaaggctgggggctgtcaagctcagttaaaaaaaaaaagatagaagagatggatccatgaggcaaagtacagccccaggctaatcccacgatcacccgacttcatgtccaagagtggcttctcaccttcattagccagttcacaattttcatggagtttttctacctgcactagcaaaaacttcaaggaaaatacatattaataaatctaagcaaagtgaccagaagacagagcaatcaggagaccctttgcatccagcagaagaggaactgctaagtatttacatctccacagagaagaatttctgttgggttttaattgaaccccaagaaccacatgattcttcaaccattattgggaagatcattttcttaggtctggttttaactggctttttatttgggaattcatttatgtttatataaaatgccaagcataacatgaaaagtggttacaggactattctaagggagagacagaatggacaccaaaaatattccaatgttcttgtgaatcttttccttgcaccaggacaaaaaaaaaaagaagtgaaaagaagaaaggaggaggggcataatcagagtcagtaaagacaactgctatttttatctatcgtagctgttgcagtcaaatgggaagcaatttccaacattcaactatggagctggtacttacatggaaatagaagttgcctagtgtttgttgctggcaaagagttatcagagaggttaaatatataaaagggaaaagagtcagatacaggttcttcttcctactttaggttttccactgtgtgtgcaaatgatactccctggtggtgtgcagatgcctcaaagctatcctcacaccacaagggagaggagcgagatcctgctgtcctggagaagtgcagagttagaacagctgtggccacttgcatccaatcatcaatcttgaatcacagggactctttcttaagtaaacattatacctggccgggcacggtggctcacgcctgtaatcccagcactttgggatgccaaagtgggcatatcatctgaggtcaggagttcaagaccagcctggccaacatggcaaaactccgtctttatgaaaaatacaaaaattagccaggcatggtggcaggcgcctgtaatcccagctaattgggaggctgaggctggagaatcccttgaatctaggaggcagaggttgcagtgagctgagatcgtgccattgcactccagcctgggtgacaagagtaaaactctgtctcaaaaaaaaaaaattatacctacattctcttcttatcagagaaaaaaatctacagtgagcttttcaaaaagtttttacaaactttttgccatttaatttcagttaggagttttccctacttctgacttagttgaggggaaatgttcataacatgtttataacatgtttatgtgtgttagttggtgggggtgtattactttgccatgccatttgtttcctccatgcgtaacttaatccagactttcacaccttatagGAAAACCAGTGAGTCTGGAGAGCTGCATGGGCTCACAACTGAGGAGGAATTTGTAGAAGGGATATACAAAGTGGAAATAGACACCAAATCTTACTGGAAGGCACTTGGCATCTCCCCATTCCATGAGCATGCAGAGgtgagtatacagaccttcgagggttgttttggttttggtttttgcttttggcattccaggaaatgcacagttttactcagtgtaccacagaaatgtcctaaggaaggtgatgaatgaccaaaggttccctttcctattatacaagaaaaaattcacaacactctgagaagcaaatttctttttgactttgatgaaaatccacttagtaacatgacttgaacttacatgaaactactcatagtctattcattccactttatatgaatattgatgtatctgctgttgaaataatagtttatgaggcagccctccagaccccacgtagagtgtatgtaacaagagatgcaccattttatttctcgaaaacccgtaacattcttcattccaaaacacatctggcttctcggaggtctggacaagtgattcttggcaacacatacctatagagacaataaaatcaaagtaataatggcaacacaatagataacatttaccaagcatacaccatgtggcagacacaattataagtgttttccatatttaacctacttaatcctcaggaataagccactgaggtcagtcctattattatccccatcttatagatgaagaaaatgaggcaccaggaagtcaaataacttgtcaaaggtcacaagactaggaaatacacaagtagaaatgtttacaattaaggcccaggctgggtttgccctcagttctgctatgcctcgcattatgccccaggaaactttttcccttgtgaaagccaagcttaaaaaaagaaaagccacatttgtaacgtgctctgttcccctgcctatggtgaggatcttcaaacagttatacatggacccagtccccctgccttctccttaatttcttaagtcatttgaaacagatggctgtcatggaaatagaatccagacatgttggtcagagttaaagatcaactaattccatcaaaaatagctcggcatgaaagggaactattctctggcttagtcatggatgagactttcaattgctataaagtggttcctttattagacaatgttaccagggaaacaacaggggtttgtttgacttctggggcccacaagtcaacaagagagccccatctaccaaggagcatgtccctgactacccctcagccagcagcaagacatggaccccagtcagggcaggagcagggtttcggcggcgcccagcacaagacattgcccctagagtctcagcccctaccctcgagtaatagatctgcctacctgagactgttgtttgcccaagagctgggtctcagcctgatgggaaccatataaaaaggttcactgacatactgcccacatgttgttctctttcattagatcttagcttccttgtctgctcttcattcttgcagtattcattcaacaaacattaaaaaaaaaaaaaagcattctatgtgtggaacactctgctagatgctgtggatttagaaatgaaaatacatcccgacccttggaatggaagggaaaggactgaagtaagacagattaagcaggaccgtcagcccagcttgaagcccagataaatacggagaacaagagagagcgagtagtgagagatgagtcccaatgcctcactttggtgacgggtgcgtggtgggcttcatgcagcttcttctgataaatgcctccttcagaactggtcaactctaccttggccagtgacccaggtggtcatagtagatttaccaagggaaaatggaaacttttattaggagctcttaggcctcttcacttcatggatttttttttcctttttttttgagatggagttttgccctgtcacccaggctggaatgcagtggtgcaatctcagctcactgcaacctccgcctcccaggttcaagcaattctcctgcctcagcctcccgagtagctgggactacaggtgtgcgccaccacaccaggctaatttttgtattttttgtaaagacaggttttcaccacgttggccaggctggtctgaactccagacctcaggtgattcacctgtctcagcctcccaaagtgctgggattacaggtgtgagccaccgtgcccggctacttcatggatttttgattacagattatgcctcttacaatttttaagaagaatcaagtgggctgaaggtcaatgtcaccataagacaaaagacatttttattagttgattctagggaattggccttaaggggagccctttcttcctaagagattcttaggtgattctcacttcctcttgccccagtattatttttgtttttggtatggctcactcagatccttttttcctcctatccctaagtaatccgggtttctttttcccatatttagaacaaaatgtatttatgcagagtgtgtccaaacctcaacccaaggcctgtatacaaaataaatcaaattaaacacatctttactgtcttctacctctttcctgacctcaatatatcccaacttgcctcactctgagaaccaaggctgtcccagcacctgagtcgcagatattctactgatttgacagaactgtgtgactatctggaacagcattttgatccacaatttgcccagttacaaagcttaaatgagctctagtgcatgcatatatatttcaaaattccaccatgatcttccacactctgtattgtaaatagagccctgtaatgcttttacttcgtatttcattgcttgttatacataaaaatatacttttcttcttcatgttagaaaatgcaaagaataggagggtgggggaatctctgggcttggagacaggagacttgccttcctactatggttccatcagaatgtagactgggacaatacaataattcaagtctggtttgctcatctgtaaattgggaagaatgtttccagctccagaatgctaaatctctaagtctgtggttggcagccactattgcagcagctcttcaatgactcaatgcagttttgcattctccctaccttttttttctaaaaccaataaaatagatacagcctttaggctttctgggatttcccttagtcaagctagggtcatcctgactttcggcgtgaatttgcaaaacaagacctgactctgtactcctgctctaaggactgtgcatggttccaaaggcttagcttgccagcatatttgagctttttccttctgttcaaactgttccaaaatataaaagaataaaattaattaagttggcactggacttccggtggtcagtcatgtgtgtcatctgtcacgtttttcgggctctggtggaaatggatctgtctgtcttctctcatagGTGGTATTCACAGCCAACGACTCCGGCCCCCGCCGCTACACCATTGCCGCCCTGCTGAGCCCCTACTCCTATTCCACCACGGCTGTCGTCACCAATCCCAAGGAATGAGGGACTTCTCCTCCAGTGGACCTGAAGGACGAGGGATGGGATTTCATGTAACCAAGAGTATTCCATTTTTACTAAAGCAGTGTTTTCACCTCATATGCTATGTTAGAAGTCCAGGCAGAGACAATAAAACATTCCTGTGAAAGGCA";
    const TTR_201_JSON: &str = r#"{"end":31598821,"object_type":"Transcript","is_canonical":1,"length":616,"db_type":"core","id":"ENST00000237014","Translation":{"version":4,"species":"homo_sapiens","start":31591903,"length":147,"id":"ENSP00000237014","db_type":"core","Parent":"ENST00000237014","end":31598675,"object_type":"Translation"},"species":"homo_sapiens","display_name":"TTR-201","start":31591877,"version":8,"seq_region_name":"18","assembly_name":"GRCh38","logic_name":"ensembl_havana_transcript_homo_sapiens","Exon":[{"species":"homo_sapiens","start":31591877,"version":2,"assembly_name":"GRCh38","seq_region_name":"18","end":31591971,"object_type":"Exon","db_type":"core","id":"ENSE00001836564","strand":1},{"start":31592896,"species":"homo_sapiens","seq_region_name":"18","assembly_name":"GRCh38","version":1,"end":31593026,"object_type":"Exon","id":"ENSE00003556666","db_type":"core","strand":1},{"id":"ENSE00000796939","db_type":"core","strand":1,"end":31595255,"object_type":"Exon","version":1,"seq_region_name":"18","assembly_name":"GRCh38","species":"homo_sapiens","start":31595120},{"end":31598821,"object_type":"Exon","db_type":"core","id":"ENSE00001827041","strand":1,"start":31598568,"species":"homo_sapiens","seq_region_name":"18","assembly_name":"GRCh38","version":2}],"strand":1,"Parent":"ENSG00000118271","source":"ensembl_havana","UTR":[{"assembly_name":"GRCh38","seq_region_name":"18","start":31591877,"source":"ensembl_havana","type":"five_prime_utr","species":"homo_sapiens","db_type":"core","id":"ENST00000237014","strand":1,"Parent":"ENST00000237014","end":31591902,"object_type":"five_prime_UTR"},{"type":"three_prime_utr","species":"homo_sapiens","source":"ensembl_havana","start":31598676,"seq_region_name":"18","assembly_name":"GRCh38","object_type":"three_prime_UTR","end":31598821,"strand":1,"Parent":"ENST00000237014","id":"ENST00000237014","db_type":"core"}],"biotype":"protein_coding"}"#;
    const TTR_V30M_DEL_PROTEIN: &str =
        "MASHRLLLLCLAGLVFVSEAGPTGTGESKCPLMVKVLDAVRGSPAINVACMCSERLLMTPGSHLPLGKPVSLESCMGSQLRRNL*";
    const TTR_V30M_INS_PROTEIN: &str = "MASHRLLLLCLAGLVFVSEAGPTGTGESKCPLMVKVLDAVRGSPAINVAGACVQKGC*";
    #[test]
    fn test_snp() {
        let transcript = serde_json::from_str::<super::Transcript>(TTR_201_JSON).unwrap();
        let genomic_seq: GenomicSequence = GenomicSequence {
            query: "".to_owned(),
            id: "".to_owned(),
            desc: None,
            seq: TTR_GENOME_SEQ.to_owned(),
        };
        const END: u32 = 31592974;
        const START: u32 = 31592974;
        const VARIANT_ALLELE: &str = "A";
        let consequences =
            super::make_consequences(&genomic_seq, &transcript, START, END, VARIANT_ALLELE);
        let super::Consequences::Coding {
            edited_genomic_sequence,
            edited_protein_sequence,
            unedited_protein_sequence,
        } = consequences
        else {
            panic!()
        };
        const V30M_TTR: &str = "MASHRLLLLCLAGLVFVSEAGPTGTGESKCPLMVKVLDAVRGSPAINVAMHVFRKAADDTWEPFASGKTSESGELHGLTTEEEFVEGIYKVEIDTKSYWKALGISPFHEHAEVVFTANDSGPRRYTIAALLSPYSYSTTAVVTNPKE*";
        assert_eq!(&edited_protein_sequence.protein_sequence, V30M_TTR);
    }
    #[test]
    fn test_del() {
        let transcript = serde_json::from_str::<super::Transcript>(TTR_201_JSON).unwrap();
        let genomic_seq: GenomicSequence = GenomicSequence {
            query: "".to_owned(),
            id: "".to_owned(),
            desc: None,
            seq: TTR_GENOME_SEQ.to_owned(),
        };
        const END: u32 = 31592974;
        const START: u32 = 31592974;
        const VARIANT_ALLELE: &str = "-";
        let consequences =
            super::make_consequences(&genomic_seq, &transcript, START, END, VARIANT_ALLELE);
        let super::Consequences::Coding {
            edited_genomic_sequence,
            edited_protein_sequence,
            unedited_protein_sequence,
        } = consequences
        else {
            panic!()
        };
        assert_eq!(
            &edited_protein_sequence.protein_sequence,
            TTR_V30M_DEL_PROTEIN
        );
    }
    #[test]
    fn test_ins() {
        let transcript = serde_json::from_str::<super::Transcript>(TTR_201_JSON).unwrap();
        let genomic_seq: GenomicSequence = GenomicSequence {
            query: "".to_owned(),
            id: "".to_owned(),
            desc: None,
            seq: TTR_GENOME_SEQ.to_owned(),
        };
        const END: u32 = 31592974;
        const START: u32 = 31592975;
        const VARIANT_ALLELE: &str = "G";
        let consequences =
            super::make_consequences(&genomic_seq, &transcript, START, END, VARIANT_ALLELE);
        let super::Consequences::Coding {
            edited_genomic_sequence,
            edited_protein_sequence,
            unedited_protein_sequence,
        } = consequences
        else {
            panic!()
        };
        assert_eq!(
            &edited_protein_sequence.protein_sequence,
            TTR_V30M_INS_PROTEIN
        );
    }
    const JSON: &str = {
        r##"{"ENST00000457901":{"source":"havana","logic_name":"havana_homo_sapiens","seq_region_name":"2","id":"ENST00000457901","version":1,"strand":1,"start":21221185,"assembly_name":"GRCh38","UTR":[],"db_type":"core","object_type":"Transcript","length":504,"biotype":"lncRNA","is_canonical":0,"Exon":[{"seq_region_name":"2","version":1,"id":"ENSE00001774670","start":21221185,"strand":1,"assembly_name":"GRCh38","end":21221294,"species":"homo_sapiens","db_type":"core","object_type":"Exon"},{"seq_region_name":"2","version":1,"id":"ENSE00001710215","strand":1,"start":21263693,"assembly_name":"GRCh38","db_type":"core","object_type":"Exon","end":21264086,"species":"homo_sapiens"}],"end":21264086,"species":"homo_sapiens","Parent":"ENSG00000233005"},"ENST00000368926":{"source":"ensembl_havana","logic_name":"ensembl_havana_transcript_homo_sapiens","Translation":{"species":"homo_sapiens","end":151050458,"db_type":"core","object_type":"Translation","Parent":"ENST00000368926","length":341,"id":"ENSP00000357922","version":5,"start":151047848},"seq_region_name":"1","id":"ENST00000368926","version":6,"start":151047751,"strand":1,"assembly_name":"GRCh38","UTR":[{"assembly_name":"GRCh38","end":151047847,"species":"homo_sapiens","Parent":"ENST00000368926","db_type":"core","object_type":"five_prime_UTR","type":"five_prime_utr","source":"ensembl_havana","seq_region_name":"1","id":"ENST00000368926","strand":1,"start":151047751},{"seq_region_name":"1","id":"ENST00000368926","source":"ensembl_havana","type":"three_prime_utr","strand":1,"start":151050459,"assembly_name":"GRCh38","Parent":"ENST00000368926","object_type":"three_prime_UTR","db_type":"core","end":151051420,"species":"homo_sapiens"}],"db_type":"core","object_type":"Transcript","length":2085,"biotype":"protein_coding","is_canonical":1,"Exon":[{"species":"homo_sapiens","end":151048852,"object_type":"Exon","db_type":"core","assembly_name":"GRCh38","start":151047751,"strand":1,"version":6,"id":"ENSE00001448297","seq_region_name":"1"},{"assembly_name":"GRCh38","db_type":"core","object_type":"Exon","species":"homo_sapiens","end":151051420,"id":"ENSE00001712848","version":2,"seq_region_name":"1","start":151050438,"strand":1}],"end":151051420,"species":"homo_sapiens","Parent":"ENSG00000143443","display_name":"C1orf56-201"},"ENST00000491825":{"source":"havana","logic_name":"havana_homo_sapiens","seq_region_name":"1","version":1,"id":"ENST00000491825","start":151055583,"strand":-1,"assembly_name":"GRCh38","UTR":[],"object_type":"Transcript","db_type":"core","length":837,"biotype":"protein_coding_CDS_not_defined","is_canonical":0,"Exon":[{"strand":-1,"start":151059479,"seq_region_name":"1","version":1,"id":"ENSE00001871667","end":151059773,"species":"homo_sapiens","object_type":"Exon","db_type":"core","assembly_name":"GRCh38"},{"strand":-1,"start":151056656,"seq_region_name":"1","version":1,"id":"ENSE00001809953","end":151056786,"species":"homo_sapiens","object_type":"Exon","db_type":"core","assembly_name":"GRCh38"},{"start":151055583,"strand":-1,"version":1,"id":"ENSE00001943806","seq_region_name":"1","species":"homo_sapiens","end":151055993,"db_type":"core","object_type":"Exon","assembly_name":"GRCh38"}],"end":151059773,"species":"homo_sapiens","Parent":"ENSG00000197622","display_name":"CDC42SE1-205"},"ENST00000670105":{"start":21221169,"strand":1,"version":1,"id":"ENST00000670105","seq_region_name":"2","logic_name":"havana_tagene_homo_sapiens","source":"havana_tagene","object_type":"Transcript","db_type":"core","UTR":[],"assembly_name":"GRCh38","biotype":"lncRNA","length":642,"Parent":"ENSG00000233005","species":"homo_sapiens","end":21264078,"Exon":[{"assembly_name":"GRCh38","object_type":"Exon","db_type":"core","end":21221294,"species":"homo_sapiens","seq_region_name":"2","id":"ENSE00003869071","version":1,"start":21221169,"strand":1},{"assembly_name":"GRCh38","object_type":"Exon","db_type":"core","species":"homo_sapiens","end":21258899,"id":"ENSE00003853909","version":1,"seq_region_name":"2","start":21258770,"strand":1},{"seq_region_name":"2","version":1,"id":"ENSE00003878245","start":21263693,"strand":1,"assembly_name":"GRCh38","db_type":"core","object_type":"Exon","end":21264078,"species":"homo_sapiens"}],"is_canonical":0},"ENST00000622592":{"biotype":"lncRNA","length":859,"Exon":[{"object_type":"Exon","db_type":"core","end":10342669,"species":"homo_sapiens","assembly_name":"GRCh38","strand":-1,"start":10342616,"seq_region_name":"21","id":"ENSE00003717148","version":1},{"species":"homo_sapiens","end":10340478,"db_type":"core","object_type":"Exon","assembly_name":"GRCh38","strand":-1,"start":10340414,"version":1,"id":"ENSE00003714071","seq_region_name":"21"},{"strand":-1,"start":10338455,"id":"ENSE00003720003","version":1,"seq_region_name":"21","db_type":"core","object_type":"Exon","species":"homo_sapiens","end":10338566,"assembly_name":"GRCh38"},{"species":"homo_sapiens","end":10329038,"object_type":"Exon","db_type":"core","assembly_name":"GRCh38","strand":-1,"start":10328411,"version":1,"id":"ENSE00003716827","seq_region_name":"21"}],"is_canonical":1,"Parent":"ENSG00000277693","species":"homo_sapiens","end":10342669,"version":1,"id":"ENST00000622592","seq_region_name":"21","logic_name":"havana_homo_sapiens","source":"havana","start":10328411,"strand":-1,"assembly_name":"GRCh38","db_type":"core","object_type":"Transcript","UTR":[]},"ENST00000465135":{"assembly_name":"GRCh38","db_type":"core","object_type":"Transcript","UTR":[],"seq_region_name":"1","logic_name":"havana_homo_sapiens","version":1,"id":"ENST00000465135","source":"havana","strand":1,"start":151048569,"Exon":[{"strand":1,"start":151048569,"seq_region_name":"1","id":"ENSE00001827584","version":1,"object_type":"Exon","db_type":"core","end":151048852,"species":"homo_sapiens","assembly_name":"GRCh38"},{"seq_region_name":"1","version":1,"id":"ENSE00001847261","start":151050438,"strand":1,"assembly_name":"GRCh38","object_type":"Exon","db_type":"core","end":151050484,"species":"homo_sapiens"},{"id":"ENSE00001829671","version":1,"seq_region_name":"1","start":151051885,"strand":1,"assembly_name":"GRCh38","db_type":"core","object_type":"Exon","species":"homo_sapiens","end":151051986}],"is_canonical":0,"Parent":"ENSG00000143443","display_name":"C1orf56-202","end":151051986,"species":"homo_sapiens","biotype":"protein_coding_CDS_not_defined","length":433},"ENST00000470278":{"Exon":[{"end":151059773,"species":"homo_sapiens","db_type":"core","object_type":"Exon","assembly_name":"GRCh38","strand":-1,"start":151059479,"seq_region_name":"1","version":1,"id":"ENSE00001871667"},{"start":151056656,"strand":-1,"seq_region_name":"1","id":"ENSE00001809953","version":1,"end":151056786,"species":"homo_sapiens","object_type":"Exon","db_type":"core","assembly_name":"GRCh38"},{"assembly_name":"GRCh38","species":"homo_sapiens","end":151055993,"db_type":"core","object_type":"Exon","version":1,"id":"ENSE00003678094","seq_region_name":"1","strand":-1,"start":151055677},{"seq_region_name":"1","version":1,"id":"ENSE00003609119","start":151055016,"strand":-1,"assembly_name":"GRCh38","end":151055126,"species":"homo_sapiens","db_type":"core","object_type":"Exon"},{"db_type":"core","object_type":"Exon","species":"homo_sapiens","end":151054321,"assembly_name":"GRCh38","start":151054237,"strand":-1,"version":1,"id":"ENSE00001944752","seq_region_name":"1"}],"is_canonical":0,"Parent":"ENSG00000197622","display_name":"CDC42SE1-203","species":"homo_sapiens","end":151059773,"biotype":"protein_coding_CDS_not_defined","length":939,"assembly_name":"GRCh38","db_type":"core","object_type":"Transcript","UTR":[],"id":"ENST00000470278","version":5,"logic_name":"havana_homo_sapiens","seq_region_name":"1","source":"havana","strand":-1,"start":151054237},"ENST00000404930":{"assembly_name":"GRCh38","db_type":"core","object_type":"Transcript","UTR":[],"seq_region_name":"6","logic_name":"havana_homo_sapiens","version":1,"id":"ENST00000404930","source":"havana","start":105666326,"strand":1,"Exon":[{"db_type":"core","object_type":"Exon","species":"homo_sapiens","end":105667998,"assembly_name":"GRCh38","start":105666326,"strand":1,"id":"ENSE00001552310","version":1,"seq_region_name":"6"}],"is_canonical":1,"Parent":"ENSG00000219088","end":105667998,"species":"homo_sapiens","biotype":"processed_pseudogene","length":1673},"ENST00000434805":{"is_canonical":1,"Exon":[{"version":1,"id":"ENSE00001756882","seq_region_name":"1","start":35350722,"strand":1,"assembly_name":"GRCh38","db_type":"core","object_type":"Exon","species":"homo_sapiens","end":35351607}],"end":35351607,"species":"homo_sapiens","display_name":"RPL5P4-201","Parent":"ENSG00000229994","length":886,"biotype":"processed_pseudogene","assembly_name":"GRCh38","UTR":[],"object_type":"Transcript","db_type":"core","source":"havana","seq_region_name":"1","logic_name":"havana_homo_sapiens","version":1,"id":"ENST00000434805","strand":1,"start":35350722},"ENST00000314607":{"species":"homo_sapiens","end":35422058,"display_name":"ZMYM4-201","Parent":"ENSG00000146463","is_canonical":1,"Exon":[{"strand":1,"start":35268709,"seq_region_name":"1","version":2,"id":"ENSE00001670766","object_type":"Exon","db_type":"core","end":35269085,"species":"homo_sapiens","assembly_name":"GRCh38"},{"seq_region_name":"1","id":"ENSE00003615141","version":1,"start":35325360,"strand":1,"assembly_name":"GRCh38","db_type":"core","object_type":"Exon","end":35325405,"species":"homo_sapiens"},{"start":35358925,"strand":1,"seq_region_name":"1","id":"ENSE00001765619","version":1,"end":35359446,"species":"homo_sapiens","db_type":"core","object_type":"Exon","assembly_name":"GRCh38"},{"species":"homo_sapiens","end":35361255,"db_type":"core","object_type":"Exon","assembly_name":"GRCh38","strand":1,"start":35361194,"id":"ENSE00001425567","version":1,"seq_region_name":"1"},{"seq_region_name":"1","version":1,"id":"ENSE00001417767","start":35361619,"strand":1,"assembly_name":"GRCh38","end":35361789,"species":"homo_sapiens","db_type":"core","object_type":"Exon"},{"start":35370029,"strand":1,"version":1,"id":"ENSE00003478009","seq_region_name":"1","object_type":"Exon","db_type":"core","species":"homo_sapiens","end":35370113,"assembly_name":"GRCh38"},{"db_type":"core","object_type":"Exon","end":35370627,"species":"homo_sapiens","assembly_name":"GRCh38","start":35370372,"strand":1,"seq_region_name":"1","version":1,"id":"ENSE00001429123"},{"db_type":"core","object_type":"Exon","species":"homo_sapiens","end":35381433,"assembly_name":"GRCh38","strand":1,"start":35381259,"version":1,"id":"ENSE00001125401","seq_region_name":"1"},{"start":35381546,"strand":1,"version":1,"id":"ENSE00001066929","seq_region_name":"1","species":"homo_sapiens","end":35381758,"object_type":"Exon","db_type":"core","assembly_name":"GRCh38"},{"version":1,"id":"ENSE00001066912","seq_region_name":"1","strand":1,"start":35385442,"assembly_name":"GRCh38","db_type":"core","object_type":"Exon","species":"homo_sapiens","end":35385592},{"version":1,"id":"ENSE00001125375","seq_region_name":"1","strand":1,"start":35386074,"assembly_name":"GRCh38","species":"homo_sapiens","end":35386189,"object_type":"Exon","db_type":"core"},{"end":35387278,"species":"homo_sapiens","object_type":"Exon","db_type":"core","assembly_name":"GRCh38","strand":1,"start":35387003,"seq_region_name":"1","version":1,"id":"ENSE00001066919"},{"start":35387454,"strand":1,"id":"ENSE00001616388","version":1,"seq_region_name":"1","species":"homo_sapiens","end":35387604,"db_type":"core","object_type":"Exon","assembly_name":"GRCh38"},{"assembly_name":"GRCh38","species":"homo_sapiens","end":35389082,"db_type":"core","object_type":"Exon","id":"ENSE00001125352","version":1,"seq_region_name":"1","strand":1,"start":35388910},{"assembly_name":"GRCh38","species":"homo_sapiens","end":35390098,"object_type":"Exon","db_type":"core","id":"ENSE00001066932","version":1,"seq_region_name":"1","strand":1,"start":35389948},{"db_type":"core","object_type":"Exon","end":35392352,"species":"homo_sapiens","assembly_name":"GRCh38","start":35392212,"strand":1,"seq_region_name":"1","version":1,"id":"ENSE00001125338"},{"strand":1,"start":35392647,"seq_region_name":"1","id":"ENSE00001125330","version":1,"end":35392684,"species":"homo_sapiens","object_type":"Exon","db_type":"core","assembly_name":"GRCh38"},{"strand":1,"start":35393595,"seq_region_name":"1","version":1,"id":"ENSE00003474598","end":35393739,"species":"homo_sapiens","object_type":"Exon","db_type":"core","assembly_name":"GRCh38"},{"seq_region_name":"1","version":1,"id":"ENSE00003482654","strand":1,"start":35396552,"assembly_name":"GRCh38","end":35396670,"species":"homo_sapiens","db_type":"core","object_type":"Exon"},{"end":35397545,"species":"homo_sapiens","object_type":"Exon","db_type":"core","assembly_name":"GRCh38","strand":1,"start":35397377,"seq_region_name":"1","version":1,"id":"ENSE00001066928"},{"seq_region_name":"1","id":"ENSE00001125294","version":1,"strand":1,"start":35398413,"assembly_name":"GRCh38","end":35398466,"species":"homo_sapiens","object_type":"Exon","db_type":"core"},{"version":1,"id":"ENSE00001125285","seq_region_name":"1","strand":1,"start":35398864,"assembly_name":"GRCh38","species":"homo_sapiens","end":35399043,"object_type":"Exon","db_type":"core"},{"seq_region_name":"1","version":1,"id":"ENSE00001125276","strand":1,"start":35399482,"assembly_name":"GRCh38","end":35399576,"species":"homo_sapiens","object_type":"Exon","db_type":"core"},{"species":"homo_sapiens","end":35405194,"db_type":"core","object_type":"Exon","assembly_name":"GRCh38","start":35405023,"strand":1,"id":"ENSE00001066906","version":1,"seq_region_name":"1"},{"db_type":"core","object_type":"Exon","species":"homo_sapiens","end":35405468,"assembly_name":"GRCh38","strand":1,"start":35405373,"id":"ENSE00003688215","version":1,"seq_region_name":"1"},{"object_type":"Exon","db_type":"core","species":"homo_sapiens","end":35408159,"assembly_name":"GRCh38","start":35408008,"strand":1,"id":"ENSE00001066914","version":1,"seq_region_name":"1"},{"db_type":"core","object_type":"Exon","species":"homo_sapiens","end":35414083,"assembly_name":"GRCh38","strand":1,"start":35413972,"version":1,"id":"ENSE00001066931","seq_region_name":"1"},{"seq_region_name":"1","id":"ENSE00001125229","version":1,"start":35415466,"strand":1,"assembly_name":"GRCh38","end":35415714,"species":"homo_sapiens","db_type":"core","object_type":"Exon"},{"start":35418443,"strand":1,"seq_region_name":"1","id":"ENSE00000955938","version":1,"db_type":"core","object_type":"Exon","end":35418572,"species":"homo_sapiens","assembly_name":"GRCh38"},{"species":"homo_sapiens","end":35422058,"db_type":"core","object_type":"Exon","assembly_name":"GRCh38","strand":1,"start":35419470,"version":1,"id":"ENSE00001626559","seq_region_name":"1"}],"length":7366,"biotype":"protein_coding","UTR":[{"start":35268709,"strand":1,"source":"ensembl_havana","type":"five_prime_utr","seq_region_name":"1","id":"ENST00000314607","end":35269046,"species":"homo_sapiens","Parent":"ENST00000314607","object_type":"five_prime_UTR","db_type":"core","assembly_name":"GRCh38"},{"id":"ENST00000314607","seq_region_name":"1","type":"three_prime_utr","source":"ensembl_havana","start":35419678,"strand":1,"assembly_name":"GRCh38","db_type":"core","object_type":"three_prime_UTR","Parent":"ENST00000314607","species":"homo_sapiens","end":35422058}],"object_type":"Transcript","db_type":"core","assembly_name":"GRCh38","start":35268709,"strand":1,"source":"ensembl_havana","version":11,"id":"ENST00000314607","seq_region_name":"1","Translation":{"species":"homo_sapiens","end":35419677,"db_type":"core","object_type":"Translation","Parent":"ENST00000314607","length":1548,"version":6,"id":"ENSP00000322915","start":35269047},"logic_name":"ensembl_havana_transcript_homo_sapiens"},"ENST00000441447":{"db_type":"core","object_type":"Transcript","UTR":[{"end":35269085,"species":"homo_sapiens","Parent":"ENST00000441447","object_type":"five_prime_UTR","db_type":"core","assembly_name":"GRCh38","strand":1,"start":35269033,"type":"five_prime_utr","source":"havana","seq_region_name":"1","id":"ENST00000441447"},{"source":"havana","type":"five_prime_utr","seq_region_name":"1","id":"ENST00000441447","start":35295927,"strand":1,"assembly_name":"GRCh38","end":35295988,"species":"homo_sapiens","Parent":"ENST00000441447","object_type":"five_prime_UTR","db_type":"core"},{"strand":1,"start":35325360,"id":"ENST00000441447","seq_region_name":"1","source":"havana","type":"five_prime_utr","db_type":"core","object_type":"five_prime_UTR","Parent":"ENST00000441447","species":"homo_sapiens","end":35325405,"assembly_name":"GRCh38"},{"strand":1,"start":35358925,"id":"ENST00000441447","seq_region_name":"1","type":"five_prime_utr","source":"havana","db_type":"core","object_type":"five_prime_UTR","Parent":"ENST00000441447","species":"homo_sapiens","end":35358935,"assembly_name":"GRCh38"}],"assembly_name":"GRCh38","strand":1,"start":35269033,"seq_region_name":"1","Translation":{"end":35359229,"species":"homo_sapiens","Parent":"ENST00000441447","db_type":"core","object_type":"Translation","length":98,"id":"ENSP00000397524","version":1,"start":35358936},"logic_name":"havana_homo_sapiens","version":1,"id":"ENST00000441447","source":"havana","Parent":"ENSG00000146463","display_name":"ZMYM4-202","end":35359229,"species":"homo_sapiens","Exon":[{"db_type":"core","object_type":"Exon","species":"homo_sapiens","end":35269085,"assembly_name":"GRCh38","strand":1,"start":35269033,"version":1,"id":"ENSE00001792125","seq_region_name":"1"},{"strand":1,"start":35295927,"seq_region_name":"1","version":1,"id":"ENSE00001649553","end":35295988,"species":"homo_sapiens","object_type":"Exon","db_type":"core","assembly_name":"GRCh38"},{"object_type":"Exon","db_type":"core","end":35325405,"species":"homo_sapiens","assembly_name":"GRCh38","strand":1,"start":35325360,"seq_region_name":"1","id":"ENSE00003616828","version":1},{"id":"ENSE00001716452","version":1,"seq_region_name":"1","strand":1,"start":35358925,"assembly_name":"GRCh38","species":"homo_sapiens","end":35359229,"db_type":"core","object_type":"Exon"}],"is_canonical":0,"biotype":"protein_coding","length":466},"ENST00000435237":{"id":"ENST00000435237","version":1,"logic_name":"havana_homo_sapiens","seq_region_name":"2","source":"havana","strand":1,"start":21221175,"assembly_name":"GRCh38","db_type":"core","object_type":"Transcript","UTR":[],"biotype":"lncRNA","length":488,"Exon":[{"end":21221294,"species":"homo_sapiens","db_type":"core","object_type":"Exon","assembly_name":"GRCh38","start":21221175,"strand":1,"seq_region_name":"2","version":1,"id":"ENSE00001630951"},{"assembly_name":"GRCh38","end":21527548,"species":"homo_sapiens","db_type":"core","object_type":"Exon","seq_region_name":"2","id":"ENSE00001776775","version":1,"strand":1,"start":21527507},{"assembly_name":"GRCh38","species":"homo_sapiens","end":21529243,"object_type":"Exon","db_type":"core","version":1,"id":"ENSE00001640430","seq_region_name":"2","strand":1,"start":21529213},{"seq_region_name":"2","version":1,"id":"ENSE00001723765","strand":1,"start":21799357,"assembly_name":"GRCh38","db_type":"core","object_type":"Exon","end":21799478,"species":"homo_sapiens"},{"strand":1,"start":21951632,"seq_region_name":"2","version":1,"id":"ENSE00001714190","object_type":"Exon","db_type":"core","end":21951689,"species":"homo_sapiens","assembly_name":"GRCh38"},{"assembly_name":"GRCh38","object_type":"Exon","db_type":"core","end":21970959,"species":"homo_sapiens","seq_region_name":"2","id":"ENSE00001677049","version":1,"start":21970845,"strand":1}],"is_canonical":1,"Parent":"ENSG00000233005","species":"homo_sapiens","end":21970959},"ENST00000402318":{"display_name":"ANKRD20A7P-201","Parent":"ENSG00000236816","species":"homo_sapiens","end":42920095,"Exon":[{"version":1,"id":"ENSE00001647847","seq_region_name":"9","start":42852675,"strand":1,"assembly_name":"GRCh38","species":"homo_sapiens","end":42852877,"db_type":"core","object_type":"Exon"},{"strand":1,"start":42856198,"id":"ENSE00001707608","version":1,"seq_region_name":"9","species":"homo_sapiens","end":42856312,"object_type":"Exon","db_type":"core","assembly_name":"GRCh38"},{"id":"ENSE00001691736","version":1,"seq_region_name":"9","start":42856475,"strand":1,"assembly_name":"GRCh38","object_type":"Exon","db_type":"core","species":"homo_sapiens","end":42856648},{"seq_region_name":"9","version":1,"id":"ENSE00001671768","start":42860536,"strand":1,"assembly_name":"GRCh38","end":42860643,"species":"homo_sapiens","db_type":"core","object_type":"Exon"},{"end":42861818,"species":"homo_sapiens","object_type":"Exon","db_type":"core","assembly_name":"GRCh38","strand":1,"start":42861684,"seq_region_name":"9","id":"ENSE00001757785","version":1},{"assembly_name":"GRCh38","species":"homo_sapiens","end":42864465,"object_type":"Exon","db_type":"core","id":"ENSE00001725046","version":1,"seq_region_name":"9","start":42864410,"strand":1},{"seq_region_name":"9","version":1,"id":"ENSE00001654421","start":42871003,"strand":1,"assembly_name":"GRCh38","end":42871033,"species":"homo_sapiens","object_type":"Exon","db_type":"core"},{"species":"homo_sapiens","end":42873790,"object_type":"Exon","db_type":"core","assembly_name":"GRCh38","start":42873721,"strand":1,"id":"ENSE00001643276","version":1,"seq_region_name":"9"},{"end":42877817,"species":"homo_sapiens","db_type":"core","object_type":"Exon","assembly_name":"GRCh38","strand":1,"start":42877733,"seq_region_name":"9","id":"ENSE00001714906","version":1},{"end":42880658,"species":"homo_sapiens","object_type":"Exon","db_type":"core","assembly_name":"GRCh38","start":42880630,"strand":1,"seq_region_name":"9","version":1,"id":"ENSE00001756444"},{"seq_region_name":"9","id":"ENSE00001641577","version":1,"start":42880750,"strand":1,"assembly_name":"GRCh38","object_type":"Exon","db_type":"core","end":42880822,"species":"homo_sapiens"},{"version":1,"id":"ENSE00001738184","seq_region_name":"9","start":42886778,"strand":1,"assembly_name":"GRCh38","db_type":"core","object_type":"Exon","species":"homo_sapiens","end":42886848},{"strand":1,"start":42890906,"version":1,"id":"ENSE00001753332","seq_region_name":"9","db_type":"core","object_type":"Exon","species":"homo_sapiens","end":42891069,"assembly_name":"GRCh38"},{"db_type":"core","object_type":"Exon","species":"homo_sapiens","end":42892665,"assembly_name":"GRCh38","start":42892485,"strand":1,"id":"ENSE00001635026","version":1,"seq_region_name":"9"},{"assembly_name":"GRCh38","end":42894754,"species":"homo_sapiens","db_type":"core","object_type":"Exon","seq_region_name":"9","id":"ENSE00001738231","version":2,"strand":1,"start":42893833},{"object_type":"Exon","db_type":"core","end":42896450,"species":"homo_sapiens","assembly_name":"GRCh38","start":42896301,"strand":1,"seq_region_name":"9","id":"ENSE00003878119","version":1},{"strand":1,"start":42901238,"id":"ENSE00001724094","version":1,"seq_region_name":"9","object_type":"Exon","db_type":"core","species":"homo_sapiens","end":42901355,"assembly_name":"GRCh38"},{"seq_region_name":"9","version":1,"id":"ENSE00002260074","strand":1,"start":42903343,"assembly_name":"GRCh38","db_type":"core","object_type":"Exon","end":42903472,"species":"homo_sapiens"},{"assembly_name":"GRCh38","species":"homo_sapiens","end":42911013,"db_type":"core","object_type":"Exon","version":1,"id":"ENSE00001803648","seq_region_name":"9","strand":1,"start":42910732},{"assembly_name":"GRCh38","object_type":"Exon","db_type":"core","end":42911815,"species":"homo_sapiens","seq_region_name":"9","id":"ENSE00001676679","version":1,"strand":1,"start":42911604},{"seq_region_name":"9","id":"ENSE00001673441","version":1,"start":42913734,"strand":1,"assembly_name":"GRCh38","end":42913956,"species":"homo_sapiens","db_type":"core","object_type":"Exon"},{"assembly_name":"GRCh38","db_type":"core","object_type":"Exon","end":42917416,"species":"homo_sapiens","seq_region_name":"9","version":1,"id":"ENSE00001637688","start":42917367,"strand":1},{"seq_region_name":"9","id":"ENSE00001742137","version":1,"start":42918583,"strand":1,"assembly_name":"GRCh38","db_type":"core","object_type":"Exon","end":42918858,"species":"homo_sapiens"},{"start":42919688,"strand":1,"seq_region_name":"9","id":"ENSE00001761673","version":1,"object_type":"Exon","db_type":"core","end":42920095,"species":"homo_sapiens","assembly_name":"GRCh38"}],"is_canonical":1,"biotype":"transcribed_unprocessed_pseudogene","length":4266,"db_type":"core","object_type":"Transcript","UTR":[],"assembly_name":"GRCh38","start":42852675,"strand":1,"id":"ENST00000402318","version":3,"seq_region_name":"9","logic_name":"havana_homo_sapiens","source":"havana"},"ENST00000342888":{"length":1859,"biotype":"lncRNA","end":18757894,"species":"homo_sapiens","display_name":"FAM230E-201","Parent":"ENSG00000182824","is_canonical":1,"Exon":[{"start":18733914,"strand":1,"seq_region_name":"22","version":1,"id":"ENSE00002220770","end":18734054,"species":"homo_sapiens","object_type":"Exon","db_type":"core","assembly_name":"GRCh38"},{"id":"ENSE00001596564","version":1,"seq_region_name":"22","start":18736032,"strand":1,"assembly_name":"GRCh38","species":"homo_sapiens","end":18736090,"object_type":"Exon","db_type":"core"},{"assembly_name":"GRCh38","db_type":"core","object_type":"Exon","species":"homo_sapiens","end":18736583,"id":"ENSE00001763149","version":1,"seq_region_name":"22","start":18736559,"strand":1},{"assembly_name":"GRCh38","species":"homo_sapiens","end":18739616,"db_type":"core","object_type":"Exon","id":"ENSE00001670476","version":1,"seq_region_name":"22","strand":1,"start":18739530},{"version":1,"id":"ENSE00001729746","seq_region_name":"22","strand":1,"start":18744510,"assembly_name":"GRCh38","species":"homo_sapiens","end":18744575,"db_type":"core","object_type":"Exon"},{"version":1,"id":"ENSE00001693656","seq_region_name":"22","start":18746574,"strand":1,"assembly_name":"GRCh38","db_type":"core","object_type":"Exon","species":"homo_sapiens","end":18746635},{"assembly_name":"GRCh38","object_type":"Exon","db_type":"core","species":"homo_sapiens","end":18747134,"id":"ENSE00001744915","version":1,"seq_region_name":"22","strand":1,"start":18747101},{"assembly_name":"GRCh38","db_type":"core","object_type":"Exon","species":"homo_sapiens","end":18752825,"id":"ENSE00001642671","version":1,"seq_region_name":"22","strand":1,"start":18751869},{"assembly_name":"GRCh38","db_type":"core","object_type":"Exon","end":18757894,"species":"homo_sapiens","seq_region_name":"22","version":1,"id":"ENSE00001608479","strand":1,"start":18757467}],"start":18733914,"strand":1,"source":"havana","logic_name":"havana_homo_sapiens","seq_region_name":"22","id":"ENST00000342888","version":3,"UTR":[],"object_type":"Transcript","db_type":"core","assembly_name":"GRCh38"},"ENST00000483763":{"UTR":[],"object_type":"Transcript","db_type":"core","assembly_name":"GRCh38","start":151052114,"strand":-1,"source":"havana","seq_region_name":"1","logic_name":"havana_homo_sapiens","version":5,"id":"ENST00000483763","end":151059574,"species":"homo_sapiens","display_name":"CDC42SE1-204","Parent":"ENSG00000197622","is_canonical":0,"Exon":[{"species":"homo_sapiens","end":151059574,"object_type":"Exon","db_type":"core","assembly_name":"GRCh38","start":151059479,"strand":-1,"id":"ENSE00001944190","version":1,"seq_region_name":"1"},{"assembly_name":"GRCh38","species":"homo_sapiens","end":151055993,"object_type":"Exon","db_type":"core","id":"ENSE00001954293","version":1,"seq_region_name":"1","strand":-1,"start":151055016},{"db_type":"core","object_type":"Exon","end":151054321,"species":"homo_sapiens","assembly_name":"GRCh38","start":151054231,"strand":-1,"seq_region_name":"1","id":"ENSE00003460320","version":1},{"start":151052114,"strand":-1,"seq_region_name":"1","id":"ENSE00001809619","version":1,"object_type":"Exon","db_type":"core","end":151053327,"species":"homo_sapiens","assembly_name":"GRCh38"}],"length":2379,"biotype":"retained_intron"},"ENST00000668653":{"Parent":"ENSG00000291166","species":"homo_sapiens","end":42914197,"Exon":[{"assembly_name":"GRCh38","db_type":"core","object_type":"Exon","end":42845420,"species":"homo_sapiens","seq_region_name":"9","version":1,"id":"ENSE00003863616","strand":1,"start":42845354},{"version":1,"id":"ENSE00004021096","seq_region_name":"9","strand":1,"start":42856198,"assembly_name":"GRCh38","db_type":"core","object_type":"Exon","species":"homo_sapiens","end":42856312},{"version":1,"id":"ENSE00004021053","seq_region_name":"9","start":42856475,"strand":1,"assembly_name":"GRCh38","db_type":"core","object_type":"Exon","species":"homo_sapiens","end":42856648},{"version":1,"id":"ENSE00003858241","seq_region_name":"9","strand":1,"start":42856952,"assembly_name":"GRCh38","species":"homo_sapiens","end":42857107,"object_type":"Exon","db_type":"core"},{"object_type":"Exon","db_type":"core","end":42860643,"species":"homo_sapiens","assembly_name":"GRCh38","start":42860536,"strand":1,"seq_region_name":"9","id":"ENSE00004021015","version":1},{"species":"homo_sapiens","end":42861818,"object_type":"Exon","db_type":"core","assembly_name":"GRCh38","start":42861684,"strand":1,"id":"ENSE00004021200","version":1,"seq_region_name":"9"},{"end":42864465,"species":"homo_sapiens","object_type":"Exon","db_type":"core","assembly_name":"GRCh38","strand":1,"start":42864410,"seq_region_name":"9","version":1,"id":"ENSE00004021130"},{"version":1,"id":"ENSE00004020983","seq_region_name":"9","start":42871003,"strand":1,"assembly_name":"GRCh38","object_type":"Exon","db_type":"core","species":"homo_sapiens","end":42871033},{"strand":1,"start":42873721,"id":"ENSE00004020959","version":1,"seq_region_name":"9","species":"homo_sapiens","end":42873790,"db_type":"core","object_type":"Exon","assembly_name":"GRCh38"},{"strand":1,"start":42877736,"seq_region_name":"9","id":"ENSE00003886473","version":1,"object_type":"Exon","db_type":"core","end":42877817,"species":"homo_sapiens","assembly_name":"GRCh38"},{"start":42880630,"strand":1,"id":"ENSE00004021196","version":1,"seq_region_name":"9","species":"homo_sapiens","end":42880658,"db_type":"core","object_type":"Exon","assembly_name":"GRCh38"},{"assembly_name":"GRCh38","object_type":"Exon","db_type":"core","species":"homo_sapiens","end":42880822,"id":"ENSE00004020956","version":1,"seq_region_name":"9","start":42880750,"strand":1},{"assembly_name":"GRCh38","species":"homo_sapiens","end":42886848,"db_type":"core","object_type":"Exon","id":"ENSE00004021157","version":1,"seq_region_name":"9","strand":1,"start":42886778},{"seq_region_name":"9","version":1,"id":"ENSE00004021188","start":42890906,"strand":1,"assembly_name":"GRCh38","object_type":"Exon","db_type":"core","end":42891069,"species":"homo_sapiens"},{"id":"ENSE00004020938","version":1,"seq_region_name":"9","start":42892485,"strand":1,"assembly_name":"GRCh38","species":"homo_sapiens","end":42892665,"object_type":"Exon","db_type":"core"},{"assembly_name":"GRCh38","species":"homo_sapiens","end":42894754,"object_type":"Exon","db_type":"core","id":"ENSE00004021158","version":1,"seq_region_name":"9","start":42893833,"strand":1},{"species":"homo_sapiens","end":42896450,"object_type":"Exon","db_type":"core","assembly_name":"GRCh38","strand":1,"start":42896301,"id":"ENSE00004024184","version":1,"seq_region_name":"9"},{"assembly_name":"GRCh38","species":"homo_sapiens","end":42903472,"db_type":"core","object_type":"Exon","id":"ENSE00004021749","version":1,"seq_region_name":"9","strand":1,"start":42903343},{"end":42911013,"species":"homo_sapiens","db_type":"core","object_type":"Exon","assembly_name":"GRCh38","strand":1,"start":42910732,"seq_region_name":"9","id":"ENSE00004021288","version":1},{"start":42911604,"strand":1,"id":"ENSE00004021025","version":1,"seq_region_name":"9","species":"homo_sapiens","end":42911815,"db_type":"core","object_type":"Exon","assembly_name":"GRCh38"},{"strand":1,"start":42913734,"seq_region_name":"9","version":1,"id":"ENSE00003857699","db_type":"core","object_type":"Exon","end":42914197,"species":"homo_sapiens","assembly_name":"GRCh38"}],"is_canonical":1,"biotype":"lncRNA","length":3672,"db_type":"core","object_type":"Transcript","UTR":[],"assembly_name":"GRCh38","strand":1,"start":42845354,"id":"ENST00000668653","version":1,"logic_name":"havana_homo_sapiens","seq_region_name":"9","source":"havana"},"ENST00000492796":{"assembly_name":"GRCh38","db_type":"core","object_type":"Transcript","UTR":[],"seq_region_name":"1","logic_name":"havana_homo_sapiens","id":"ENST00000492796","version":5,"source":"havana","strand":-1,"start":151052946,"Exon":[{"db_type":"core","object_type":"Exon","end":151059615,"species":"homo_sapiens","assembly_name":"GRCh38","start":151059479,"strand":-1,"seq_region_name":"1","version":1,"id":"ENSE00001875320"},{"version":1,"id":"ENSE00001513184","seq_region_name":"1","strand":-1,"start":151056651,"assembly_name":"GRCh38","object_type":"Exon","db_type":"core","species":"homo_sapiens","end":151056786},{"species":"homo_sapiens","end":151055993,"db_type":"core","object_type":"Exon","assembly_name":"GRCh38","start":151055677,"strand":-1,"id":"ENSE00003678094","version":1,"seq_region_name":"1"},{"assembly_name":"GRCh38","object_type":"Exon","db_type":"core","end":151055126,"species":"homo_sapiens","seq_region_name":"1","version":1,"id":"ENSE00003609119","start":151055016,"strand":-1},{"assembly_name":"GRCh38","db_type":"core","object_type":"Exon","end":151054321,"species":"homo_sapiens","seq_region_name":"1","id":"ENSE00003460320","version":1,"strand":-1,"start":151054231},{"id":"ENSE00001888685","version":1,"seq_region_name":"1","start":151052946,"strand":-1,"assembly_name":"GRCh38","species":"homo_sapiens","end":151053327,"db_type":"core","object_type":"Exon"}],"is_canonical":0,"display_name":"CDC42SE1-206","Parent":"ENSG00000197622","end":151059615,"species":"homo_sapiens","biotype":"protein_coding_CDS_not_defined","length":1174},"ENST00000666959":{"is_canonical":0,"Exon":[{"start":42892482,"strand":1,"id":"ENSE00003851892","version":1,"seq_region_name":"9","object_type":"Exon","db_type":"core","species":"homo_sapiens","end":42892665,"assembly_name":"GRCh38"},{"strand":1,"start":42893833,"id":"ENSE00004021158","version":1,"seq_region_name":"9","species":"homo_sapiens","end":42894754,"object_type":"Exon","db_type":"core","assembly_name":"GRCh38"},{"assembly_name":"GRCh38","end":42896450,"species":"homo_sapiens","object_type":"Exon","db_type":"core","seq_region_name":"9","version":1,"id":"ENSE00004024184","start":42896301,"strand":1},{"object_type":"Exon","db_type":"core","end":42903472,"species":"homo_sapiens","assembly_name":"GRCh38","start":42903343,"strand":1,"seq_region_name":"9","id":"ENSE00004021749","version":1},{"start":42911604,"strand":1,"seq_region_name":"9","id":"ENSE00004021025","version":1,"end":42911815,"species":"homo_sapiens","db_type":"core","object_type":"Exon","assembly_name":"GRCh38"},{"version":1,"id":"ENSE00004021019","seq_region_name":"9","start":42913734,"strand":1,"assembly_name":"GRCh38","object_type":"Exon","db_type":"core","species":"homo_sapiens","end":42913956},{"seq_region_name":"9","id":"ENSE00003865979","version":1,"start":42917367,"strand":1,"assembly_name":"GRCh38","db_type":"core","object_type":"Exon","end":42917412,"species":"homo_sapiens"},{"strand":1,"start":42918583,"seq_region_name":"9","version":1,"id":"ENSE00003872866","end":42918636,"species":"homo_sapiens","db_type":"core","object_type":"Exon","assembly_name":"GRCh38"},{"assembly_name":"GRCh38","db_type":"core","object_type":"Exon","species":"homo_sapiens","end":42920381,"id":"ENSE00003860663","version":1,"seq_region_name":"9","strand":1,"start":42920263},{"version":1,"id":"ENSE00003852778","seq_region_name":"9","strand":1,"start":42950090,"assembly_name":"GRCh38","species":"homo_sapiens","end":42950197,"db_type":"core","object_type":"Exon"},{"seq_region_name":"9","version":1,"id":"ENSE00003878804","strand":1,"start":42950566,"assembly_name":"GRCh38","object_type":"Exon","db_type":"core","end":42950827,"species":"homo_sapiens"}],"end":42950827,"species":"homo_sapiens","Parent":"ENSG00000291166","length":2410,"biotype":"lncRNA","assembly_name":"GRCh38","UTR":[],"object_type":"Transcript","db_type":"core","source":"havana","seq_region_name":"9","logic_name":"havana_homo_sapiens","version":1,"id":"ENST00000666959","start":42892482,"strand":1},"ENST00000540998":{"species":"homo_sapiens","end":151059649,"Parent":"ENSG00000197622","display_name":"CDC42SE1-207","is_canonical":0,"Exon":[{"object_type":"Exon","db_type":"core","end":151059649,"species":"homo_sapiens","assembly_name":"GRCh38","strand":-1,"start":151059479,"seq_region_name":"1","id":"ENSE00002253899","version":1},{"assembly_name":"GRCh38","object_type":"Exon","db_type":"core","end":151056786,"species":"homo_sapiens","seq_region_name":"1","id":"ENSE00001513184","version":1,"start":151056651,"strand":-1},{"strand":-1,"start":151055677,"version":1,"id":"ENSE00003519252","seq_region_name":"1","species":"homo_sapiens","end":151055993,"db_type":"core","object_type":"Exon","assembly_name":"GRCh38"},{"version":1,"id":"ENSE00003604581","seq_region_name":"1","strand":-1,"start":151055016,"assembly_name":"GRCh38","db_type":"core","object_type":"Exon","species":"homo_sapiens","end":151055126},{"seq_region_name":"1","id":"ENSE00003684539","version":1,"strand":-1,"start":151054231,"assembly_name":"GRCh38","db_type":"core","object_type":"Exon","end":151054321,"species":"homo_sapiens"},{"assembly_name":"GRCh38","end":151053327,"species":"homo_sapiens","db_type":"core","object_type":"Exon","seq_region_name":"1","id":"ENSE00002319572","version":1,"start":151050973,"strand":-1}],"length":3181,"biotype":"protein_coding","UTR":[{"strand":-1,"start":151059479,"seq_region_name":"1","id":"ENST00000540998","type":"five_prime_utr","source":"ensembl","Parent":"ENST00000540998","db_type":"core","object_type":"five_prime_UTR","end":151059649,"species":"homo_sapiens","assembly_name":"GRCh38"},{"end":151056786,"species":"homo_sapiens","Parent":"ENST00000540998","object_type":"five_prime_UTR","db_type":"core","assembly_name":"GRCh38","strand":-1,"start":151056651,"source":"ensembl","type":"five_prime_utr","seq_region_name":"1","id":"ENST00000540998"},{"start":151055731,"strand":-1,"id":"ENST00000540998","seq_region_name":"1","type":"five_prime_utr","source":"ensembl","object_type":"five_prime_UTR","db_type":"core","Parent":"ENST00000540998","species":"homo_sapiens","end":151055993,"assembly_name":"GRCh38"},{"db_type":"core","object_type":"three_prime_UTR","Parent":"ENST00000540998","species":"homo_sapiens","end":151054246,"assembly_name":"GRCh38","strand":-1,"start":151054231,"id":"ENST00000540998","seq_region_name":"1","source":"ensembl","type":"three_prime_utr"},{"strand":-1,"start":151050973,"id":"ENST00000540998","seq_region_name":"1","type":"three_prime_utr","source":"ensembl","object_type":"three_prime_UTR","db_type":"core","Parent":"ENST00000540998","species":"homo_sapiens","end":151053327,"assembly_name":"GRCh38"}],"db_type":"core","object_type":"Transcript","assembly_name":"GRCh38","start":151050973,"strand":-1,"source":"ensembl","version":5,"id":"ENST00000540998","Translation":{"length":79,"id":"ENSP00000445647","version":1,"start":151054247,"end":151055730,"species":"homo_sapiens","Parent":"ENST00000540998","db_type":"core","object_type":"Translation"},"logic_name":"ensembl_homo_sapiens","seq_region_name":"1"},"ENST00000439374":{"assembly_name":"GRCh38","object_type":"Transcript","db_type":"core","UTR":[{"seq_region_name":"1","id":"ENST00000439374","type":"five_prime_utr","source":"havana","start":151070264,"strand":-1,"assembly_name":"GRCh38","Parent":"ENST00000439374","object_type":"five_prime_UTR","db_type":"core","end":151070325,"species":"homo_sapiens"},{"Parent":"ENST00000439374","db_type":"core","object_type":"five_prime_UTR","end":151068380,"species":"homo_sapiens","assembly_name":"GRCh38","strand":-1,"start":151068319,"seq_region_name":"1","id":"ENST00000439374","type":"five_prime_utr","source":"havana"},{"strand":-1,"start":151067132,"id":"ENST00000439374","seq_region_name":"1","source":"havana","type":"five_prime_utr","object_type":"five_prime_UTR","db_type":"core","Parent":"ENST00000439374","species":"homo_sapiens","end":151067334,"assembly_name":"GRCh38"},{"Parent":"ENST00000439374","object_type":"five_prime_UTR","db_type":"core","end":151059773,"species":"homo_sapiens","assembly_name":"GRCh38","strand":-1,"start":151059479,"seq_region_name":"1","id":"ENST00000439374","source":"havana","type":"five_prime_utr"},{"Parent":"ENST00000439374","db_type":"core","object_type":"five_prime_UTR","end":151055993,"species":"homo_sapiens","assembly_name":"GRCh38","start":151055731,"strand":-1,"seq_region_name":"1","id":"ENST00000439374","source":"havana","type":"five_prime_utr"},{"end":151054246,"species":"homo_sapiens","Parent":"ENST00000439374","db_type":"core","object_type":"three_prime_UTR","assembly_name":"GRCh38","start":151054231,"strand":-1,"source":"havana","type":"three_prime_utr","seq_region_name":"1","id":"ENST00000439374"},{"assembly_name":"GRCh38","species":"homo_sapiens","end":151053327,"db_type":"core","object_type":"three_prime_UTR","Parent":"ENST00000439374","type":"three_prime_utr","source":"havana","id":"ENST00000439374","seq_region_name":"1","strand":-1,"start":151050971}],"id":"ENST00000439374","version":6,"Translation":{"length":79,"version":1,"id":"ENSP00000475845","start":151054247,"species":"homo_sapiens","end":151055730,"object_type":"Translation","db_type":"core","Parent":"ENST00000439374"},"seq_region_name":"1","logic_name":"havana_homo_sapiens","source":"havana","start":151050971,"strand":-1,"Exon":[{"seq_region_name":"1","version":1,"id":"ENSE00001786108","start":151070264,"strand":-1,"assembly_name":"GRCh38","end":151070325,"species":"homo_sapiens","db_type":"core","object_type":"Exon"},{"start":151068319,"strand":-1,"version":1,"id":"ENSE00001748787","seq_region_name":"1","species":"homo_sapiens","end":151068380,"db_type":"core","object_type":"Exon","assembly_name":"GRCh38"},{"strand":-1,"start":151067132,"seq_region_name":"1","version":1,"id":"ENSE00001695783","end":151067334,"species":"homo_sapiens","object_type":"Exon","db_type":"core","assembly_name":"GRCh38"},{"assembly_name":"GRCh38","end":151059773,"species":"homo_sapiens","db_type":"core","object_type":"Exon","seq_region_name":"1","version":1,"id":"ENSE00001871667","strand":-1,"start":151059479},{"assembly_name":"GRCh38","end":151055993,"species":"homo_sapiens","object_type":"Exon","db_type":"core","seq_region_name":"1","id":"ENSE00003519252","version":1,"start":151055677,"strand":-1},{"species":"homo_sapiens","end":151055126,"object_type":"Exon","db_type":"core","assembly_name":"GRCh38","strand":-1,"start":151055016,"version":1,"id":"ENSE00003604581","seq_region_name":"1"},{"object_type":"Exon","db_type":"core","species":"homo_sapiens","end":151054321,"assembly_name":"GRCh38","start":151054231,"strand":-1,"version":1,"id":"ENSE00003684539","seq_region_name":"1"},{"assembly_name":"GRCh38","end":151053327,"species":"homo_sapiens","object_type":"Exon","db_type":"core","seq_region_name":"1","id":"ENSE00001847145","version":1,"start":151050971,"strand":-1}],"is_canonical":0,"Parent":"ENSG00000197622","display_name":"CDC42SE1-202","species":"homo_sapiens","end":151070325,"biotype":"protein_coding","length":3498},"ENST00000616952":{"strand":-1,"start":10328936,"logic_name":"havana_homo_sapiens","seq_region_name":"21","id":"ENST00000616952","version":1,"source":"havana","db_type":"core","object_type":"Transcript","UTR":[],"assembly_name":"GRCh38","biotype":"lncRNA","length":225,"Parent":"ENSG00000277693","end":10342737,"species":"homo_sapiens","Exon":[{"start":10342616,"strand":-1,"seq_region_name":"21","version":1,"id":"ENSE00003740779","object_type":"Exon","db_type":"core","end":10342737,"species":"homo_sapiens","assembly_name":"GRCh38"},{"seq_region_name":"21","id":"ENSE00003726864","version":1,"start":10328936,"strand":-1,"assembly_name":"GRCh38","object_type":"Exon","db_type":"core","end":10329038,"species":"homo_sapiens"}],"is_canonical":0}}"##
    };
    #[test]
    fn test_transcript_() {
        let json =
            serde_json::from_str::<std::collections::BTreeMap<String, super::Transcript>>(JSON)
                .unwrap();
        for v in json.values() {
            println!("{:?}", v);
        }
    }
}
