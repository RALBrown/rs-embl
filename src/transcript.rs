use itertools::Itertools;
use regex::Regex;
use serde::{Deserialize, Serialize};

use crate::{
    sequence::{CdnaSequence, GenomicSequence},
    Client,
};

const LAST_EJC_REGEX: &str = r".+([A-Z][a-z]+[A-Z]+)$";

/**

*/
#[derive(Debug, Clone, Serialize, Deserialize, Hash, PartialEq, Eq)]
pub struct Transcript {
    pub id: String,
    #[serde(default)]
    pub display_name: String,
    pub start: u32,
    pub end: u32,
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
        upstream = &seq.seq[..(start - transcript.start) as usize];
        downstream = &seq.seq[(end - transcript.start + 1) as usize..];
    } else {
        upstream = &seq.seq[..(transcript.end - end) as usize];
        downstream = &seq.seq[(transcript.end - start + 1) as usize..];
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
    const TTR_V30M_DEL_PROTEIN: &str = "MASHRLLLLCLAGLVFVSEAGPTGTGESKCPLMVKVLDAVRGSPAINVACMCSERLLMTPGSHLPLGKPVSLESCMGSQLRRNL*";
    const TTR_V30M_INS_PROTEIN: &str = "MASHRLLLLCLAGLVFVSEAGPTGTGESKCPLMVKVLDAVRGSPAINVAGACVQKGC*";
    #[test]
    fn test_snp(){
        let transcript = serde_json::from_str::<super::Transcript>(TTR_201_JSON).unwrap();
        let genomic_seq: GenomicSequence = GenomicSequence{
            query: "".to_owned(),
            id: "".to_owned(),
            desc: None,
            seq: TTR_GENOME_SEQ.to_owned(),
        };
        const END: u32 = 31592974;
        const START: u32 = 31592974;
        const VARIANT_ALLELE: &str = "A";
        let consequences = super::make_consequences(
            &genomic_seq,
            &transcript, 
            START, 
            END, 
            VARIANT_ALLELE);
        let super::Consequences::Coding{
            edited_genomic_sequence,
            edited_protein_sequence,
            unedited_protein_sequence,} = consequences else {panic!()};
        const V30M_TTR: &str = "MASHRLLLLCLAGLVFVSEAGPTGTGESKCPLMVKVLDAVRGSPAINVAMHVFRKAADDTWEPFASGKTSESGELHGLTTEEEFVEGIYKVEIDTKSYWKALGISPFHEHAEVVFTANDSGPRRYTIAALLSPYSYSTTAVVTNPKE*";
        assert_eq!(&edited_protein_sequence.protein_sequence, V30M_TTR);
    }
    #[test]
    fn test_del(){
        let transcript = serde_json::from_str::<super::Transcript>(TTR_201_JSON).unwrap();
        let genomic_seq: GenomicSequence = GenomicSequence{
            query: "".to_owned(),
            id: "".to_owned(),
            desc: None,
            seq: TTR_GENOME_SEQ.to_owned(),
        };
        const END: u32 = 31592974;
        const START: u32 = 31592974;
        const VARIANT_ALLELE: &str = "-";
        let consequences = super::make_consequences(
            &genomic_seq,
            &transcript, 
            START, 
            END, 
            VARIANT_ALLELE);
        let super::Consequences::Coding{
            edited_genomic_sequence,
            edited_protein_sequence,
            unedited_protein_sequence,} = consequences else {panic!()};
        assert_eq!(&edited_protein_sequence.protein_sequence, TTR_V30M_DEL_PROTEIN);
    }
    #[test]
    fn test_ins(){
        let transcript = serde_json::from_str::<super::Transcript>(TTR_201_JSON).unwrap();
        let genomic_seq: GenomicSequence = GenomicSequence{
            query: "".to_owned(),
            id: "".to_owned(),
            desc: None,
            seq: TTR_GENOME_SEQ.to_owned(),
        };
        const END: u32 = 31592974;
        const START: u32 = 31592975;
        const VARIANT_ALLELE: &str = "G";
        let consequences = super::make_consequences(
            &genomic_seq,
            &transcript, 
            START, 
            END, 
            VARIANT_ALLELE);
        let super::Consequences::Coding{
            edited_genomic_sequence,
            edited_protein_sequence,
            unedited_protein_sequence,} = consequences else {panic!()};
        assert_eq!(&edited_protein_sequence.protein_sequence, TTR_V30M_INS_PROTEIN);
    }
}