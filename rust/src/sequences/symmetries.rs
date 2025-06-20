use crate::sequences::{williamson::*, equivalence::*};


#[derive(Clone)]
pub enum Symmetry{ // enum for the different types of Quaternion Sequences
    I, II, III, IV
}

#[derive(Clone)]
pub enum RowsumPairing{
    XY, XZ, XW
}

#[derive(Clone)]
pub enum SequenceType{ // enum for the different types of Quadruplets of sequences
    Williamson, WilliamsonType, QuaternionType, ItoType, ExtraTypeI, ExtraTypeII, ExtraTypeIII
}

impl SequenceType {
    // Returns a list of equivalence operations for the given sequence type
    pub fn equivalences(&self) -> Vec<fn(&QuadSeq) -> Vec<QuadSeq>> {
        match self {
            Self::QuaternionType => vec![equivalent_double_negate, equivalent_uniform_shift, equivalent_double_reorder, equivalent_alternated_negation, equivalent_automorphism, equivalent_reverse],
            Self::WilliamsonType => vec![equivalent_negate, equivalent_uniform_shift, equivalent_reorder, equivalent_alternated_negation, equivalent_automorphism, equivalent_reverse],
            Self::Williamson => vec![equivalent_negate, equivalent_uniform_half_shift, equivalent_reorder, equivalent_alternated_negation, equivalent_automorphism],
            _ => vec![]
        }
    }
}

impl ToString for SequenceType {
    fn to_string(&self) -> String {
        match self {
            SequenceType::Williamson => "ws".to_string(),
            SequenceType::WilliamsonType => "wts".to_string(),
            SequenceType::QuaternionType => "qts".to_string(),
            SequenceType::ItoType => "its".to_string(),
            SequenceType::ExtraTypeI => "et1".to_string(),
            SequenceType::ExtraTypeII => "et2".to_string(),
            SequenceType::ExtraTypeIII => "et3".to_string()
        }
    }
}
