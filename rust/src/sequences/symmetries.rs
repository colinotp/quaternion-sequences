use crate::sequences::{williamson::*, equivalence::*};
use std::collections::HashSet;


#[derive(Clone)]
pub enum Symmetry{ // enum for the different types of Quaternion Sequences
    I, II, III, IV
}

#[derive(Clone)]
pub enum RowsumPairing{
    WX, WY, WZ
}

#[derive(Clone)]
pub enum SequenceType{ // enum for the different types of Quadruplets of sequences
    Williamson, WilliamsonType, QuaternionType, ItoType, ExtraTypeI, ExtraTypeII, ExtraTypeIII, Hadamard
}

impl SequenceType {
    // Returns a list of equivalence operations for the given sequence type
    pub fn equivalences(&self) -> Vec<fn(&QuadSeq, SequenceType, bool) -> HashSet<QuadSeq>> {
        match self {
            Self::QuaternionType => vec![equivalent_double_negate, equivalent_uniform_shift, equivalent_double_reorder, equivalent_even_alternated_negation, equivalent_automorphism, equivalent_reverse, equivalent_negate_swap],
            Self::WilliamsonType => vec![equivalent_negate, equivalent_uniform_shift, equivalent_reorder, equivalent_even_alternated_negation, equivalent_automorphism, equivalent_reverse],
            Self::Williamson => vec![equivalent_negate, equivalent_uniform_half_shift, equivalent_reorder, equivalent_even_alternated_negation, equivalent_automorphism],
            Self::Hadamard => vec![equivalent_double_negate, equivalent_disjoint_swaps],
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
            SequenceType::ExtraTypeIII => "et3".to_string(),
            SequenceType::Hadamard => "hm".to_string()
        }
    }
}
