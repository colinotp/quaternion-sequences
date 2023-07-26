use std::collections::HashMap;

use crate::sequences::fourier::iter_over_filtered_couples;

use super::williamson::{SequenceTag, periodic_autocorrelation, cross_correlation};


#[derive(Hash, Eq, PartialEq)]
pub struct MatchData {
    auto_correlation_values : Vec<isize>,
    cross_correlation_values : Vec<isize>,
}


impl MatchData {

    pub fn from(seq1 : &Vec<i8>, seq2 : &Vec<i8>, tags : &(SequenceTag, SequenceTag)) -> MatchData {
        MatchData { auto_correlation_values: compute_auto_correlations(seq1, seq2), cross_correlation_values: compute_cross_correlations(seq1, seq2, tags) }
    }

    pub fn new(auto_correlation_values : Vec<isize>, cross_correlation_values : Vec<isize>) -> MatchData {
        MatchData { auto_correlation_values, cross_correlation_values }
    }

}

pub fn compute_auto_correlations(seq1 : &Vec<i8>, seq2 : &Vec<i8>) -> Vec<isize> {
    let mut res = vec![];
    for offset in 1..seq1.len() {
        res.push(periodic_autocorrelation(seq1, offset) + periodic_autocorrelation(seq2, offset));
    }
    res
}

pub fn compute_cross_correlations(seq1 : &Vec<i8>, seq2 : &Vec<i8>, tags : &(SequenceTag, SequenceTag)) -> Vec<isize> {

    let compute_crossc_with_offset = 
        match tags {
            (SequenceTag::W, _) | (SequenceTag::X, SequenceTag::Y) | (SequenceTag::Y, SequenceTag::Z) | (SequenceTag::Z, SequenceTag::X) => {
                |s1, s2, offset| cross_correlation(s1, s2, offset) - cross_correlation(s2, s1, offset)
            }
            (_, SequenceTag::W) | (SequenceTag::Y, SequenceTag::X) | (SequenceTag::Z, SequenceTag::Y) | (SequenceTag::X, SequenceTag::Z) => {
                |s1, s2, offset| cross_correlation(s2, s1, offset) - cross_correlation(s1, s2, offset)
            }
            _ => {panic!("incorrect tags entered !")}
    };

    let mut res = vec![];
    for offset in 0..seq1.len() {
        res.push(compute_crossc_with_offset(seq1, seq2, offset));
    }
    res
}


pub fn compute_complementary_auto_correlations(seq1 : &Vec<i8>, seq2 : &Vec<i8>) -> Vec<isize> {
    let mut res = vec![];
    
    for offset in 1..seq1.len() {
        res.push(-(periodic_autocorrelation(seq1, offset) + periodic_autocorrelation(seq2, offset)));
    }
    res
}

pub fn compute_complementary_cross_correlations(seq1 : &Vec<i8>, seq2 : &Vec<i8>, tags : &(SequenceTag, SequenceTag)) -> Vec<isize> {

    let compute_crossc_with_offset = 
        match tags {
            (SequenceTag::W, _) | (SequenceTag::X, SequenceTag::Y) | (SequenceTag::Y, SequenceTag::Z) | (SequenceTag::Z, SequenceTag::X) => {
                |s1, s2, offset| cross_correlation(s1, s2, offset) - cross_correlation(s2, s1, offset)
            }
            (_, SequenceTag::W) | (SequenceTag::Y, SequenceTag::X) | (SequenceTag::Z, SequenceTag::Y) | (SequenceTag::X, SequenceTag::Z) => {
                |s1, s2, offset| cross_correlation(s2, s1, offset) - cross_correlation(s1, s2, offset)
            }
            _ => {panic!("incorrect tags entered :{:?}, {:?}", tags.0, tags.1)}
    };

    let mut res = vec![];
    for offset in 0..seq1.len() {
        res.push(-compute_crossc_with_offset(seq1, seq2, offset));
    }
    res
}




pub fn verify_cross_correlation(sequences : &[&Vec<i8>; 4], tags : &Vec<SequenceTag>) -> bool {

    let seqx : &Vec<i8> = {
        let mut res = None;
        for i in 0..4 {
            if tags[i] == SequenceTag::X { 
                res =  Some(sequences[i]);
            } 
        }
        res.expect("Missing tag X !")
    };
    let seqy : &Vec<i8> = {
        let mut res = None;
        for i in 0..4 {
            if tags[i] == SequenceTag::Y { 
                res =  Some(sequences[i]);
            } 
        }
        res.expect("Missing tag Y !")
    };
    let seqz : &Vec<i8> = {
        let mut res = None;
        for i in 0..4 {
            if tags[i] == SequenceTag::Z { 
                res =  Some(sequences[i]);
            } 
        }
        res.expect("Missing tag Z !")
    };
    let seqw : &Vec<i8> = {
        let mut res = None;
        for i in 0..4 {
            if tags[i] == SequenceTag::W { 
                res =  Some(sequences[i]);
            } 
        }
        res.expect("Missing tag W !")
    };


    let (seq_for_cond1, seq_for_cond2) = match (tags[0].clone(), tags[1].clone()) {
        (SequenceTag::X, SequenceTag::Y) | (SequenceTag::Y, SequenceTag::X) | (SequenceTag::W, SequenceTag::Z) | (SequenceTag::Z, SequenceTag::W) => {
            ((seqw, seqx, seqy, seqz,), (seqw, seqy, seqz, seqx))
        }
        (SequenceTag::X, SequenceTag::Z) | (SequenceTag::Z, SequenceTag::X) | (SequenceTag::W, SequenceTag::Y) | (SequenceTag::Y, SequenceTag::W) => {
            ((seqw, seqx, seqy, seqz,), (seqw, seqz, seqx, seqy))
        }
        (SequenceTag::X, SequenceTag::W) | (SequenceTag::W, SequenceTag::X) | (SequenceTag::Y, SequenceTag::Z) | (SequenceTag::Z, SequenceTag::Y) => {
            ((seqw, seqy, seqz, seqx,), (seqw, seqz, seqx, seqy))
        }
        _ => {panic!("Incompatible tags !")}
    };


    for offset in 0..sequences[0].len() {
        
        if cross_correlation(seq_for_cond1.0, seq_for_cond1.1, offset) - cross_correlation(seq_for_cond1.1, seq_for_cond1.0, offset) +
           cross_correlation(seq_for_cond1.2, seq_for_cond1.3, offset) - cross_correlation(seq_for_cond1.3, seq_for_cond1.2, offset) != 0 ||
           cross_correlation(seq_for_cond2.0, seq_for_cond2.1, offset) - cross_correlation(seq_for_cond2.1, seq_for_cond2.0, offset) +
           cross_correlation(seq_for_cond2.2, seq_for_cond2.3, offset) - cross_correlation(seq_for_cond2.3, seq_for_cond2.2, offset) != 0 {
            return false;
        }
    }

    true
}





pub fn generate_matching_table<'a, 'b>(sequences1 : &'a Vec<Vec<i8>>, sequences2 : &'a Vec<Vec<i8>>, tags : &'b (SequenceTag, SequenceTag), p : usize) -> HashMap<MatchData, Vec<(&'a Vec<i8>, &'a Vec<i8>)>> {

    let mut match_table : HashMap<MatchData, Vec<(&Vec<i8>, &Vec<i8>)>> = HashMap::new();

    for (seq1, seq2) in iter_over_filtered_couples(&sequences1, &sequences2, 4.*p as f64) {
        // We iterate over the couples of sequences, but we filter out some with the dft checks

        let match_data = MatchData::from(seq1, seq2, tags);

        match match_table.get_mut(&match_data) {
            Some(list) => {list.push((seq1, seq2))}
            None => {match_table.insert(match_data, vec![]);}
        }
    }

    match_table
}