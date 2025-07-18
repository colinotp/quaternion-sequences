use std::collections::HashMap;

use num_complex::Complex;

use crate::sequences::{fourier::{dft_sequence, inverse_dft, iter_over_filtered_couples}, sequence::{transpose, seq_multiply_pointwise_complex}};

use super::williamson::{SequenceTag, periodic_autocorrelation, cross_correlation};


#[derive(Hash, Eq, PartialEq)]
pub struct MatchData {
    auto_correlation_values : Vec<isize>,
    cross_correlation_values : Vec<isize>,
}


impl MatchData {

    pub fn from(seq1 : &Vec<i8>, seq2 : &Vec<i8>, tags : &(SequenceTag, SequenceTag)) -> MatchData {
        MatchData { auto_correlation_values: compute_auto_correlation_pair(seq1, seq2), cross_correlation_values: compute_cross_correlations(seq1, seq2, tags) }
    }

    pub fn new(auto_correlation_values : Vec<isize>, cross_correlation_values : Vec<isize>) -> MatchData {
        MatchData { auto_correlation_values, cross_correlation_values }
    }

}

// Compute vector of autocorrelation values for a given sequence. Omits first entry
// No longer used in favour of compute_auto_correlation_psd
pub fn compute_auto_correlation(seq : &Vec<i8>) -> Vec<isize> {
    let mut res = vec![];
    for offset in 1..=(seq.len() / 2) {
        res.push(periodic_autocorrelation(seq, offset));
    }

    res
}

// Compute vector for sum of autocorrelations of two sequences
pub fn compute_auto_correlation_pair(seq1 : &Vec<i8>, seq2 : &Vec<i8>) -> Vec<isize> {
    let auto1 = compute_auto_correlation(&seq1);
    let auto2 = compute_auto_correlation(&seq2);
    
    let mut res = vec![];
    for offset in 0..(seq1.len() / 2) {
        res.push(auto1[offset] + auto2[offset]);
    }
    res
}

pub fn compute_auto_correlation_dft(psd_vec : &Vec<f64>, len : usize) -> Vec<isize> {
    let complex_psd : Vec<Complex<f64>> = psd_vec.iter().map(|&elm| Complex::new(elm as f64, 0.0)).collect();
    let mut res : Vec<isize> = inverse_dft(&complex_psd, len).into_iter().map(|elm| elm.round() as isize).collect();
    res.remove(0);

    res
}

pub fn compute_auto_correlation_pair_dft(psd_vec1 : &Vec<f64>, len1 : usize, psd_vec2 : &Vec<f64>, len2 : usize) -> Vec<isize> {
    let auto1 = compute_auto_correlation_dft(&psd_vec1, len1);
    let auto2 = compute_auto_correlation_dft(&psd_vec2, len2);

    // Happens with sequences of length 1
    if auto1.len() == 0 {
        return vec![];
    }

    let mut res = vec![];

    for offset in 0..=(auto1.len() / 2) {
        res.push(auto1[offset] + auto2[offset]);
    }
    res
}

pub fn compute_cross_correlations(seq1 : &Vec<i8>, seq2 : &Vec<i8>, tags : &(SequenceTag, SequenceTag)) -> Vec<isize> {

    let compute_crossc_with_offset = 
        match tags {
            (SequenceTag::Z, _) | (SequenceTag::W, SequenceTag::X) | (SequenceTag::X, SequenceTag::Y) | (SequenceTag::Y, SequenceTag::W) => {
                |s1, s2, offset| cross_correlation(s1, s2, offset) - cross_correlation(s2, s1, offset)
            }
            (_, SequenceTag::Z) | (SequenceTag::X, SequenceTag::W) | (SequenceTag::Y, SequenceTag::X) | (SequenceTag::W, SequenceTag::Y) => {
                |s1, s2, offset| cross_correlation(s2, s1, offset) - cross_correlation(s1, s2, offset)
            }
            _ => {panic!("incorrect tags entered !")}
    };

    let mut res = vec![];
    for offset in 1..=(seq1.len() / 2) {
        res.push(compute_crossc_with_offset(seq1, seq2, offset));
    }
    res
}

// Computes crosscorrelation sums on a given side of the equation, using DFT
pub fn compute_cross_correlations_dft(seq1 : &Vec<i8>, seq2 : &Vec<i8>, tags : &(SequenceTag, SequenceTag)) -> Vec<isize> {
    let mut mut_seq1 = seq1.clone();
    let mut mut_seq2 = seq2.clone();
    
    let crossc1 = compute_cross_correlations_dft_aux(&mut mut_seq1, &mut mut_seq2);
    let crossc2 = compute_cross_correlations_dft_aux(&mut mut_seq2, &mut mut_seq1);

    let cross_at_offset : Box<dyn Fn(usize) -> isize> = match tags {
        (SequenceTag::Z, _) | (SequenceTag::W, SequenceTag::X) | (SequenceTag::X, SequenceTag::Y) | (SequenceTag::Y, SequenceTag::W) => {
            Box::new(|offset| crossc2[offset] - crossc1[offset])
        }
        (_, SequenceTag::Z) | (SequenceTag::X, SequenceTag::W) | (SequenceTag::Y, SequenceTag::X) | (SequenceTag::W, SequenceTag::Y) => {
            Box::new(|offset| crossc1[offset] - crossc2[offset])
        }
        _ => {panic!("incorrect tags entered !")}
    };

    
    let mut res : Vec<isize> = vec![];
    for offset in 1..=(seq1.len() / 2) {
        res.push(cross_at_offset(offset));
    }

    res
}

// Computes crosscorrelation vector for seq1, seq2 via DFT
pub fn compute_cross_correlations_dft_aux(seq1 : &mut Vec<i8>, seq2 : &mut Vec<i8>) -> Vec<isize> {
    transpose(seq1);

    let dft1 = dft_sequence(&seq1);
    let dft2 = dft_sequence(&seq2);

    let product : Vec<Complex<f64>> = seq_multiply_pointwise_complex(&dft1, &dft2);

    inverse_dft(&product, seq1.len()).into_iter().map(|elm| elm.round() as isize).collect()
}


pub fn compute_cross_psd_pair(dft1 : Vec<Complex<f64>>, dft2 : Vec<Complex<f64>>, tags : &(SequenceTag, SequenceTag), len : usize) -> Vec<(isize, isize)> {
    let cross_psd1 = compute_cross_psd(&dft1, &dft2);
    let cross_psd2 = compute_cross_psd(&dft2, &dft1);

    let cross_at_offset : Box<dyn Fn(usize) -> Complex<isize>> = match tags {
        (SequenceTag::Z, _) | (SequenceTag::W, SequenceTag::X) | (SequenceTag::X, SequenceTag::Y) | (SequenceTag::Y, SequenceTag::W) => {
            Box::new(|offset| cross_psd1[offset] - cross_psd2[offset])
        }
        (_, SequenceTag::Z) | (SequenceTag::X, SequenceTag::W) | (SequenceTag::Y, SequenceTag::X) | (SequenceTag::W, SequenceTag::Y) => {
            Box::new(|offset| cross_psd2[offset] - cross_psd1[offset])
        }
        _ => {panic!("incorrect tags entered !")}
    };

    let mut res : Vec<(isize, isize)> = vec![];
    for offset in 1..=(len/2) {
        let cpsd = cross_at_offset(offset);
        res.push((cpsd.re, cpsd.im));
    }
    
    res
}


// Computes cross-power-spectral-density vector, with real and imaginary components rounded to the nearest integer
pub fn compute_cross_psd(dft1 : &Vec<Complex<f64>>, dft2 : &Vec<Complex<f64>>) -> Vec<Complex<isize>> {
    dft1.into_iter().zip(dft2.into_iter()).map(|(elm1, elm2)| {
        let product = elm1 * elm2.conj();
        Complex::new(product.re.round() as isize, product.im.round() as isize)
    }).collect()
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
            (SequenceTag::Z, _) | (SequenceTag::W, SequenceTag::X) | (SequenceTag::X, SequenceTag::Y) | (SequenceTag::Y, SequenceTag::W) => {
                |s1, s2, offset| cross_correlation(s1, s2, offset) - cross_correlation(s2, s1, offset)
            }
            (_, SequenceTag::Z) | (SequenceTag::X, SequenceTag::W) | (SequenceTag::Y, SequenceTag::X) | (SequenceTag::W, SequenceTag::Y) => {
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

    let seqw : &Vec<i8> = {
        let mut res = None;
        for i in 0..4 {
            if tags[i] == SequenceTag::W { 
                res =  Some(sequences[i]);
            } 
        }
        res.expect("Missing tag W !")
    };
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
    


    let (seq_for_cond1, seq_for_cond2) = match (tags[0].clone(), tags[1].clone()) {
        (SequenceTag::W, SequenceTag::X) | (SequenceTag::X, SequenceTag::W) | (SequenceTag::Z, SequenceTag::Y) | (SequenceTag::Y, SequenceTag::Z) => {
            ((seqz, seqw, seqx, seqy,), (seqz, seqx, seqy, seqw))
        }
        (SequenceTag::W, SequenceTag::Y) | (SequenceTag::Y, SequenceTag::W) | (SequenceTag::Z, SequenceTag::X) | (SequenceTag::X, SequenceTag::Z) => {
            ((seqz, seqw, seqx, seqy,), (seqz, seqy, seqw, seqx))
        }
        (SequenceTag::W, SequenceTag::Z) | (SequenceTag::Z, SequenceTag::W) | (SequenceTag::X, SequenceTag::Y) | (SequenceTag::Y, SequenceTag::X) => {
            ((seqz, seqx, seqy, seqw,), (seqz, seqy, seqw, seqx))
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