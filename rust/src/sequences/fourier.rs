
use fftw::array::AlignedVec;
use fftw::plan::*;
use fftw::types::*;
use itertools::iproduct;
use num_complex::Complex;
use std::sync::{OnceLock, Mutex};

static R2C_PLAN : OnceLock<Mutex<R2CPlan64>> = OnceLock::new();
static C2R_PLAN : OnceLock<Mutex<C2RPlan64>> = OnceLock::new();

pub fn dft_sequence(seq : &Vec<i8>) -> Vec<Complex<f64>>{
    // returns the dft of the sequence
    let n = seq.len();
    
    let mutex = R2C_PLAN.get_or_init(|| Mutex::new(R2CPlan::aligned(&[n], Flag::MEASURE).expect("Failed to create Mutex for FFTW plan")));
    let mut plan = mutex.lock().unwrap();

    let mut a = AlignedVec::new(n);
    let mut b = AlignedVec::new(n/2+1);
    for i in 0..n {
        a[i] = seq[i] as f64;
    }
    plan.r2c(&mut a, &mut b).unwrap();
    
    b.to_vec()
}

pub fn inverse_dft(freq: &Vec<Complex<f64>>, n: usize) -> Vec<f64> {
    let mutex = C2R_PLAN.get_or_init(|| Mutex::new(C2RPlan::aligned(&[n], Flag::MEASURE).expect("Failed to create Mutex for FFTW plan")));
    let mut plan = mutex.lock().unwrap();
    
    //let plan = C2R_PLAN.get_or_init(|| C2RPlan::aligned(&[n], Flag::MEASURE).expect("Failed to create FFT plan"));
    let mut a = AlignedVec::new(n / 2 + 1);
    let mut b = AlignedVec::new(n);
    for i in 0..a.len() {
        a[i] = freq[i];
    }
    plan.c2r(&mut a, &mut b).unwrap();

    let norm: f64 = n as f64;
    b.iter().map(|x| x / norm).collect()
}


pub fn iter_over_filtered_dft<'a>(sequences : &'a Vec<Vec<i8>>, bound : f64) -> impl std::iter::Iterator<Item = &'a Vec<i8>> {
    let bound_sqr = bound * bound;
    sequences.iter()
        .filter(move |seq| {
            let dft = dft_sequence(seq);
            for elm in dft {
                if elm.norm_sqr() > bound_sqr {return false;}
            }
            true
        })
}


pub fn iter_over_filtered_couples<'a>(sequences1 : &'a Vec<Vec<i8>>, sequences2 : &'a Vec<Vec<i8>>, bound : f64) -> impl std::iter::Iterator<Item = (&'a Vec<i8>, &'a Vec<i8>)> {
    let couples = iproduct!(sequences1, sequences2);
    couples
        .filter(move |(seq1, seq2)| {
            let dft1 = dft_sequence(seq1);
            let dft2 = dft_sequence(seq2);
            for (elm1, elm2) in dft1.iter().zip(dft2.iter()) {
                if elm1.norm_sqr() + elm2.norm_sqr() > bound {return false;}
            }
            true
        })
}



pub fn iter_over_enumerate_filtered_couples<'a>(sequences1 : &'a Vec<Vec<i8>>, sequences2 : &'a Vec<Vec<i8>>, bound : f64) -> impl std::iter::Iterator<Item = ((usize, &'a Vec<i8>), (usize, &'a Vec<i8>))> {
    let couples = iproduct!(sequences1.iter().enumerate(), sequences2.iter().enumerate());      // Iterator for every unique combination of sequences
    couples
        .filter(move |((_, seq1), (_, seq2))| {
            let dft1 = dft_sequence(seq1);
            let dft2 = dft_sequence(seq2);
            for (elm1, elm2) in dft1.iter().zip(dft2.iter()) {                
                if elm1.norm_sqr() + elm2.norm_sqr() > bound {return false;}
            }
            true
        })
}

pub struct SequencePairData<'a> {
    pub seq_enum1: (usize, &'a Vec<i8>),
    pub seq_enum2: (usize, &'a Vec<i8>),
    pub dft1: Vec<Complex<f64>>,
    pub norm1: Vec<f64>,
    pub dft2: Vec<Complex<f64>>,
    pub norm2: Vec<f64>,
}

pub fn iter_over_enumerate_filtered_couples_psds<'a>(sequences1: &'a Vec<Vec<i8>>, sequences2: &'a Vec<Vec<i8>>, bound: f64,) -> impl Iterator<Item = SequencePairData<'a>> {
    let couples = iproduct!(sequences1.iter().enumerate(), sequences2.iter().enumerate());

    couples.filter_map(move |(seq_enum1, seq_enum2)| {
        let dft1 = dft_sequence(seq_enum1.1);
        let dft2 = dft_sequence(seq_enum2.1);

        let norm1: Vec<f64> = dft1.iter().map(|elm| elm.norm_sqr().round()).collect();
        let norm2: Vec<f64> = dft2.iter().map(|elm| elm.norm_sqr().round()).collect();
        
        let bound_exceeded = norm1.iter().enumerate().zip(norm2.iter().enumerate()).any(|((i1, n1), (i2, n2))| i1 != 0 && i2 != 0 && n1 + n2 > bound);

        if bound_exceeded {
            None
        } else {
            Some(SequencePairData {
                seq_enum1,
                seq_enum2,
                dft1,
                norm1,
                dft2,
                norm2
            })
        }
    })
}
