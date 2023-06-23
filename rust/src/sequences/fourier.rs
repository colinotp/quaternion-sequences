
use fftw::array::AlignedVec;
use fftw::plan::*;
use fftw::types::*;
use num_complex::Complex;

pub fn dft_sequence(seq : &Vec<i8>) -> Vec<Complex<f64>>{
    // returns the dft of the sequence
    let n = seq.len();
    let mut plan: R2CPlan64 = R2CPlan::aligned(&[n], Flag::MEASURE).unwrap();
    let mut a = AlignedVec::new(n);
    let mut b = AlignedVec::new(n/2+1);
    for i in 0..n {
        a[i] = seq[i] as f64;
    }
    plan.r2c(&mut a, &mut b).unwrap();
    
    b.to_vec()
}

pub fn iter_over_filtered_dft<'a>(sequences : &'a Vec<Vec<i8>>, bound : &'a f64) -> impl std::iter::Iterator<Item = &'a Vec<i8>> {
    sequences.iter()
        .filter(|seq| {
            let dft = dft_sequence(seq);
            for elm in dft {
                if elm.norm_sqr() > *bound {return false;}
            }
            true
        })
}