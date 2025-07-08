
#[cfg(test)]
mod tests {
    use num_complex::Complex64;

    use crate::sequences::{fourier::{dft_sequence, inverse_dft}, rowsum::generate_sequences_with_rowsum};


    #[test]
    fn test_dft(){
        let eps = 0.01;
        
        let seq = vec![1];
        let _ = dft_sequence(&seq);
        
        let seq = vec![1,-1,1];
        let dft = dft_sequence(&seq);
        assert_eq!(dft[0], Complex64::new(1.,0.));
        assert!((dft[1] - Complex64::new(1.,1.732)).norm() <= eps);

        let seq = vec![1,-1,1,-1,1,-1,1,-1];
        let dft = dft_sequence(&seq);
        assert!((dft[0] - Complex64::new(0.,0.)).norm() <= eps);
        assert!((dft[1] - Complex64::new(0.,0.)).norm() <= eps);
        assert!((dft[2] - Complex64::new(0.,0.)).norm() <= eps);
        assert!((dft[3] - Complex64::new(0.,0.)).norm() <= eps);
        assert!((dft[4] - Complex64::new(8.,0.)).norm() <= eps);

        let seq = vec![1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
        let _ = dft_sequence(&seq);
    }

    #[test]
    fn test_dft_inv() {
        let seq : Vec<i8> = vec![-1,1,-1,-1];

        let dft = dft_sequence(&seq);
        let inv : Vec<i8> = inverse_dft(&dft, seq.len()).iter().map(|elm| elm.round() as i8).collect();

        assert_eq!(seq, inv);
    }

    #[test]
    fn test_filter(){
        
        let r = 2;
        let p = 5;
        let seqs = generate_sequences_with_rowsum(r, p);

        for seq in seqs{
            println!("===new seq===");
            for elm in dft_sequence(&seq){
                println!("{:?}", elm.norm_sqr());
            }
        }

    }
}