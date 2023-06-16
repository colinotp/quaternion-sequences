
#[cfg(test)]
mod tests {
    use num_complex::Complex64;

    use crate::sequences::fourier::dft_sequence;


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
}