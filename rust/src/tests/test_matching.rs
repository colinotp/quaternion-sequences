

#[cfg(test)]
mod tests {
    use crate::{find::{find_with_rowsum::sort, find_write::join_pairs}, sequences::{fourier::dft_sequence, matching::{compute_auto_correlation, compute_auto_correlation_dft, compute_auto_correlation_pair, compute_auto_correlation_pair_dft, compute_complementary_auto_correlations, compute_complementary_cross_correlations, compute_cross_correlations, compute_cross_correlations_dft_aux, verify_cross_correlation}, symmetries::SequenceType, williamson::cross_correlation}};
    use crate::sequences::williamson::SequenceTag;


    #[test]
    fn test_sort() {

        let quad = (5,3,1,2);
        let (sorted, indices) = sort(&quad);
        assert_eq!(sorted, [5,3,2,1]);
        assert_eq!(indices, [0,1,3,2]);

        let quad = (2,6,4,7);
        let (sorted, indices) = sort(&quad);
        assert_eq!(sorted, [7,6,4,2]);
        assert_eq!(indices, [3,1,2,0]);

        let quad = (2,6,4,1);
        let (sorted, indices) = sort(&quad);
        assert_eq!(sorted, [6,4,2,1]);
        assert_eq!(indices, [1,2,0,3]);
    }

    #[test]
    fn test_verify_crossc() {
        let seq_x = vec![1,-1,-1,-1,1,1,-1,1,-1,1];
        let seq_y = vec![-1,1,1,1,-1,1,-1,-1,-1,1];
        let seq_z = vec![1,1,1,1,1,-1,-1,1,-1,-1];
        let seq_w = vec![1,1,-1,1,1,1,1,-1,1,1];

        let sequences = &[&seq_x, &seq_y, &seq_z, &seq_w];
        let tags = &vec![SequenceTag::X, SequenceTag::Y, SequenceTag::Z, SequenceTag::W];
        assert!(verify_cross_correlation(sequences, tags));
        
        let sequences = &[&seq_x, &seq_z, &seq_w, &seq_y];
        let tags = &vec![SequenceTag::X, SequenceTag::Z, SequenceTag::W, SequenceTag::Y];
        assert!(verify_cross_correlation(sequences, tags));

        let sequences = &[&seq_z, &seq_x, &seq_w, &seq_y];
        let tags = &vec![SequenceTag::Z, SequenceTag::X, SequenceTag::W, SequenceTag::Y];
        assert!(verify_cross_correlation(sequences, tags));
        
        let seq_w = vec![1,1,-1,1,1,1,1,-1,-1,1];
        
        let sequences = &[&seq_z, &seq_x, &seq_w, &seq_y];
        let tags = &vec![SequenceTag::Z, SequenceTag::X, SequenceTag::W, SequenceTag::Y];
        assert!(!verify_cross_correlation(sequences, tags));
    }

    #[test]
    fn test_autocorrelation() {
        let seq_x = vec![1,-1,-1,-1,1,1,-1,1,-1,1];
        let seq_y = vec![-1,1,1,1,-1,1,-1,-1,-1,1];
        let seq_z = vec![1,1,1,1,1,-1,-1,1,-1,-1];
        let seq_w = vec![1,1,-1,1,1,1,1,-1,1,1];

        assert_eq!(compute_auto_correlation_pair(&seq_x, &seq_y), compute_complementary_auto_correlations(&seq_z, &seq_w));
        assert_eq!(compute_auto_correlation_pair(&seq_w, &seq_y), compute_complementary_auto_correlations(&seq_z, &seq_x));

    }

    #[test]
    fn test_autocorrelation_dft() {
        // Test sequences
        let seq1 = vec![1,-1,-1,-1,1,1,-1,1,-1,1];
        let seq2 = vec![-1,1,1,1,-1,1,-1,-1,-1,1];
        let seq3 = vec![1,1,1,1,1,-1,-1,1,-1,-1];
        let seq4 = vec![1,1,-1,1,1,1,1,-1,1,1];

        // Compute PSD vectors from DFT
        let psd1 : Vec<f64> = dft_sequence(&seq1).iter().map(|elm| elm.norm_sqr()).collect();
        let psd2 : Vec<f64> = dft_sequence(&seq2).iter().map(|elm| elm.norm_sqr()).collect();
        let psd3 : Vec<f64> = dft_sequence(&seq3).iter().map(|elm| elm.norm_sqr()).collect();
        let psd4 : Vec<f64> = dft_sequence(&seq4).iter().map(|elm| elm.norm_sqr()).collect();

        // Computing and filtering out second half of autocorrelation values to make output vector align with normal autocorrelation
        let auto_dft1 : Vec<isize> = compute_auto_correlation_dft(&psd1, seq1.len()).iter().enumerate().filter_map(|(i, &elm)| if i < seq1.len() / 2 {Some(elm)} else {None}).collect();
        let auto_dft2 : Vec<isize> = compute_auto_correlation_dft(&psd2, seq1.len()).iter().enumerate().filter_map(|(i, &elm)| if i < seq1.len() / 2 {Some(elm)} else {None}).collect();
        let auto_dft3 : Vec<isize> = compute_auto_correlation_dft(&psd3, seq1.len()).iter().enumerate().filter_map(|(i, &elm)| if i < seq1.len() / 2 {Some(elm)} else {None}).collect();
        let auto_dft4 : Vec<isize> = compute_auto_correlation_dft(&psd4, seq1.len()).iter().enumerate().filter_map(|(i, &elm)| if i < seq1.len() / 2 {Some(elm)} else {None}).collect();

        // Verify autocorrelation functions are equivalent
        assert_eq!(compute_auto_correlation(&seq1), auto_dft1);
        assert_eq!(compute_auto_correlation(&seq2), auto_dft2);
        assert_eq!(compute_auto_correlation(&seq3), auto_dft3);
        assert_eq!(compute_auto_correlation(&seq4), auto_dft4);

        // Verify pair autocorrelation functions are equivalent
        assert_eq!(compute_auto_correlation_pair(&seq1, &seq2), compute_auto_correlation_pair_dft(&psd1, seq1.len(), &psd2, seq2.len()));
        assert_eq!(compute_auto_correlation_pair(&seq1, &seq3), compute_auto_correlation_pair_dft(&psd1, seq1.len(), &psd3, seq3.len()));
        assert_eq!(compute_auto_correlation_pair(&seq1, &seq4), compute_auto_correlation_pair_dft(&psd1, seq1.len(), &psd4, seq4.len()));
        assert_eq!(compute_auto_correlation_pair(&seq2, &seq3), compute_auto_correlation_pair_dft(&psd2, seq2.len(), &psd3, seq3.len()));
        assert_eq!(compute_auto_correlation_pair(&seq2, &seq4), compute_auto_correlation_pair_dft(&psd2, seq2.len(), &psd4, seq4.len()));
        assert_eq!(compute_auto_correlation_pair(&seq3, &seq4), compute_auto_correlation_pair_dft(&psd3, seq3.len(), &psd4, seq4.len()));
    }

    #[test]
    fn test_crosscorrelation() {
        let seq_x = vec![1,-1,-1,-1,1,1,-1,1,-1,1];
        let seq_y = vec![-1,1,1,1,-1,1,-1,-1,-1,1];
        let seq_z = vec![1,1,1,1,1,-1,-1,1,-1,-1];
        let seq_w = vec![1,1,-1,1,1,1,1,-1,1,1];

        assert_eq!(compute_cross_correlations(&seq_x, &seq_y, &(SequenceTag::X, SequenceTag::Y)), compute_complementary_cross_correlations(&seq_z, &seq_w, &(SequenceTag::Z, SequenceTag::W)));
        assert_eq!(compute_cross_correlations(&seq_w, &seq_y, &(SequenceTag::W, SequenceTag::Y)), compute_complementary_cross_correlations(&seq_z, &seq_x, &(SequenceTag::Z, SequenceTag::X)));

    }

    #[test]
    fn test_crosscorrelation_dft() {
        let seq_x = vec![1,-1,-1,-1,1,1,-1,1,-1,1];
        let seq_y = vec![-1,1,1,1,-1,1,-1,-1,-1,1];
        let seq_z = vec![1,1,1,1,1,-1,-1,1,-1,-1];
        let seq_w = vec![1,1,-1,1,1,1,1,-1,1,1];

        let dft_xy = compute_cross_correlations_dft_aux(&seq_x, &seq_y);
        let mut xy_vec = vec![];
        for offset in 0..seq_x.len() {
            xy_vec.push(cross_correlation(&seq_x, &seq_y, offset));
        }

        assert_eq!(dft_xy, xy_vec);

        for offset in 1..seq_x.len() {
            assert_eq!(dft_xy[offset], cross_correlation(&seq_x, &seq_y, offset));
        }

        let dft_xz = compute_cross_correlations_dft_aux(&seq_x, &seq_z);
        for offset in 1..seq_x.len() {
            assert_eq!(dft_xz[offset], cross_correlation(&seq_x, &seq_z, offset));
        }

        let dft_xw = compute_cross_correlations_dft_aux(&seq_x, &seq_w);
        for offset in 1..seq_x.len() {
            assert_eq!(dft_xw[offset], cross_correlation(&seq_x, &seq_w, offset));
        }
        
        let dft_yz = compute_cross_correlations_dft_aux(&seq_y, &seq_z);
        for offset in 1..seq_x.len() {
            assert_eq!(dft_yz[offset], cross_correlation(&seq_y, &seq_z, offset));
        }
        
        let dft_yw = compute_cross_correlations_dft_aux(&seq_y, &seq_w);
        for offset in 1..seq_x.len() {
            assert_eq!(dft_yw[offset], cross_correlation(&seq_y, &seq_w, offset));
        }
        
        let dft_zw = compute_cross_correlations_dft_aux(&seq_z, &seq_w);
        for offset in 1..seq_x.len() {
            assert_eq!(dft_zw[offset], cross_correlation(&seq_z, &seq_w, offset));
        }
    }


    #[test]
    fn test_join() {
        join_pairs(7, SequenceType::QuaternionType);
    }

}