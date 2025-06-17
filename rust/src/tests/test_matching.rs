

#[cfg(test)]
mod tests {
    use crate::{find::{find_with_rowsum::sort, find_write::join_pairs}, sequences::{matching::{compute_auto_correlations, compute_complementary_auto_correlations, compute_complementary_cross_correlations, compute_cross_correlations, verify_cross_correlation}, symmetries::SequenceType}};
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

        assert_eq!(compute_auto_correlations(&seq_x, &seq_y), compute_complementary_auto_correlations(&seq_z, &seq_w));
        assert_eq!(compute_auto_correlations(&seq_w, &seq_y), compute_complementary_auto_correlations(&seq_z, &seq_x));

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
    fn test_join() {
        join_pairs(7, SequenceType::QuaternionType);
    }

}