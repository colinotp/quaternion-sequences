
#[cfg(test)]
mod tests {

    use crate::sequences::williamson::{QuadSeq, QUADRUPLETS, periodic_autocorrelation, cross_correlation};

    #[test]
    fn test_conversion() {
        let mut will = QuadSeq::new(16);
        for (i, q) in QUADRUPLETS.iter().enumerate(){
            will.set_sequence_value(q, i);
        }

        let res = will.to_qs().to_string();
        let str_test = "[+-iIjJkKqQxXyYzZ]";
        assert_eq!(res, str_test);
    }

    #[test]
    fn test_correlation() {
        let seq = vec![1,1,-1,1,1,-1,1];

        assert_eq!(periodic_autocorrelation(&seq, 0), 7);
        assert_eq!(periodic_autocorrelation(&seq, 2), -1);


        let seq2 = vec![1,-1,1,1,-1,1,-1];

        assert_eq!(cross_correlation(&seq, &seq2, 1), 5);
        assert_eq!(cross_correlation(&seq, &seq2, 2), -3);

        
        let seq = vec![1,1,-1,1,1,-1];

        assert_eq!(periodic_autocorrelation(&seq, 0), 6);
        assert_eq!(periodic_autocorrelation(&seq, 2), -2);


        let seq2 = vec![1,-1,1,1,-1,1];

        assert_eq!(cross_correlation(&seq, &seq2, 1), 6);
        assert_eq!(cross_correlation(&seq, &seq2, 2), -2);
    }


    
    #[test]
    fn test_periodic_complementary() {
        let mut will = QuadSeq::new(2);
        assert!(!will.is_periodic_complementary());

        let val = [QUADRUPLETS[0], QUADRUPLETS[5]];
        for (i,q) in val.iter().enumerate(){
            will.set_sequence_value(q, i);
        }

        assert!(will.is_periodic_complementary());


        
        let mut will = QuadSeq::new(4);
        assert!(!will.is_periodic_complementary());

        let val = [QUADRUPLETS[0], QUADRUPLETS[0], QUADRUPLETS[1], QUADRUPLETS[0]];
        for (i,q) in val.iter().enumerate(){
            will.set_sequence_value(q, i);
        }

        assert!(will.is_periodic_complementary());
    }


    #[test]
    fn test_symmetric() {
        let mut will = QuadSeq::new(4);
        assert!(will.is_symmetric());

        let val = [QUADRUPLETS[0], QUADRUPLETS[0], QUADRUPLETS[15], QUADRUPLETS[0]];
        for (i,q) in val.iter().enumerate(){
            will.set_sequence_value(q, i);
        }

        assert!(will.is_symmetric());

        will.set_sequence_value(&QUADRUPLETS[3], 3);
        assert!(!will.is_symmetric());
    }


    #[test]
    fn test_amicable() {
        let mut will = QuadSeq::new(4);
        assert!(will.is_amicable());

        let val = [QUADRUPLETS[0], QUADRUPLETS[0], QUADRUPLETS[15], QUADRUPLETS[0]];
        for (i,q) in val.iter().enumerate(){
            will.set_sequence_value(q, i);
        }

        assert!(will.is_amicable());

        will.set_sequence_value(&QUADRUPLETS[3], 3);
        assert!(!will.is_amicable());

    }


    #[test]
    fn test_symmetric_implies_amicable() {
        let size = 5;
        let mut will = QuadSeq::new(size);
    
        aux_recursive(&mut will, size, 1);
    }
    
    fn aux_recursive(will : &mut QuadSeq, size : usize, index : usize) {
    
        if index >= will.search_size(){
            assert!(!will.is_symmetric() || will.is_amicable());
            return;
        }
    
        for value_to_test in QUADRUPLETS.iter(){
            let mut will1 = will.clone();
            will1.set_sequence_value(value_to_test, index);
            aux_recursive(&mut will1, size, index+1);
        }
    }


    #[test]
    fn test_equivalence_cross_correlation_property() {
        let size = 5;
        let mut will = QuadSeq::new(size);
    
        aux_recursive2(&mut will, size, 1);
    }
    
    fn aux_recursive2(will : &mut QuadSeq, size : usize, index : usize) {
    
        if index >= will.search_size(){
            assert!(!(will.verify_cross_correlation() ^ will.is_amicable())); // XNOR operation
            return;
        }
    
        for value_to_test in QUADRUPLETS.iter(){
            let mut will1 = will.clone();
            will1.set_sequence_value(value_to_test, index);
            aux_recursive2(&mut will1, size, index+1);
        }
    }

    #[test]
    fn test_williamson_type() {
        
        let size = 8;
        // rowsum: 0,0,2,6
        let seq_x = vec![-1,-1,-1, 1,1,-1,1, 1];
        let seq_y = vec![-1,-1,-1,-1,1,-1,1,-1];
        let seq_z = vec![-1,-1,-1,-1,1,-1,1,-1];
        let seq_w = vec![-1,-1,-1, 1,1,-1,1, 1];

        let mut will = QuadSeq::new(size);
        will.set_all_values((&seq_x, &seq_y, &seq_z, &seq_w));
        
        assert!(will.to_qs().is_perfect());
    }
}


