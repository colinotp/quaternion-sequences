

#[cfg(test)]
mod tests {
    use crate::sequences::{williamson::{Williamson, QUADRUPLETS}, equivalence::*};

    #[test]
    fn test_coprime() {
        assert!(coprime(5,7));
        assert!(!coprime(16,12));
        assert!(!coprime(42,15));
        assert!(coprime(33,28));
        assert!(coprime(12,11));

        for elm in COPRIMES.iter() {
            println!("{:?}", elm);
        }
    }



    #[test]
    fn equ_reorder() {

        
        let size = 7;
        let seq_x = vec![-1,1,1,-1, 1, 1, 1];
        let seq_y = vec![-1,1,1,-1,-1,-1,-1];
        let seq_z = vec![-1,1,1,-1,-1,-1,-1];
        let seq_w = vec![-1,1,1,-1, 1,-1, 1];

        let mut will = Williamson::new(size);
        will.set_all_values(&(seq_x, seq_y, seq_z, seq_w));

        let equivalent = equivalent_reorder(&will);
        assert_eq!(equivalent.len(), 24);

        for seq in equivalent {
            println!("{}", seq.to_qs().to_string_raw());
            assert!(seq.to_qs().is_perfect());
        }
    }


    #[test]
    fn equ_shift() {

        let size = 10;
        // rowsum: 0,0,2,6
        let seq_x = vec![1,-1,-1,-1,1,1,-1,1,-1,1];
        let seq_y = vec![-1,1,1,1,-1,1,-1,-1,-1,1];
        let seq_z = vec![1,1,1,1,1,-1,-1,1,-1,-1];
        let seq_w = vec![1,1,-1,1,1,1,1,-1,1,1];

        let mut will = Williamson::new(size);
        will.set_all_values(&(seq_x, seq_y, seq_z, seq_w));

        let equivalent = equivalent_shift(&will);
        assert!(equivalent.len() == 10);

        for seq in equivalent {
            println!("{}", seq.to_qs().to_string_raw());
            assert!(seq.to_qs().is_perfect());
        }
    }

    #[test]
    fn equ_negate() {

        let size = 10;
        // rowsum: 0,0,2,6
        let seq_x = vec![1,-1,-1,-1,1,1,-1,1,-1,1];
        let seq_y = vec![-1,1,1,1,-1,1,-1,-1,-1,1];
        let seq_z = vec![1,1,1,1,1,-1,-1,1,-1,-1];
        let seq_w = vec![1,1,-1,1,1,1,1,-1,1,1];

        let mut will = Williamson::new(size);
        will.set_all_values(&(seq_x, seq_y, seq_z, seq_w));

        let equivalent = equivalent_negate(&will);
        assert!(equivalent.len() == 16);

        for seq in equivalent {
            println!("{}", seq.to_qs().to_string_raw());
            assert!(seq.to_qs().is_perfect());
        }
    }

    #[test]
    fn equ_alt_negate() {

        let size = 8;
        let seq_x = vec![-1,-1, 1,-1,-1,1, 1,1];
        let seq_y = vec![-1,-1,-1,-1,-1,1,-1,1];
        let seq_z = vec![-1,-1,-1,-1,-1,1,-1,1];
        let seq_w = vec![-1,-1, 1,-1,-1,1, 1,1];

        let mut will = Williamson::new(size);
        will.set_all_values(&(seq_x, seq_y, seq_z, seq_w));

        let equivalent = equivalent_alternated_negation(&will);
        assert!(equivalent.len() == 2);

        for seq in equivalent {
            assert!(seq.to_qs().is_perfect());
        }
    }

    
    fn test_equivalence_alternate_negation() {
        let size = 6;
        let mut will = Williamson::new(size);
    
        alt_recursive(&mut will, size, 1);
    }
    
    fn alt_recursive(will : &mut Williamson, size : usize, index : usize) {
    
        if index >= will.search_size(){
            if will.to_qs().is_perfect() {
                let equivalent = equivalent_alternated_negation(&will);
                assert!(equivalent.len() == 2);
                println!("{}, {}", equivalent[0].to_qs().to_string_raw(), equivalent[1].to_qs().to_string_raw());
                assert!(equivalent[1].to_qs().is_perfect());
            }
            return;
        }
    
        for value_to_test in QUADRUPLETS.iter(){
            let mut will1 = will.clone();
            will1.set_sequence_value(value_to_test, index);
            alt_recursive(&mut will1, size, index+1);
        }
    }


    fn factorial(n : usize) -> usize {
        if n <= 1 {
            n
        }
        else {
            n*factorial(n-1)
        }
    }


    #[test]
    fn equ_perm() {

        let size = 7;
        let seq_x = vec![-1,1,1,-1, 1, 1, 1];
        let seq_y = vec![-1,1,1,-1,-1,-1,-1];
        let seq_z = vec![-1,1,1,-1,-1,-1,-1];
        let seq_w = vec![-1,1,1,-1, 1,-1, 1];

        let mut will = Williamson::new(size);
        will.set_all_values(&(seq_x, seq_y, seq_z, seq_w));

        let equivalent = equivalent_automorphism(&will);
        // assert_eq!(equivalent.len(), factorial(size));

        for seq in equivalent {
            println!("{}", seq.to_qs().to_string_raw());
            assert!(seq.to_qs().is_perfect());
        }
    
    }
}