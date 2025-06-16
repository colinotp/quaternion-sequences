

#[cfg(test)]
mod tests {

    use crate::sequences::{williamson::{QuadSeq, QUADRUPLETS}, equivalence::*};
    use crate::sequences::sequence::*;
    use crate::find::find_unique::reduce_to_equivalence;
    use crate::read_lines;

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
        // TODO : fix
        
        let size = 7;
        let seq_x = vec![-1,1,1,-1, 1, 1, 1];
        let seq_y = vec![-1,1,1,-1,-1,-1,-1];
        let seq_z = vec![-1,1,1,-1,-1,-1,-1];
        let seq_w = vec![-1,1,1,-1, 1,-1, 1];

        let mut will = QuadSeq::new(size);
        will.set_all_values((&seq_x, &seq_y, &seq_z, &seq_w));

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

        let mut will = QuadSeq::new(size);
        will.set_all_values((&seq_x, &seq_y, &seq_z, &seq_w));

        let equivalent = equivalent_uniform_shift(&will);
        assert!(equivalent.len() == 10);

        for seq in equivalent {
            println!("{}", seq.to_qs().to_string_raw());
            assert!(seq.to_qs().is_perfect());
        }
    }

    #[test]
    fn equ_negate() {
        // TODO : fix

        let size = 10;
        // rowsum: 0,0,2,6
        let seq_x = vec![1,-1,-1,-1,1,1,-1,1,-1,1];
        let seq_y = vec![-1,1,1,1,-1,1,-1,-1,-1,1];
        let seq_z = vec![1,1,1,1,1,-1,-1,1,-1,-1];
        let seq_w = vec![1,1,-1,1,1,1,1,-1,1,1];

        let mut will = QuadSeq::new(size);
        will.set_all_values((&seq_x, &seq_y, &seq_z, &seq_w));

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

        let mut will = QuadSeq::new(size);
        will.set_all_values((&seq_x, &seq_y, &seq_z, &seq_w));

        let equivalent = equivalent_alternated_negation(&will);
        assert!(equivalent.len() == 2);

        for seq in equivalent {
            assert!(seq.to_qs().is_perfect());
        }
    }

    #[test]
    fn test_equivalence_class() {
        let size = 5;
        let mut will = QuadSeq::new(size);
    
        alt_recursive(&mut will, size, 1);
    }
    
    fn alt_recursive(will : &mut QuadSeq, size : usize, index : usize) {
    
        if index >= will.search_size(){
            if will.to_qs().is_perfect() {
                let equivalent = generate_equivalence_class(&will);
                for seq in equivalent {
                    assert!(seq.to_qs().is_perfect());
                }
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

        let mut will = QuadSeq::new(size);
        will.set_all_values((&seq_x, &seq_y, &seq_z, &seq_w));

        let equivalent = equivalent_automorphism(&will);
        // assert_eq!(equivalent.len(), factorial(size));

        for seq in equivalent {
            println!("{}", seq.to_qs().to_string_raw());
            assert!(seq.to_qs().is_perfect());
        }
    
    }

    #[test]
    fn test_class() {

        let size = 7;
        let seq_x = vec![-1,1,1,-1, 1, 1, 1];
        let seq_y = vec![-1,1,1,-1,-1,-1,-1];
        let seq_z = vec![-1,1,1,-1,-1,-1,-1];
        let seq_w = vec![-1,1,1,-1, 1,-1, 1];

        let mut will = QuadSeq::new(size);
        will.set_all_values((&seq_x, &seq_y, &seq_z, &seq_w));

        for elm in generate_equivalence_class(&will) {
            println!("{}", elm.to_qs().to_string_raw());
            assert!(elm.to_qs().is_perfect());
        }
    }


    #[test]
    fn test_class_size() {
        let i=16;

        println!("{}", &("./results/pairs/wts/find_".to_string() + &i.to_string() + &"/result.seq"));
        if let Ok(lines) = read_lines(&("./results/pairs/wts/find_".to_string() + &i.to_string() + &"/result.seq")) {
            // Consumes the iterator, returns an (Optional) String

            let mut sequences = vec![];
            for line in lines {
                if let Ok(pqs) = line {
                    let seq = QuadSeq::from_pqs(&QS::from_str(&pqs));
                    sequences.push(seq);
                }
            }

            let reducted_sequences = reduce_to_equivalence(&sequences);

            println!("{}", sequences.len());
            println!("{}", reducted_sequences.len());
        }

        assert!(false);
    }


    #[test]
    fn test_wts_qts() {
        for n in 17..18{
            
            let filepath="results/pairs/wts/find_".to_string() + &n.to_string() + &"/result.seq";

            let mut sequences = vec![];

            for line in read_lines(filepath).expect("Invalid file") {
                sequences.push(QuadSeq::from_pqs(&QS::from_str(&line.expect("error reading line"))));
            }

            let mut total = vec![];

            for seq in sequences {
                let mut v = generate_equivalence_class(&seq).into_iter().collect();
                total.append(&mut v);
            }

            println!("{}",total.len());
            for seq in total {
                assert!(seq.verify_cross_correlation());
                assert!(seq.is_amicable());
            }

        }
    }
}