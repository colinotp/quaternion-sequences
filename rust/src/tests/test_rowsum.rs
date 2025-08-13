


#[cfg(test)]
mod tests {
    use std::{collections::HashSet, env};

    use crate::{find::find_unique::reduce_to_equivalence, read_lines, sequences::{equivalence::{filter_by_rowsums, generate_equivalent_quad_seqs}, rowsum::*, sequence::QS, symmetries::SequenceType, williamson::QuadSeq}};

    #[test]
    fn test_prop_5() {
        // A few specific test cases

        let w = vec![-1,-1,-1,1,1,1];   // rowsum = 0
        let x = vec![-1,-1,1,1,1,1];    // rowsum = 2
        let y = vec![-1,1,1,1,1,1];     // rowsum = 4
        let z = vec![1,1,1,1,1,1];      // rowsum = 6

        let mut seq : QuadSeq = QuadSeq::new(6);
        seq.set_all_values((&w,&x,&y,&z));
        println!("{}", seq.to_string());

        assert!(has_sorted_rowsums(&seq));


        let w = vec![-1,-1,1,1,1,1];        // rowsum = 2
        let x = vec![-1,-1,1,1,1,1];        // rowsum = 2
        let y = vec![-1,1,1,1,1,1];         // rowsum = 4
        let z = vec![-1,-1,-1,-1,-1,-1];    // rowsum = -6

        let mut seq : QuadSeq = QuadSeq::new(6);
        seq.set_all_values((&w,&x,&y,&z));
        println!("{}", seq.to_string());

        assert!(!has_sorted_rowsums(&seq));


        let w = vec![-1,-1,1,1,1,1];        // rowsum = 2
        let x = vec![-1,-1,1,1,1,1];        // rowsum = 2
        let y = vec![-1,1,1,1,1,1];         // rowsum = 4
        let z = vec![1,1,1,1,1,1];    // rowsum = 6

        let mut seq : QuadSeq = QuadSeq::new(6);
        seq.set_all_values((&w,&x,&y,&z));
        println!("{}", seq.to_string());

        assert!(has_sorted_rowsums(&seq));

        
        // Verify that a list filtered by Proposition 5 filters no sequences that are not filtered by NS
        let mut seqs = vec![];
        let pathname = "results/pairs/qts/find_9/result.seq";
        let seqtype = SequenceType::Hadamard;

        println!("{:?}",env::current_dir());
        println!("{pathname}");
        for line_res in read_lines(&pathname).expect("error reading the file") {
            let line = line_res.expect("Error reading line");
            println!("{}", &line);
            seqs.push(QS::from_str(&line.to_string()));
        }

        let quad_seqs : Vec<QuadSeq> = seqs.into_iter().map(|s| QuadSeq::from_pqs(&s)).collect();

        for quad_seq in &quad_seqs {
            quad_seq.verify(seqtype.clone());
        }

        // List to be filtered
        let all = generate_equivalent_quad_seqs(&quad_seqs, seqtype.clone());

        let reduced = reduce_to_equivalence(&all, seqtype.clone()); // Sequences reduced by NS
        let prop5 = filter_by_rowsums(&all); // Sequences filtered by prop 5
        let prop5_reduced = reduce_to_equivalence(&prop5, seqtype.clone()); // Reducing the list filtered by prop 5

        // Turn both into sets for comparison
        let mut reduced_set = HashSet::new();
        let mut prop5_set = HashSet::new();

        reduced_set.extend(reduced);
        prop5_set.extend(prop5_reduced);

        assert_eq!(reduced_set, prop5_set, "Proposition 5 filtering too many sequences");

    }

    #[test]
    fn test_rowsum_gen() {
        for p in 1..23 {
            println!("\nn={}", p);
            for rs in generate_rowsums(p) {
                println!("{:?}", rs);
                assert_eq!(rs.0*rs.0 + rs.1*rs.1 + rs.2*rs.2 + rs.3*rs.3, 4*p as isize);
            }
        }
    }

    #[test]
    fn test_rowsum(){

        let seq = vec![1;18];
        assert_eq!(rowsum(seq), 18);

        let seq = vec![1,-1,-1,-1,1,-1,1,1,-1,-1,1,1];
        assert_eq!(rowsum(seq), 0);

    }

    fn choose(n : usize, k : usize) -> usize{
        if k == 0 {return 1}
        return (n * choose(n - 1, k - 1)) / k
    }


    #[test]
    fn test_generate_sequence(){
        
        let tests = vec![(4,3), (6,4), (2,1), (4,0), (7,3)];
        
        for (n,k) in tests{
            let mut seq : Vec<i8> = vec![-1;n];
            let res = gen_seq_rec(&mut seq, k, 0);
            assert_eq!(res.len(), choose(n,k));
        }
    }


    #[test]
    fn test_rowsum_generate(){
        let tests = vec![(4,2), (6,4), (7,-3), (4,0)];
        
        for (n,r) in tests{
            let res = generate_sequences_with_rowsum(r, n);
            assert_eq!(res.len(), choose(n,(n as isize + r) as usize/2));
        }

        
        let res = generate_sequences_with_rowsum(4, 3);
        assert_eq!(res.len(), 0);
    }


    #[test]
    fn test_four_squares(){

        let list = sum_of_four_squares(4*19);

        assert_eq!(list.len(), 5);
        assert_eq!(list[0], (0,2,6,6));
        assert_eq!(list[1], (1,1,5,7));
        assert_eq!(list[2], (1,5,5,5));
        assert_eq!(list[3], (2,2,2,8));
        assert_eq!(list[4], (3,3,3,7));

    }

    #[test]
    fn test_four_squares2(){

        let list = sum_of_four_squares(4*17);

        for elm in list {
            println!("{:?}", elm);
        }

    }


    #[test]
    fn test_gen_quads(){

        let quads = generate_other_quadruplets(&(1,1,5,7), SequenceType::QuaternionType);
        assert_eq!(quads.len(), 2);
        assert_eq!(quads[0], (1,1,5,7));
        assert_eq!(quads[1], (-1,1,5,7));
        
        let quads = generate_other_quadruplets(&(1,3,3,7), SequenceType::QuaternionType);
        assert_eq!(quads.len(), 2);
        assert_eq!(quads[0], (1,3,3,7));
        assert_eq!(quads[1], (-1,3,3,7));

        let quads = generate_other_quadruplets(&(3,3,3,8), SequenceType::QuaternionType);
        assert_eq!(quads.len(), 2);
        assert_eq!(quads[0], (3,3,3,8));
        assert_eq!(quads[1], (-3,3,3,8));

        let quads = generate_other_quadruplets(&(3,3,5,5), SequenceType::QuaternionType);
        assert_eq!(quads.len(), 2);
        assert_eq!(quads[0], (3,3,5,5));
        assert_eq!(quads[1], (-3,3,5,5));
    }



    #[test]
    fn test_possible_rowsums(){

        let rowsums = generate_rowsums(17);
        
        assert_eq!(rowsums.len(), 4);
        assert_eq!(rowsums[0], (1,3,3,7));
        assert_eq!(rowsums[1], (-1,3,3,7));
        assert_eq!(rowsums[2], (3,3,5,5));
        assert_eq!(rowsums[3], (-3,3,5,5));
    }



    #[test]
    fn test_aux(){
        
        let res = generate_sequences_with_rowsum(7, 11);

        for elm in res {
            println!("{:?}", elm);
        }
        
    }
}