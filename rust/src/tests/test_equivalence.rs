

#[cfg(test)]
mod tests {

    use std::{collections::{HashMap, HashSet}, env, time::Instant};

    use crate::sequences::{equivalence::*, symmetries::SequenceType, williamson::{QuadSeq, QUADRUPLETS}};
    use crate::sequences::sequence::*;
    use crate::find::find_unique::reduce_to_equivalence;
    use crate::read_lines;

    #[test]
    fn test_dual_half_shift() {
        let mut q = QuadSeq::new(4);

        let a = vec![1,1,-1,-1];
        let b = vec![1,1,-1,-1];
        let c = vec![1,1,-1,-1];
        let d = vec![1,1,-1,-1];

        q.set_all_values((&a,&b,&c,&d));

        let shifted = equivalent_dual_half_shift(&q, SequenceType::QuaternionType, true);

        for seq in shifted {
            println!("{}\n", seq.to_string());
        }
    }

    #[test]
    fn test_qt_canon() {
        let a = vec![-1,-1,-1,-1,-1,1,-1,1];
        let b = vec![-1,-1,-1,-1,-1,1,-1,1];
        let c = vec![-1,-1,1,-1,-1,1,1,1];
        let d = vec![-1,1,1,1,-1,-1,1,-1];

        let mut qs = QuadSeq::new(8);
        qs.set_all_values((&a, &b, &c, &d));

        let sym = generate_symmetry_group(8, SequenceType::QuaternionType, &vec![equivalent_alternated_negation, equivalent_automorphism, equivalent_uniform_shift]);

        println!("Canonical form:\n{}", qt_canonical(&qs, &sym).to_string());
    }

    #[test]
    fn reduce_equiv() {
        let mut seqs = vec![];
        let pathname = "results/pairs/qts/find_9/result.seq";
        let seqtype = SequenceType::QuaternionType;

        println!("{:?}",env::current_dir());
        println!("{pathname}");
        for line_res in read_lines(&pathname).expect("error reading the file") {
            let line = line_res.expect("Error reading line");
            println!("{}", &line);
            seqs.push(QS::from_str(&line.to_string()));
        }

        let quad_seqs : Vec<QuadSeq> = seqs.into_iter().map(|s| QuadSeq::from_pqs(&s)).collect();

        for quad_seq in &quad_seqs {
            quad_seq.verify(seqtype);
        }

        let mut set1in = HashSet::new();
        let mut set2in = HashSet::new();
        set1in.extend(quad_seqs.clone().into_iter());
        set2in.extend(quad_seqs.clone().into_iter());
        let vec1in : Vec<QuadSeq> = set1in.into_iter().collect();
        let vec2in : Vec<QuadSeq> = set2in.into_iter().collect();

        assert_ne!(vec1in, vec2in, "sets arranged the same way");

        let vec1out = reduce_to_equivalence(&vec1in, seqtype, &vec![equivalent_negate_swap]);
        let vec2out = reduce_to_equivalence(&vec2in, seqtype, &vec![equivalent_negate_swap]);
        let mut set1out = HashSet::new();
        let mut set2out = HashSet::new();
        set1out.extend(vec1out.into_iter());
        set2out.extend(vec2out.into_iter());

        assert_eq!(set1out, set2out);
    }

    #[test]
    fn equ_runtimes() {
        let w = vec![-1,-1,-1,-1,-1,-1,-1,-1,1,1,-1,1,1];
        let x = vec![-1,-1,1,1,1,1,-1,-1,-1,1,-1,1,-1];
        let y = vec![-1,1,-1,1,1,-1,1,-1,1,1,1,1,1];
        let z = vec![1,-1,-1,1,1,-1,-1,1,1,-1,1,-1,1];

        let mut qs = QuadSeq::new(3);
        qs.set_all_values((&w, &x, &y, &z));

        assert!(qs.verify(SequenceType::QuaternionType));

        let iterations = 10000;

        let mut times = HashMap::new();

        let avg_negate_swap = (0..iterations).map(|_| {
            let time = Instant::now();
            equivalent_negate_swap(&qs, SequenceType::QuaternionType, true);
            time.elapsed().as_nanos()
        }).sum::<u128>() / iterations;
        times.insert("Negate swap", avg_negate_swap);

        let avg_reorder = (0..iterations).map(|_| {
            let time = Instant::now();
            equivalent_reorder(&qs, SequenceType::QuaternionType, true);
            time.elapsed().as_nanos()
        }).sum::<u128>() / iterations;
        times.insert("Reorder", avg_reorder);

        let avg_double_reorder = (0..iterations).map(|_| {
            let time = Instant::now();
            equivalent_double_reorder(&qs, SequenceType::QuaternionType, true);
            time.elapsed().as_nanos()
        }).sum::<u128>() / iterations;
        times.insert("Double reorder", avg_double_reorder);

        let avg_uni_half_shift = (0..iterations).map(|_| {
            let time = Instant::now();
            equivalent_uniform_half_shift(&qs, SequenceType::QuaternionType, true);
            time.elapsed().as_nanos()
        }).sum::<u128>() / iterations;
        times.insert("Uniform half shift", avg_uni_half_shift);

        let avg_uni_shift = (0..iterations).map(|_| {
            let time = Instant::now();
            equivalent_uniform_shift(&qs, SequenceType::QuaternionType, true);
            time.elapsed().as_nanos()
        }).sum::<u128>() / iterations;
        times.insert("Uniform shift", avg_uni_shift);

        let avg_reverse = (0..iterations).map(|_| {
            let time = Instant::now();
            equivalent_reverse(&qs, SequenceType::QuaternionType, true);
            time.elapsed().as_nanos()
        }).sum::<u128>() / iterations;
        times.insert("Reverse", avg_reverse);

        let avg_negate = (0..iterations).map(|_| {
            let time = Instant::now();
            equivalent_negate(&qs, SequenceType::QuaternionType, true);
            time.elapsed().as_nanos()
        }).sum::<u128>() / iterations;
        times.insert("Negate", avg_negate);

        let avg_double_negate = (0..iterations).map(|_| {
            let time = Instant::now();
            equivalent_double_negate(&qs, SequenceType::QuaternionType, true);
            time.elapsed().as_nanos()
        }).sum::<u128>() / iterations;
        times.insert("Double negate", avg_double_negate);

        let avg_alt_negate = (0..iterations).map(|_| {
            let time = Instant::now();
            equivalent_alternated_negation(&qs, SequenceType::QuaternionType, true);
            time.elapsed().as_nanos()
        }).sum::<u128>() / iterations;
        times.insert("Alternated negation", avg_alt_negate);

        let avg_even_alt_negate = (0..iterations).map(|_| {
            let time = Instant::now();
            equivalent_even_alternated_negation(&qs, SequenceType::QuaternionType, true);
            time.elapsed().as_nanos()
        }).sum::<u128>() / iterations;
        times.insert("Even alternated negation", avg_even_alt_negate);

        let avg_automorphism = (0..iterations).map(|_| {
            let time = Instant::now();
            equivalent_automorphism(&qs, SequenceType::QuaternionType, true);
            time.elapsed().as_nanos()
        }).sum::<u128>() / iterations;
        times.insert("Automorphism", avg_automorphism);

        let avg_disjoint_swap = (0..iterations).map(|_| {
            let time = Instant::now();
            equivalent_disjoint_swaps(&qs, SequenceType::QuaternionType, true);
            time.elapsed().as_nanos()
        }).sum::<u128>() / iterations;
        times.insert("Disjoint swaps", avg_disjoint_swap);

        let mut sorted : Vec<_> = times.iter().collect();
        sorted.sort_by(|a, b| b.1.cmp(a.1));

        for data in sorted {
            println!("{:<25} avg runtime: {}ns", data.0, data.1);
        }
    }

    #[test]
    fn test_equiv() {
        let a = vec![-1,-1,-1];
        let b = vec![-1,-1,1];
        let c = vec![-1,-1,1];
        let d = vec![-1,-1,1];

        let mut qs = QuadSeq::new(3);
        qs.set_all_values((&a, &b, &c, &d));

        let symmetries = generate_symmetry_group(3, SequenceType::QuaternionType, &SequenceType::QuaternionType.equivalences());
        let equiv = generate_equivalence_class_fast(&qs, &symmetries);
        let equiv_old = generate_equivalence_class(&qs, SequenceType::QuaternionType, &SequenceType::QuaternionType.equivalences(), false);

        for elm in equiv.clone() {
            assert!(equiv_old.contains(&elm), "New fn returns seq not found in unoptimized build: {}", elm.to_string());
        }
        for elm in equiv_old.clone() {
            assert!(equiv.contains(&elm), "Old fn returns seq not found in optimized build: {}", elm.to_string());
        }
    }

    #[test]
    fn test_disjoint_swap() {
        let a = vec![-1,-1,-1];
        let b = vec![-1,-1,1];
        let c = vec![-1,1,-1];
        let d = vec![-1,1,1];

        let mut qs = QuadSeq::new(3);
        qs.set_all_values((&a, &b, &c, &d));

        let h_equ = equivalent_disjoint_swaps(&qs, SequenceType::QuaternionType, false);
        println!("Orig seq:\n{}\nEquivalence class:", qs.to_string());
        for seq in h_equ.clone() {
            println!("{}", seq.to_string());
        }

        assert_eq!(h_equ.len(), 4);
    }

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

        let equivalent = equivalent_double_reorder(&will, SequenceType::QuaternionType, false);
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

        let equivalent = equivalent_uniform_shift(&will, SequenceType::QuaternionType, false);
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

        let equivalent = equivalent_double_negate(&will, SequenceType::QuaternionType, false);
        assert!(equivalent.len() == 16);

        for seq in equivalent {
            println!("{}", seq.to_qs().to_string_raw());
            assert!(seq.to_qs().is_perfect());
        }
    }

    #[test]
    fn equ_negate_swap() {
        let w = vec![-1,-1,-1];
        let x = vec![-1,-1,1];
        let y = vec![-1,-1,1];
        let z = vec![1,1,-1];

        let size = x.len();

        let mut qts = QuadSeq::new(size);
        qts.set_all_values((&w,&x,&y,&z));

        let class = generate_equivalence_class(&qts, SequenceType::QuaternionType, &vec![equivalent_negate_swap], false);
        println!("orig seq: {}\nequivalence class:", qts.to_string());

        for seq in class.iter() {
            println!("{}", seq.to_string());
            assert!(seq.verify(SequenceType::QuaternionType));

            let neg_swap = equivalent_negate_swap(&seq, SequenceType::QuaternionType, false);
            for alt_seq in neg_swap {
                assert!(class.contains(&alt_seq));
            }
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

        let equivalent = equivalent_even_alternated_negation(&will, SequenceType::QuaternionType, false);
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
                let equivalent = generate_equivalence_class(&will, SequenceType::QuaternionType, &SequenceType::QuaternionType.equivalences(), false);
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

        let equivalent = equivalent_automorphism(&will, SequenceType::QuaternionType, false);
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

        for elm in generate_equivalence_class(&will, SequenceType::QuaternionType, &SequenceType::QuaternionType.equivalences(), false) {
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

            let reducted_sequences = reduce_to_equivalence(&sequences, SequenceType::QuaternionType, &SequenceType::QuaternionType.equivalences());

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
                let mut v = generate_equivalence_class(&seq, SequenceType::QuaternionType, &SequenceType::QuaternionType.equivalences(), false).into_iter().collect();
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