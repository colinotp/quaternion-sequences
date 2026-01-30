#[cfg(test)]
mod tests {

    use cgmath::Quaternion;

    use crate::sequences::{sequence::QS, matrices::{QHM, HM, OpMat}, williamson::QuadSeq, symmetries::SequenceType, matrix_equivalence::reduce_to_hadamard_equivalence, equivalence::generate_equivalent_quad_seqs};

    #[test]
    fn qhmtest() {
        let qs = QS::from_str(&"+JJ+x".to_string());
        let qhm = QHM::from_pqs(qs);

        assert!(qhm.verify());
    }

    #[test]
    fn test_qts_from_hm() {
        let mut qts = QuadSeq::new(3);
        let w = vec![-1,-1,-1];
        let x = vec![-1,-1,1];
        let y = vec![-1,-1,1];
        let z = vec![-1,-1,1];
        qts.set_all_values((&w, &x, &y, &z));

        let hm = HM::from_williamson(&qts, SequenceType::QuaternionType);

        let qts_new = hm.get_qts();

        assert_eq!(qts, qts_new);
    }

    #[test]
    fn matrix_from_pqs() {

        let mut pqs : QS = QS::new(5, None);
        let values = vec![Quaternion::new(1.,0.,0.,0.),
                          Quaternion::new(0.,1.,0.,0.),
                          Quaternion::new(0.,1.,0.,0.),
                          Quaternion::new(1.,0.,0.,0.),
                          Quaternion::new(-0.5,-0.5,-0.5,-0.5)];

        pqs.set_values(values);

        let mut qhm = QHM::from_pqs(pqs);
        println!("{}", qhm.to_string());

        qhm.dephase();
        println!("{}", qhm.to_string());

    }


    #[test]
    fn matrix_from_seq() {
        let seq = vec![1,-1,-1,1,-1];

        let hm = HM::from_sequence(&seq);

        println!("{}", hm.to_string());

        assert!(!hm.verify());
    }

    #[test]
    fn matrix_op() {
        let size = 5;
        let seq = vec![1,-1,-1,1,-1];

        let hm = HM::from_sequence(&seq);

        println!("id :\n{}", hm.to_string());
    
        let mut hm2 = HM::new(size);
        hm2.copy_block_to(&hm, 0, 0, &OpMat::MINUS);

        println!("minus: \n{}", hm2.to_string());
    
        assert!(!hm.verify());
    }

    #[test]
    fn matrix_from_will() {
        let size = 3;
        let seq_x = vec![-1,-1, 1];
        let seq_y = vec![-1,-1, 1];
        let seq_z = vec![-1,-1, 1];
        let seq_w = vec![-1,-1,-1];

        let mut will = QuadSeq::new(size);
        will.set_all_values((&seq_x, &seq_y, &seq_z, &seq_w));
        assert!(will.verify(SequenceType::QuaternionType));

        let hm = HM::from_williamson(&will, SequenceType::QuaternionType);

        println!("{}", hm.to_string());

        assert!(hm.verify());
    }


    #[test]
    fn matrix_equ2() {
        let size = 3;
        let seq_x = vec![-1,-1, 1];
        let seq_y = vec![-1,-1, 1];
        let seq_z = vec![-1,-1, 1];
        let seq_w = vec![-1,-1,-1];

        let mut will = QuadSeq::new(size);
        will.set_all_values((&seq_x, &seq_y, &seq_z, &seq_w));
        assert!(will.verify(SequenceType::QuaternionType));

        let wills = vec![will];
        let all = generate_equivalent_quad_seqs(&wills, SequenceType::QuaternionType);

        println!("{}", all.len());

        let liste: Vec<HM> = all.iter().map(|w| HM::from_williamson(w,SequenceType::QuaternionType)).collect();

        let equ = reduce_to_hadamard_equivalence(&liste);

        println!("{}", equ.len());

        for elm in equ {
            println!("{}", elm.to_string());
            assert!(elm.verify());
        }
    }

}
