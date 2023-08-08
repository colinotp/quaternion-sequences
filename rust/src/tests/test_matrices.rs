#[cfg(test)]
mod tests {

    use cgmath::Quaternion;

    use crate::sequences::{sequence::QS, matrices::{QHM, HM, OpMat}, williamson::Williamson, symmetries::SequenceType, matrix_equivalence::reduce_to_hadamard_equivalence, equivalence::generate_equivalence_classes};



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
    
    }

    #[test]
    fn matrix_from_will() {
        let size = 3;
        let seq_x = vec![-1,-1, 1];
        let seq_y = vec![-1,-1,-1];
        let seq_z = vec![-1,-1,-1];
        let seq_w = vec![-1,-1,-1];

        let mut will = Williamson::new(size);
        will.set_all_values((&seq_x, &seq_y, &seq_z, &seq_w));

        let hm = HM::from_williamson(&will, SequenceType::WilliamsonType);

        println!("{}", hm.to_string());

    }


    #[test]
    fn matrix_equ2() {
        // ! WAY TOO LONG !

        let size = 3;
        let seq_x = vec![-1,-1,1];
        let seq_y = vec![-1,-1,-1];
        let seq_z = vec![-1,-1,-1];
        let seq_w = vec![-1,-1,-1];

        let mut will1 = Williamson::new(size);
        will1.set_all_values((&seq_x, &seq_y, &seq_z, &seq_w));


        let size = 3;
        let seq_x = vec![1,1,-1];
        let seq_y = vec![1,1,1];
        let seq_z = vec![1,1,1];
        let seq_w = vec![-1,-1,-1];

        let mut will2 = Williamson::new(size);
        will2.set_all_values((&seq_x, &seq_y, &seq_z, &seq_w));

        let wills = vec![will1, will2];
        let all = generate_equivalence_classes(&wills);

        println!("{}", all.len());

        let liste: Vec<HM> = all.iter().map(|w| HM::from_williamson(w,SequenceType::WilliamsonType)).collect();

        let equ = reduce_to_hadamard_equivalence(&liste);

        println!("{}", equ.len());

        for elm in equ {
            println!("{}", elm.to_string());
        }

    }

}