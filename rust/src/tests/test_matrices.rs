#[cfg(test)]
mod tests {
    use cgmath::Quaternion;

    use crate::sequences::{sequence::QS, matrices::{QHM, HM}, williamson::Williamson, symmetries::SequenceType, matrix_equivalence::generate_equivalence_class};



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
    fn matrix_from_will() {
        let size = 3;
        let seq_x = vec![-1,-1, 1];
        let seq_y = vec![-1,-1,-1];
        let seq_z = vec![-1,-1,-1];
        let seq_w = vec![-1,-1,-1];

        let mut will = Williamson::new(size);
        will.set_all_values((&seq_x, &seq_y, &seq_z, &seq_w));

        let hm = HM::from_williamson(will, SequenceType::WilliamsonType);

        println!("{}", hm.to_string());

    }


    #[test]
    fn matrix_equivalence() {
        // ! WAY TOO LONG !

        let size = 1;
        let seq_x = vec![-1];
        let seq_y = vec![-1];
        let seq_z = vec![-1];
        let seq_w = vec![-1];

        let mut will = Williamson::new(size);
        will.set_all_values((&seq_x, &seq_y, &seq_z, &seq_w));

        let hm = HM::from_williamson(will, SequenceType::WilliamsonType);

        let equ = generate_equivalence_class(&hm);
        for mat in &equ {
            println!("{}", mat.to_string());
        }

        println!("found {} equivalent matrices", equ.len());

    }

}