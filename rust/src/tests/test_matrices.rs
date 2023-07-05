#[cfg(test)]
mod tests {
    use cgmath::Quaternion;

    use crate::sequences::{sequence::QS, matrices::QHM};



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
}