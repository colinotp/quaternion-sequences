
#[cfg(test)]
mod tests {

    use crate::sequences::{mum::*, sequence::*, matrices::QHM, williamson::Williamson};


    #[test]
    fn test_instance_operator() {

        println!("{}",OP0.to_string());
        println!("{}",OP1.to_string());
        println!("{}",OPX.to_string());
        println!("{}",OPY.to_string());
        println!("{}",OPZ.to_string());

        println!("{}",quaternion_to_operator(&Q24[4]).to_string());

    }

    #[test]
    fn test_operation_operator() {

        println!("{}", OPX.tensor(&OPY).to_string());
        println!("{}", OPY.tensor(&OPX).to_string());
        println!("{}", quaternion_to_operator(&QQ).to_string());
        println!("{}", quaternion_to_operator(&QQ).conjugate_transpose().to_string());

    }

    #[test]
    fn test_mum() {        
        let size = 2;
        let seq_x = vec![-1, 1];
        let seq_y = vec![-1,-1];
        let seq_z = vec![-1,-1];
        let seq_w = vec![-1, 1];

        let mut will = Williamson::new(size);
        will.set_all_values((&seq_x, &seq_y, &seq_z, &seq_w));

        let pqs = will.to_qs();
        let qhm = QHM::from_pqs(pqs).dephased();
        let mum = MUM::from_qhm(&qhm);

        println!("{}", mum.to_string());
    }

}