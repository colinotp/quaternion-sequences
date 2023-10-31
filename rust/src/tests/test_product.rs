

#[cfg(test)]
mod tests {
    use crate::sequences::{product::{sequence_product, mum_product}, williamson::Williamson, matrices::QHM, mum::{MUM, HMUO}};


    #[test]
    fn test_sequence_product(){
        
        let size = 7;
        let seq_x = vec![-1,1,1,-1, 1, 1, 1];
        let seq_y = vec![-1,1,1,-1,-1,-1,-1];
        let seq_z = vec![-1,1,1,-1,-1,-1,-1];
        let seq_w = vec![-1,1,1,-1, 1,-1, 1];

        let mut will1 = Williamson::new(size);
        will1.set_all_values((&seq_x, &seq_y, &seq_z, &seq_w));
        let qs1 = will1.to_qs();
        
        let size = 10;
        let seq_x = vec![1,-1,-1,-1,1,1,-1,1,-1,1];
        let seq_y = vec![-1,1,1,1,-1,1,-1,-1,-1,1];
        let seq_z = vec![1,1,1,1,1,-1,-1,1,-1,-1];
        let seq_w = vec![1,1,-1,1,1,1,1,-1,1,1];

        let mut will2 = Williamson::new(size);
        will2.set_all_values((&seq_x, &seq_y, &seq_z, &seq_w));
        let qs2 = will2.to_qs();

        let product = sequence_product(&qs1, &qs2);

        println!("{}", product.to_string_raw());
        assert!(product.is_perfect());
    }


    #[test]
    fn test_mum_product() {

        
        let size = 3;
        let seq_x = vec![-1,-1, -1];
        let seq_y = vec![-1,-1, 1];
        let seq_z = vec![-1,-1, 1];
        let seq_w = vec![-1,-1, 1];

        let mut will1 = Williamson::new(size);
        will1.set_all_values((&seq_x, &seq_y, &seq_z, &seq_w));
        let qs1 = will1.to_qs();
        let qhm1 = QHM::from_pqs(qs1).dephased();
        let mum1= MUM::from_qhm(&qhm1);


        let size = 2;
        let seq_x = vec![-1, 1];
        let seq_y = vec![-1,-1];
        let seq_z = vec![-1,-1];
        let seq_w = vec![-1, 1];

        let mut will2 = Williamson::new(size);
        will2.set_all_values((&seq_x, &seq_y, &seq_z, &seq_w));
        let qs2 = will2.to_qs();
        let qhm2 = QHM::from_pqs(qs2).dephased();
        let mum2= MUM::from_qhm(&qhm2);

        let product = mum_product(&mum1, &mum2);

        println!("length: {}, matrix size: {}", product.length(), product.matrix_size());
    }




    #[test]
    fn test_different_products() {

        // ++Q
        let size = 3;
        let seq_x = vec![-1,-1, -1];
        let seq_y = vec![-1,-1, 1];
        let seq_z = vec![-1,-1, 1];
        let seq_w = vec![-1,-1, 1];
        
        let mut will1 = Williamson::new(size);
        will1.set_all_values((&seq_x, &seq_y, &seq_z, &seq_w));
        let qs1 = will1.to_qs();
        
        // +YIQ
        let size = 4;
        let seq_x = vec![-1,-1, 1,-1];
        let seq_y = vec![-1,-1,-1, 1];
        let seq_z = vec![-1, 1,-1, 1];
        let seq_w = vec![-1, 1, 1, 1];
        
        let mut will2 = Williamson::new(size);
        will2.set_all_values((&seq_x, &seq_y, &seq_z, &seq_w));
        let qs2 = will2.to_qs();

        let qs_product = sequence_product(&qs1, &qs2);

        println!("{}", qs1.to_string());
        println!("{}", qs2.to_string());
        println!("{}", qs_product.to_string());

        let qhm1 = QHM::from_pqs(qs1).dephased();
        let mum1 = MUM::from_qhm(&qhm1);
        let hmuo1 = HMUO::from_qhm(&qhm1);
        let qhm2 = QHM::from_pqs(qs2).dephased();
        let mum2 = MUM::from_qhm(&qhm2);
        let hmuo2 = HMUO::from_qhm(&qhm2);


        println!("{}", qhm1.to_string());
        println!("{}", qhm2.to_string());
        println!("{}", hmuo1.to_string());
        println!("{}", hmuo2.to_string());
        //println!("{}", mum1.to_string());
        //println!("{}", mum2.to_string());

        let qhm_product = QHM::from_pqs(qs_product).dephased();
        let hmuo_product = HMUO::from_qhm(&qhm_product);
        
        println!("{}", qhm_product.to_string());
        println!("{}", hmuo_product.to_string());

        let hmuo_tensor = hmuo1.tensor(&hmuo2);

        println!("{}", hmuo_tensor.to_string());

        let from_mum_product = mum_product(&mum1, &mum2);
        let from_qs_product = MUM::from_qhm(&qhm_product);

        println!("MUM product result: length: {}, matrix size: {}", from_mum_product.length(), from_mum_product.matrix_size());
        println!("PQS product result: length: {}, matrix size: {}", from_qs_product.length(), from_qs_product.matrix_size());

        //println!("{}",from_mum_product.to_string());
        //println!("{}",from_qs_product.to_string());


        assert!(false);
    }


}