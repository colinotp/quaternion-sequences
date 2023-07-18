use super::williamson::Williamson;




pub fn sequence_product(will1 : &Williamson, will2 : &Williamson) -> Williamson{

    let p = will1.size();
    let q = will2.size();

    let mut result = Williamson::new(p*q);

    let (x1,y1,z1,w1) = will1.sequences();
    let (x2,y2,z2,w2) = will2.sequences();

    for i in 0..p {
        for j in 0..q {
            let index = i*q+j;
            let value = (x1[index]*x2[index], y1[index]*y2[index], z1[index]*z2[index], w1[index]*w2[index]);
            result.set_sequence_value(&value, index);
        }
    }

    result
}