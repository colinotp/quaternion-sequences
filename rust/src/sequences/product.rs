use super::{williamson::Williamson, sequence::QS, mum::MUM};




pub fn williamson_product(will1 : &Williamson, will2 : &Williamson) -> Williamson{
    // ! Doesn't work !

    let p = will1.size();
    let q = will2.size();

    let mut result = Williamson::new(p*q);

    let (x1,y1,z1,w1) = will1.sequences();
    let (x2,y2,z2,w2) = will2.sequences();

    for i in 0..p {
        for j in 0..q {
            let index = i*q+j;
            let value = (x1[index % p]*x2[index % q], y1[index % p]*y2[index % q], z1[index % p]*z2[index % q], w1[index % p]*w2[index % q]);
            result.set_sequence_value(&value, index);
        }
    }

    result
}

pub fn sequence_product(qs1 : &QS, qs2 : &QS) -> QS{

    let p = qs1.size();
    let q = qs2.size();

    let mut result = QS::new(p*q, None);

    let v1 = qs1.values();
    let v2 = qs2.values();

    for i in 0..p {
        for j in 0..q {
            let index = i*q+j;
            let value = v1[index % p]*v2[index % q];
            result.set_value(value, index);
        }
    }

    result
}

pub fn mum_product(mum1 : &MUM, mum2 : &MUM) -> MUM{

    let p = mum1.length();
    let q = mum2.length();

    let mut values = vec![];

    let v1 = mum1.values();
    let v2 = mum2.values();

    for i in 0..p {
        for j in 0..q {
            let index = i*q+j;
            let value = v1[index % p].tensor(&v2[index % q]);
            values.push(value);
        }
    }

    MUM::new(mum1.matrix_size() * mum2.matrix_size(), &values)
}
