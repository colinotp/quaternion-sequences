use super::{sequence::QS, mum::MUM};


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
