use cgmath::Quaternion;

use crate::sequences::{equivalence::generate_equivalence_class, symmetries::SequenceType};

use super::sequence::{QS, QPLUS, Q24};


#[derive(Eq, PartialEq, PartialOrd, Ord, Clone, Hash, Debug)]
pub enum SequenceTag { // enum for choosing a specific sequence
    X, Y, Z, W
}

pub fn tag_to_string(tag : &SequenceTag) -> String{
    match *tag {
        SequenceTag::X => {"X".to_string()}
        SequenceTag::Y => {"Y".to_string()}
        SequenceTag::Z => {"Z".to_string()}
        SequenceTag::W => {"W".to_string()}
    }
}

pub static QUADRUPLETS : [(i8,i8,i8,i8); 16] = [ // a table of all quadruplets possible
    (-1,-1,-1,-1),
    (1,1,1,1),
    (1,-1,-1,1), 
    (-1,1,1,-1), 
    (1,1,-1,-1), 
    (-1,-1,1,1), 
    (1,-1,1,-1), 
    (-1,1,-1,1), 
    (1,-1,-1,-1), 
    (-1,1,1,1), 
    (1,1,-1,1), 
    (-1,-1,1,-1), 
    (1,1,1,-1), 
    (-1,-1,-1,1), 
    (1,-1,1,1), 
    (-1,1,-1,-1)
];

pub fn quaternion_to_quad(quat : &Quaternion<f32>) -> (i8, i8, i8, i8) {

    let mut iterator = Q24.iter().enumerate().filter(|(_,q)| *q == quat);
    let index = iterator.next();

    match index {
        None => {panic!("Unrecognized quaternion!")}
        Some((i,_)) => {QUADRUPLETS[i]}
    }
}





#[derive(Clone, PartialEq, Eq, Hash)]
pub struct QuadSeq {
    size : usize,
    a: Vec<i8>,
    b: Vec<i8>,
    c: Vec<i8>,
    d: Vec<i8>,
}




impl QuadSeq{

    pub fn new(size: usize) -> QuadSeq {
        QuadSeq {
            size,
            a : vec![-1; size],
            b : vec![-1; size],
            c : vec![-1; size],
            d : vec![-1; size],
        }
    }

    pub const fn search_size(&self) -> usize{
        // returns the size of the search
        self.size
    }

    pub const fn size(&self) -> usize{
        // returns the size of the search
        self.size
    }

    pub fn sequence(&self, seqtag : SequenceTag) -> Vec<i8>{
        match seqtag {
            SequenceTag::X => {self.a.clone()}
            SequenceTag::Y => {self.b.clone()}
            SequenceTag::Z => {self.c.clone()}
            SequenceTag::W => {self.d.clone()}
        }
    }

    pub fn sequences(&self) -> (Vec<i8>,Vec<i8>,Vec<i8>,Vec<i8>){
        (self.a.clone(), self.b.clone(), self.c.clone(), self.d.clone())
    }

    pub fn values(&self, index : usize) -> (i8, i8, i8, i8) {
        // gets a quadruplet of values at a specific index
        (self.a[index], self.b[index], self.c[index], self.d[index])
    }

    pub fn set_all_values(&mut self, values : (&Vec<i8>,&Vec<i8>,&Vec<i8>,&Vec<i8>)){
        // sets the sequences to specific values
        let (a,b,c,d) = values;
        self.a = a.clone();
        self.b = b.clone();
        self.c = c.clone();
        self.d = d.clone();
    }

    pub fn set_sequence(&mut self, seq : &Vec<i8>, tag : &SequenceTag){
        // sets a value for a specific sequence
        match &tag {
            SequenceTag::X => self.a = seq.clone(),
            SequenceTag::Y => self.b = seq.clone(),
            SequenceTag::Z => self.c = seq.clone(),
            SequenceTag::W => self.d = seq.clone()
        }
    }

    pub fn set_single_value(&mut self, value : i8, tag : &SequenceTag, index: usize){
        // sets a value for a specific sequence and specific index
        match &tag {
            SequenceTag::X => self.a[index] = value,
            SequenceTag::Y => self.b[index] = value,
            SequenceTag::Z => self.c[index] = value,
            SequenceTag::W => self.d[index] = value
        }
    }

    pub fn set_sequence_value(&mut self, value : &(i8, i8, i8, i8), index: usize){
        // sets a value for all 4 sequences on a specific index
        // EXPECTED ORDER IS X Y Z W
        self.a[index] = value.0;
        self.b[index] = value.1;
        self.c[index] = value.2;
        self.d[index] = value.3;
    }


    pub fn from_pqs(pqs : &QS) -> QuadSeq {

        let mut qts = QuadSeq::new(pqs.size());

        let (mut seqx, mut seqy, mut seqz, mut seqw) = (vec![], vec![], vec![], vec![]);

        for elm in pqs.values() {
            let (x,y,z,w) = quaternion_to_quad(&elm);
            seqx.push(x);
            seqy.push(y);
            seqz.push(z);
            seqw.push(w);
        }

        qts.set_all_values((&seqx, &seqy, &seqz, &seqw));

        qts
    }

    // Verifies that self is valid seqtype
    pub fn verify(&self, seqtype : SequenceType) -> bool {
        match seqtype {
            SequenceType::QuaternionType => {self.verify_qts()},
            SequenceType::WilliamsonType => {self.verify_wts()},
            SequenceType::Williamson => {self.verify_ws()},
            _ => {false}
        }
    }

    pub fn verify_qts(&self) -> bool {
        if !(self.is_periodic_complementary() && self.verify_cross_correlation()) {
            return false;
        }
        true
    }

    pub fn verify_wts(&self) -> bool {
        if !(self.verify_qts() && self.is_amicable()) {
            return false;
        }
        true
    }

    pub fn verify_ws(&self) -> bool {
        if !(self.is_periodic_complementary() && self.is_symmetric()) {
            return false;
        }

        true
    }

    pub fn is_periodic_complementary(&self) -> bool{
        // tests if the sequences are periodic complementary
        for offset in 1..=((self.size-1)) {
            if periodic_autocorrelation(&self.a, offset) + periodic_autocorrelation(&self.b, offset) + periodic_autocorrelation(&self.c, offset) + periodic_autocorrelation(&self.d, offset) != 0 {
                return false;
            }
        }
        true
    }

    pub fn is_amicable(&self) -> bool { // This function is a stronger version of the condition in verify_cross_correlation
        for offset in 1..self.size {
            if !(cross_correlation(&self.a, &self.b, offset) == cross_correlation(&self.b, &self.a, offset) &&
               cross_correlation(&self.a, &self.c, offset) == cross_correlation(&self.c, &self.a, offset) &&
               cross_correlation(&self.a, &self.d, offset) == cross_correlation(&self.d, &self.a, offset) &&
               cross_correlation(&self.b, &self.c, offset) == cross_correlation(&self.c, &self.b, offset) &&
               cross_correlation(&self.b, &self.d, offset) == cross_correlation(&self.d, &self.b, offset) &&
               cross_correlation(&self.c, &self.d, offset) == cross_correlation(&self.d, &self.c, offset))
               {
                return false;
            }
        }

        true
    }

    pub fn verify_cross_correlation(&self) -> bool { // This function is a weaker version of the condition in is_amicable
        for offset in 0..self.size {
            if !(cross_correlation(&self.a, &self.b, offset) - cross_correlation(&self.b, &self.a, offset) == cross_correlation(&self.d, &self.c, offset) - cross_correlation(&self.c, &self.d, offset) &&
               cross_correlation(&self.a, &self.c, offset) - cross_correlation(&self.c, &self.a, offset) == cross_correlation(&self.b, &self.d, offset) - cross_correlation(&self.d, &self.b, offset) &&
               cross_correlation(&self.a, &self.d, offset) - cross_correlation(&self.d, &self.a, offset) == cross_correlation(&self.c, &self.b, offset) - cross_correlation(&self.b, &self.c, offset))
               {
                return false;
            }
        }

        true
    }

    pub fn is_symmetric(&self) -> bool {
        // tests if the sequence is symmetric
        let n = self.size;
        for t in 1..=((self.size)/2) { // Trying half the values is sufficient
            if self.values(t) != self.values(n-t) {
                return false;
            }
        }
        true
    }

    // Checks if self is equivalent to another sequence
    pub fn equivalent_to(&self, quad_seq : QuadSeq, seqtype : SequenceType) -> bool {
        generate_equivalence_class(&self, seqtype).contains(&quad_seq)
    }

    pub fn to_qs(&self) -> QS {
        // transforms the sequence to a Quaternion Sequence
        let mut qs = QS::new(self.size, None);

        for i in 0..self.size {
            match (self.a[i], self.b[i], self.c[i], self.d[i]) {
                (1,1,1,1) => qs.set_value(QPLUS[1].clone(), i),
                (1,1,1,-1) => qs.set_value(QPLUS[12].clone(), i),
                (1,1,-1,1) => qs.set_value(QPLUS[10].clone(), i),
                (1,1,-1,-1) => qs.set_value(QPLUS[4].clone(), i),
                (1,-1,1,1) => qs.set_value(QPLUS[14].clone(), i),
                (1,-1,1,-1) => qs.set_value(QPLUS[6].clone(), i),
                (1,-1,-1,1) => qs.set_value(QPLUS[2].clone(), i),
                (1,-1,-1,-1) => qs.set_value(QPLUS[8].clone(), i),
                (-1,1,1,1) => qs.set_value(QPLUS[9].clone(), i),
                (-1,1,1,-1) => qs.set_value(QPLUS[3].clone(), i),
                (-1,1,-1,1) => qs.set_value(QPLUS[7].clone(), i),
                (-1,1,-1,-1) => qs.set_value(QPLUS[15].clone(), i),
                (-1,-1,1,1) => qs.set_value(QPLUS[5].clone(), i),
                (-1,-1,1,-1) => qs.set_value(QPLUS[11].clone(), i),
                (-1,-1,-1,1) => qs.set_value(QPLUS[13].clone(), i),
                (-1,-1,-1,-1) => qs.set_value(QPLUS[0].clone(), i),
                _ => {panic!();}
            }
        }

        qs
    }

}


pub fn periodic_autocorrelation(seq : &Vec<i8>, offset : usize) -> isize {
    // computes the periodic auto correlation fo the sequence
    let mut res = 0;
    let n = seq.len();
    for i in 0..n {
        res += (seq[i]*seq[(i + offset) % n]) as isize;
    }
    
    res
}


pub fn cross_correlation(seq1 : &Vec<i8>, seq2 : &Vec<i8>, offset : usize) -> isize {
    // computes the periodic cross correlation fo the sequences
    assert!(seq1.len() == seq2.len());

    let n = seq1.len();
    let mut res = 0;
    for i in 0..n {
        res += (seq1[i]*seq2[(i + n - offset) % n]) as isize;
    }
    
    res
}



fn element_to_string(elem : i8) -> String {
    match elem {
        1 => "+".to_string(),
        -1 => "-".to_string(),
        _ => {panic!("invalid entry")}
    }
}

impl ToString for QuadSeq{
    fn to_string(&self) -> String {
        let mut res_str: String = "[\n".to_owned();

        let sequences = vec![&self.a, &self.b, &self.c, &self.d];
        for seq in sequences.iter() {
            res_str.push_str("  [");
            for i in 0..self.size{
                res_str.push_str(&element_to_string(seq[i]));
            }
            res_str.push_str("]\n");
        }

        res_str.push_str("]");
        res_str.to_string()
    }

    
}


pub fn seq_less_than(seq1 : &Vec<i8>, seq2 : &Vec<i8>) -> bool {

    let mut index = 0;
    while index < seq1.len() {
        if seq1[index] < seq2[index] {
            return true;
        }
        else if seq1[index] > seq2[index] {
            return false;
        }
        else {
            index += 1;
        }
    }

    false
}
