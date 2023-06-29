use super::sequence::{QS, QPLUS};


#[derive(Eq,PartialEq,PartialOrd,Ord)]
pub enum SequenceTag { // enum for choosing a specific sequence
    A, B, C, D
}

pub fn tag_to_string(tag : &SequenceTag) -> String{
    match *tag {
        SequenceTag::A => {"A".to_string()}
        SequenceTag::B => {"B".to_string()}
        SequenceTag::C => {"C".to_string()}
        SequenceTag::D => {"D".to_string()}
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




#[derive(Clone)]
pub struct Williamson {
    size : usize,
    a: Vec<i8>,
    b: Vec<i8>,
    c: Vec<i8>,
    d: Vec<i8>,
}




impl Williamson{

    pub fn new(size: usize) -> Williamson {
        Williamson {
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

    pub fn values(&self, index : usize) -> (i8, i8, i8, i8) {
        // gets a quadruplet of values at a specific index
        (self.a[index], self.b[index], self.c[index], self.d[index])
    }

    pub fn set_all_values(&mut self, values : (Vec<i8>,Vec<i8>,Vec<i8>,Vec<i8>)){
        // sets the sequences to specific values
        let (a,b,c,d) = values;
        self.a = a;
        self.b = b;
        self.c = c;
        self.d = d;
    }

    pub fn set_sequence(&mut self, seq : Vec<i8>, tag : &SequenceTag){
        // sets a value for a specific sequence
        match &tag {
            SequenceTag::A => self.a = seq,
            SequenceTag::B => self.b = seq,
            SequenceTag::C => self.c = seq,
            SequenceTag::D => self.d = seq
        }
    }

    pub fn set_single_value(&mut self, value : i8, tag : &SequenceTag, index: usize){
        // sets a value for a specific sequence and specific index
        match &tag {
            SequenceTag::A => self.a[index] = value,
            SequenceTag::B => self.b[index] = value,
            SequenceTag::C => self.c[index] = value,
            SequenceTag::D => self.d[index] = value
        }
    }

    pub fn set_sequence_value(&mut self, value : &(i8, i8, i8, i8), index: usize){
        // sets a value for all 4 sequences on a specific index
        self.a[index] = value.0;
        self.b[index] = value.1;
        self.c[index] = value.2;
        self.d[index] = value.3;
    }


    pub fn is_periodic_complementary(&self) -> bool{
        // tests if the sequences are periodic complementary
        for offset in 1..=((self.size+1)/2) { // we only have to check first half, because the second is symmetric to the first 
            if periodic_autocorrelation(&self.a, offset) + periodic_autocorrelation(&self.b, offset) + periodic_autocorrelation(&self.c, offset) + periodic_autocorrelation(&self.d, offset) != 0 {
                return false;
            }
        }
        true
    }

    pub fn is_amicable(&self) -> bool { // This function and the verify_cross_correlation function are equivalent
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

    pub fn verify_cross_correlation(&self) -> bool { // This function and the is_amicable function are equivalent
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


pub fn periodic_autocorrelation(seq : &Vec<i8>, offset : usize) -> i32 {
    // computes the periodic auto correlation fo the sequence
    let mut res = 0;
    let n = seq.len();
    for i in 0..n {
        res += (seq[i]*seq[(i + offset) % n]) as i32;
    }
    
    res
}


pub fn cross_correlation(seq1 : &Vec<i8>, seq2 : &Vec<i8>, offset : usize) -> i32 {
    // computes the periodic cross correlation fo the sequences
    assert!(seq1.len() == seq2.len());

    let n = seq1.len();
    let mut res = 0;
    for i in 0..n {
        res += (seq1[i]*seq2[(i + n - offset) % n]) as i32;
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

impl ToString for Williamson{
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




