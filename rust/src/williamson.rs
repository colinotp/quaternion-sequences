use crate::sequence::{QS, QPLUS};


#[derive(Eq,PartialEq,PartialOrd,Ord)]
pub enum SequenceTag {
    A, B, C, D
}

pub static SEQUENCE_TAGS : [SequenceTag; 4] = [SequenceTag::A, SequenceTag::B, SequenceTag::C, SequenceTag::D];

pub static QUADRUPLETS : [(i8,i8,i8,i8); 16] = [
    (1,1,1,1),
    (1,1,1,-1), 
    (1,1,-1,1), 
    (1,1,-1,-1), 
    (1,-1,1,1), 
    (1,-1,1,-1), 
    (1,-1,-1,1), 
    (1,-1,-1,-1), 
    (-1,1,1,1), 
    (-1,1,1,-1), 
    (-1,1,-1,1), 
    (-1,1,-1,-1), 
    (-1,-1,1,1), 
    (-1,-1,1,-1), 
    (-1,-1,-1,1), 
    (-1,-1,-1,-1)
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
            a : vec![1; size],
            b : vec![1; size],
            c : vec![1; size],
            d : vec![1; size],
        }
    }

    pub const fn search_size(&self) -> usize{
        self.size
    }

    pub fn set_all_values(&mut self, values : (Vec<i8>,Vec<i8>,Vec<i8>,Vec<i8>)){
        let (a,b,c,d) = values;
        self.a = a;
        self.b = b;
        self.c = c;
        self.d = d;
    }

    pub fn set_single_value(&mut self, value : i8, tag : &SequenceTag, index: usize){
        match &tag {
            SequenceTag::A => self.a[index] = value,
            SequenceTag::B => self.b[index] = value,
            SequenceTag::C => self.c[index] = value,
            SequenceTag::D => self.d[index] = value
        }
    }

    pub fn set_sequence_value(&mut self, value : &(i8, i8, i8, i8), index: usize){
        self.a[index] = value.0;
        self.b[index] = value.1;
        self.c[index] = value.2;
        self.d[index] = value.3;
    }



    pub fn periodic_autocorrelation(&self,t: usize) -> i32 {

        let mut sum_res = 0;
        for i in 0..self.size{
            sum_res += (self.a[i]*self.a[(i+t)%self.size]) as i32;
            sum_res += (self.b[i]*self.b[(i+t)%self.size]) as i32;
            sum_res += (self.c[i]*self.c[(i+t)%self.size]) as i32;
            sum_res += (self.d[i]*self.d[(i+t)%self.size]) as i32;
        }
        sum_res
    }

    pub fn odd_periodic_autocorrelation(&self,t: usize) -> i32 {

        let mut sum_res = 0;
        let mut power : i32;
        for i in 0..self.size{
            power = (-1 as i32).pow(((i+t)/self.size) as u32);
            sum_res += power*(self.a[i]*self.a[(i+t)%self.size]) as i32;
            sum_res += power*(self.b[i]*self.b[(i+t)%self.size]) as i32;
            sum_res += power*(self.c[i]*self.c[(i+t)%self.size]) as i32;
            sum_res += power*(self.d[i]*self.d[(i+t)%self.size]) as i32;
        }
        sum_res
    }

    pub fn is_perfect_complementary(&self) -> bool{
        for t in 1..((self.size+1)/2) { // we only have to check first half, because the second is symmetric to the first 
            if self.periodic_autocorrelation(t) != 0 {
                return false;
            }
        }
        true
    }

    pub fn is_odd_perfect_complementary(&self) -> bool{
        for t in 1..((self.size+1)/2) { // we only have to check first half, because the second is symmetric to the first 
            if self.odd_periodic_autocorrelation(t) != 0 {
                return false;
            }
        }
        true
    }



    pub fn to_qs(&self) -> QS {
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



pub fn element_to_string(elem : i8) -> String {
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




