
use cgmath::{*, num_traits::PrimInt};

use super::symmetries::*;

pub const N : usize = 16;

pub const Q0 : Quaternion<f32> = Quaternion::new(0.,0.,0.,0.);
pub const Q1 : Quaternion<f32> = Quaternion::new(1.,0.,0.,0.);
pub const QI : Quaternion<f32> = Quaternion::new(0.,1.,0.,0.);
pub const QJ : Quaternion<f32> = Quaternion::new(0.,0.,1.,0.);
pub const QK : Quaternion<f32> = Quaternion::new(0.,0.,0.,1.);
pub const QQ : Quaternion<f32> = Quaternion::new(0.5,0.5,0.5,0.5);

lazy_static! {
    pub static ref QPLUS : Vec<Quaternion<f32>> = vec![Q1, -Q1, QI, -QI, QJ, -QJ, QK, -QK, QQ, QQ*-1., QQ*QI, QQ*-QI, QQ*QJ, QQ*-QJ, QQ*QK, QQ*-QK];
}






#[derive(Clone)]
pub struct QS {
    size: usize,
    values: Vec<Quaternion<f32>>,
    symmetry : Option<Symmetry>
}




impl QS{

    pub fn new(size: usize, symmetry : Option<Symmetry>) -> QS {
        QS {
            values : vec![Q1.clone(); size as usize],
            size,
            symmetry
        }
    }

    pub fn set_values(&mut self, values : Vec<Quaternion<f32>>){
        // replaces the whole sequence
        self.values = values;
    }

    pub fn set_value(&mut self, value : Quaternion<f32>, index: usize){
        // sets a specific value of the sequence, and applies the symmetry
        self.values[index] = value;
        match &self.symmetry {
            Some(Symmetry::I) => {self.values[self.size - 1 - index] = value.clone();}
            Some(Symmetry::II) => {self.values[self.size/2 + index] = (-1).pow(index as u32) as f32 * value.clone();}
            Some(Symmetry::III) => {self.values[self.size/2 + index] = (-1).pow((index/2) as u32) as f32 * value.clone();}
            Some(Symmetry::IV) => {self.values[self.size/2 + index] = -value.clone();}
            None => {}
        }
    }


    pub const fn search_size(&self) -> usize{
        // returns the search size. If there is a symmetry,
        // we only have to search on half of the sequence
        match self.symmetry {
            Some(_) => {(self.size+1)/2}
            None => {self.size}
        }
    }



    pub fn periodic_autocorrelation(&self,t: usize) -> Quaternion<f32> {
        // computes the periodic auto-correlation
        let mut sum_res = Q0.clone();
        for i in 0..self.size{
            sum_res += self.values[i]*(self.values[(i+t)%self.size]).conjugate();
        }
        sum_res
    }

    pub fn odd_periodic_autocorrelation(&self,t : usize) -> Quaternion<f32> {
        // computes the odd periodic auto-correlation
        let mut sum_res = Q0.clone();
        let mut power : f32;
        for i in 0..self.size{
            power = (-1).pow(((i+t)/self.size) as u32) as f32;
            sum_res += self.values[i]*(self.values[(i+t)%self.size]).conjugate() * power;
        }
        sum_res
    }


    pub fn is_perfect(&self) -> bool{
        // tests if the sequence is perfect
        if self.size == 1 {return true;}
        for t in 1..=((self.size+1)/2) { // we only have to check first half, because the second is symmetric to the first 
            if self.periodic_autocorrelation(t) != Q0 {
                return false;
            }
        }
        true
    }

    pub fn is_odd_perfect(&self) -> bool{
        // tests if the sequence is odd perfect
        if self.size == 1 {return true;}
        for t in 1..=((self.size+1)/2) { // we only have to check first half, because the second is symmetric to the first 
            if self.odd_periodic_autocorrelation(t) != Q0 {
                return false;
            }
        }
        true
    }

    pub fn is_symmetric(&self) -> bool {
        // tests if the sequence is symmetric
        let n = self.size;
        for t in 1..((self.size+1)/2) {
            if self.values[t] != self.values[n-t] {
                return false;
            }
        }
        true
    }


    pub fn to_string_raw(&self) -> String {
        let mut res_str: String = "".to_owned();

        for q in &self.values{
            res_str.push_str(&quaternion_to_string(&q));
        }

        res_str.to_string()
    }

}





pub fn quaternion_to_string(quat : &Quaternion<f32>) -> String {
    match *quat {
        q if q==Q1 => "+".to_string(),
        q if q==-Q1 => "-".to_string(),
        q if q==QI  => "i".to_string(),
        q if q==-QI => "I".to_string(),
        q if q==QJ  => "j".to_string(),
        q if q==-QJ => "J".to_string(),
        q if q==QK  => "k".to_string(),
        q if q==-QK => "K".to_string(),
        q if q==QQ  => "q".to_string(),
        q if q==QQ*-1.  => "Q".to_string(),
        q if q==QQ*QI  => "x".to_string(),
        q if q==QQ*-QI => "X".to_string(),
        q if q==QQ*QJ  => "y".to_string(),
        q if q==QQ*-QJ => "Y".to_string(),
        q if q==QQ*QK  => "z".to_string(),
        q if q==QQ*-QK => "Z".to_string(),
        _ => panic!("Invalid entry !")
    }
}


impl ToString for QS{
    fn to_string(&self) -> String {
        let mut res_str: String = "[".to_owned();

        for q in &self.values{
            res_str.push_str(&quaternion_to_string(&q));
        }

        res_str.push_str("]");
        res_str.to_string()
    }

    
}




