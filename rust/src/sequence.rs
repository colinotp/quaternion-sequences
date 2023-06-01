
use cgmath::{*, num_traits::PrimInt};

use crate::symmetries::*;

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
        self.values = values;
    }

    pub fn set_value(&mut self, value : Quaternion<f32>, index: usize){
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
        match self.symmetry {
            Some(_) => {(self.size+1)/2}
            None => {self.size}
        }
    }



    pub fn periodic_autocorrelation(&self,t: usize) -> Quaternion<f32> {

        let mut sum_res = Q0.clone();
        for i in 0..self.size{
            sum_res += self.values[i]*(self.values[(i+t)%self.size]).conjugate();
        }
        sum_res
    }

    pub fn odd_periodic_autocorrelation(&self,t : usize) -> Quaternion<f32> {

        let mut sum_res = Q0.clone();
        let mut power : f32;
        for i in 0..self.size{
            power = (-1).pow(((i+t)/self.size) as u32) as f32;
            sum_res += self.values[i]*(self.values[(i+t)%self.size]).conjugate() * power;
        }
        sum_res
    }


    pub fn is_pqs(&self) -> bool{
        for t in 1..((self.size+1)/2) { // we only have to check first half, because the second is symmetric to the first 
            if self.periodic_autocorrelation(t) != Q0 {
                return false;
            }
        }
        true
    }

    pub fn is_opqs(&self) -> bool{
        for t in 1..((self.size+1)/2) { // we only have to check first half, because the second is symmetric to the first 
            if self.odd_periodic_autocorrelation(t) != Q0 {
                return false;
            }
        }
        true
    }


    pub fn sub_opqs(self) -> Option<QS>{

        let mut size;
        let mut pqs;
        for start in 1..(self.size-1)/2{
            size = self.size - 2*start;
            pqs = QS::new(size, self.symmetry.clone());
            pqs.set_values(self.values[start..self.size-start].to_vec());
            if pqs.is_opqs(){
                return Some(pqs);
            }
        }
        None
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




