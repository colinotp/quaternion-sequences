#[macro_use]
extern crate lazy_static;

use cgmath::{*, num_traits::PrimInt};
use time::*;

const N : usize = 16;

const Q0 : Quaternion<f32> = Quaternion::new(0.,0.,0.,0.);
const Q1 : Quaternion<f32> = Quaternion::new(1.,0.,0.,0.);
const QI : Quaternion<f32> = Quaternion::new(0.,1.,0.,0.);
const QJ : Quaternion<f32> = Quaternion::new(0.,0.,1.,0.);
const QK : Quaternion<f32> = Quaternion::new(0.,0.,0.,1.);
const QQ : Quaternion<f32> = Quaternion::new(0.5,0.5,0.5,0.5);

lazy_static! {
    static ref QPLUS : Vec<Quaternion<f32>> = vec![Q1, -Q1, QI, -QI, QJ, -QJ, QK, -QK, QQ, QQ*-1., QQ*QI, QQ*-QI, QQ*QJ, QQ*-QJ, QQ*QK, QQ*-QK];
}

fn quaternion_to_string(quat : &Quaternion<f32>) -> String {
    match *quat {
        q if q==Q1 => "+".to_string(),
        q if q==-Q1 => "-".to_string(),
        q if q==QI  => "i".to_string(),
        q if q==-QI => "I".to_string(),
        q if q==QJ  => "o".to_string(),
        q if q==-QJ => "O".to_string(),
        q if q==QK  => "u".to_string(),
        q if q==-QK => "U".to_string(),
        q if q==QQ  => "q".to_string(),
        q if q==QQ*-1.  => "Q".to_string(),
        q if q==QQ*QI  => "ì".to_string(),
        q if q==QQ*-QI => "Ì".to_string(),
        q if q==QQ*QJ  => "ò".to_string(),
        q if q==QQ*-QJ => "Ò".to_string(),
        q if q==QQ*QK  => "ù".to_string(),
        q if q==QQ*-QK => "Ù".to_string(),
        _ => panic!("Invalid entry !")
    }
}




#[derive(Clone)]
struct QS {
    size: usize,
    values: Vec<Quaternion<f32>>
}

impl QS{

    fn new(size: usize) -> QS {
        QS {
            values : vec![Q1.clone(); size as usize],
            size
        }
    }

    fn set_values(&mut self, values : Vec<Quaternion<f32>>){
        self.values = values;
    }

    fn set_value(&mut self, value : Quaternion<f32>, index: usize){
        self.values[index] = value;
    }



    fn periodic_autocorrelation(&self,t: usize) -> Quaternion<f32> {

        let mut sum_res = Q0.clone();
        for i in 0..self.size{
            sum_res += self.values[i]*(self.values[(i+t)%self.size]).conjugate();
        }
        sum_res
    }

    fn odd_periodic_autocorrelation(&self,t : usize) -> Quaternion<f32> {

        let mut sum_res = Q0.clone();
        let mut power : f32;
        for i in 0..self.size{
            power = (-1).pow(((i+t)/self.size) as u32) as f32;
            sum_res += self.values[i]*(self.values[(i+t)%self.size]).conjugate() * power;
        }
        sum_res
    }


    fn is_pqs(&self) -> bool{
        for t in 1..self.size {
            if self.periodic_autocorrelation(t) != Q0 {
                return false;
            }
        }
        true
    }

    fn is_opqs(&self) -> bool{
        for t in 1..self.size {
            if self.odd_periodic_autocorrelation(t) != Q0 {
                return false;
            }
        }
        true
    }


    fn sub_opqs(self) -> Option<QS>{

        let mut size;
        let mut pqs;
        for start in 1..(self.size-1)/2{
            size = self.size - 2*start;
            pqs = QS::new(size);
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





fn find(size : usize){
    let mut count = 0;
    let mut pqs = QS::new(size);

    find_recursive(&mut pqs, size, 1, &mut count);

    //println!("total number: {count}")
}

fn find_recursive(pqs : &mut QS, size : usize, index : usize, count: &mut usize){

    if index >= (size+1)/2{
        if pqs.is_opqs(){
            *count+=1;
            //let s = pqs.to_string();
            //println!("{s}");
        }
        return
    }

    for index_value_to_test in 0..N{
        pqs.set_value(QPLUS[index_value_to_test], index);
        pqs.set_value(QPLUS[index_value_to_test], size-1-index);

        find_recursive(pqs, size, index+1, count);
    }
}





fn main() {
    let mut now;
    let mut elapsed_time;

    for i in 1..18{
        now = Instant::now();
        find(i);
        elapsed_time = now.elapsed().as_seconds_f32();

        println!("For n = {i}, the function took: {elapsed_time} seconds");
    }

}
