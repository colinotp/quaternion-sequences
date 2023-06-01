#[macro_use]
extern crate lazy_static;

use time::*;
use williamson::Williamson;

pub mod sequence;
pub mod find_naive;
pub mod find_optim;
pub mod symmetries;
pub mod williamson;
mod find_williamson;
use crate::sequence::*;
use crate::symmetries::*;

fn print_group() {
    println!("{}", quaternion_to_string(&QQ));
    println!("{}", quaternion_to_string(&(QQ*-1.)));
    println!("{}", quaternion_to_string(&(QQ*QI)));
    println!("{}", quaternion_to_string(&(QQ*QJ)));
    println!("{}", quaternion_to_string(&(QQ*QK)));
    println!("{}", quaternion_to_string(&(QQ*QI*QI)));
    println!("{}", quaternion_to_string(&(QQ*QI*QJ)));
    println!("{}", quaternion_to_string(&(QQ*QI*QK)));
}

fn find_pqs(){
    let mut now;
    let mut elapsed_time;
    let mut count;
    let symmetry = None;


    for i in 1..18{

        match symmetry {
            None | Some(Symmetry::I) => {}
            _ => {
                if i % 2 == 1 {continue}
            },
        }

        now = Instant::now();
        count = find_naive::find(i, symmetry.clone());
        elapsed_time = now.elapsed().as_seconds_f32();
    
        println!("For n = {i}, the function took: {elapsed_time} seconds and found {count} sequences");
    }
}

fn find_williamson(){
    let mut now;
    let mut elapsed_time;
    let mut count;


    for i in 1..18{

        now = Instant::now();
        count = find_williamson::find(i);
        elapsed_time = now.elapsed().as_seconds_f32();
    
        println!("For n = {i}, the function took: {elapsed_time} seconds and found {count} sequences");
    }
}

fn print_williamson(){
    let mut will = Williamson::new(6);
    will.set_single_value(-1, &williamson::SequenceTag::A, 0);
    will.set_single_value(-1, &williamson::SequenceTag::B, 1);
    will.set_single_value(-1, &williamson::SequenceTag::D, 3);
    println!("{}", will.to_string());
    println!("{}", will.to_qs().to_string());
}

fn main() {
    find_williamson();
}
