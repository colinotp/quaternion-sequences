#![allow(dead_code)]
#[macro_use]
extern crate lazy_static;

use time::*;

mod sequences;
mod tests;
mod find;
use crate::find::*;
use crate::sequences::{sequence::*, williamson::*, symmetries::*};

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

fn find_pqs(symmetry : Option<Symmetry>){
    for i in 1..18{
        find_pqs_of_type(i, &symmetry);
    }
}

fn find_pqs_of_type(i : usize, symmetry : &Option<Symmetry>){

    match symmetry {
        None | Some(Symmetry::I) => {}
        _ => {
            if i % 2 == 1 {return}
        },
    }

    let now = Instant::now();
    let count = find_optim::find(i, symmetry.clone());
    let elapsed_time = now.elapsed().as_seconds_f32();

    eprintln!("For n = {i}, the function took: {elapsed_time} seconds and found {count} sequences");
}



fn find_williamson(){
    for i in 1..7{
        find_williamson_of_size(i);
    }
}

fn find_williamson_of_size(i : usize){

    let now = Instant::now();
    let count = find_williamson::find(i);
    let elapsed_time = now.elapsed().as_seconds_f32();

    eprintln!("For n = {i}, the function took: {elapsed_time} seconds and found {count} sequences");
    
}

fn print_williamson(){
    let mut will = Williamson::new(6);
    will.set_single_value(-1, &SequenceTag::A, 0);
    will.set_single_value(-1, &SequenceTag::B, 1);
    will.set_single_value(-1, &SequenceTag::D, 3);
    println!("{}", will.to_string());
    println!("{}", will.to_qs().to_string());
}

fn main() {
    find_pqs(None);
    //find_williamson();
}
