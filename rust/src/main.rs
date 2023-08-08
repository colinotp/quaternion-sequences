#![allow(dead_code)]
#[macro_use]
extern crate lazy_static;

use std::fs::File;
use std::io::{self, BufRead, Write};
use std::path::Path;

use sequences::matrices::QHM;
use sequences::matrix_equivalence::hadamard_equivalence_from_file;
use time::*;

mod sequences;
mod tests;
mod find;
use crate::find::*;
use crate::find::find_unique::reduce_to_equivalence;
use crate::sequences::{sequence::*, symmetries::*};

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
    for i in 1..15{
        find_williamson_of_size(i);
    }
}

fn find_williamson_of_size(i : usize){

    let now = Instant::now();
    let count = find_williamson::find(i, |will| {will.is_symmetric() && will.is_periodic_complementary()});
    let elapsed_time = now.elapsed().as_seconds_f32();

    eprintln!("For n = {i}, the function took: {elapsed_time} seconds and found {count} sequences");
}

fn find_williamson_type(){
    for i in 1..15{
        find_williamson_of_size(i);
    }
}

fn find_williamson_type_of_size(i : usize){

    let now = Instant::now();
    let count = find_williamson::find(i, |will| {will.is_amicable() && will.is_periodic_complementary()});
    let elapsed_time = now.elapsed().as_seconds_f32();

    eprintln!("For n = {i}, the function took: {elapsed_time} seconds and found {count} sequences");
}


fn find_unique_williamson_type_of_size(i : usize){

    let now = Instant::now();
    let result = find_unique::find(i);
    let elapsed_time = now.elapsed().as_seconds_f32();

    eprintln!("For n = {i}, the function took: {elapsed_time} seconds");

    let s = &("./results/sequences/unique_wts/".to_string() + &i.to_string() + &".seq");
    let path = Path::new(s);
    let mut f = File::create(path).expect("Invalid file ?");
    
    f.write(result.as_bytes()).expect("Error when writing in the file");
}

fn find_write_wts(i : usize){

    let result = find_write::join_pairs(i);

    let s = &("./results/pairs/wts/find_".to_string() + &i.to_string() + &"/result.seq");
    let path = Path::new(s);
    let mut f = File::create(path).expect("Invalid file ?");
    
    let res_string = result.iter().map(|w| w.to_qs().to_string_raw() + &"\n").fold("".to_string(), |s, t| s + &t);

    f.write(res_string.as_bytes()).expect("Error when writing in the file");
}





fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, { // compact code to read a file
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

fn convert_qs_to_matrices() {
    for i in 1.. {
        println!("{}", &("./results/sequences/unique_wts/".to_string() + &i.to_string() + &".seq"));
        if let Ok(lines) = read_lines(&("./results/sequences/unique_wts/".to_string() + &i.to_string() + &".seq")) {
            // Consumes the iterator, returns an (Optional) String
            let s = &("./results/sequences/qhm/".to_string() + &i.to_string() + &".mat");
            let path = Path::new(s);
            let mut f = File::create(path).expect("Invalid file ?");

            let mut result = "".to_string();
            for line in lines {
                if let Ok(pqs) = line {
                    let mut qhm = QHM::from_pqs(QS::from_str(&pqs));
                    qhm.dephase();
                    result += &qhm.to_string();
                    result += &"\n";
                }
            }
            f.write(result.as_bytes()).expect("Error when writing in the file");

            println!("converted sequences of size {i}");
        }
        else {
            break;
        }
    }
}


fn find_matching_algorithm(p : usize) {

    let now = Instant::now();
    let sequences = find_with_rowsum::find_matching(p);
    eprintln!("The function found {} sequences before equivalences in {} seconds", sequences.len(), now.elapsed().as_seconds_f32());
    let result = reduce_to_equivalence(&sequences);
    let count = result.len();
    let elapsed_time = now.elapsed().as_seconds_f32();
    eprintln!("For n = {p}, the function took: {elapsed_time} seconds and found {count} sequences");


    let s = &("./results/sequences/matches/".to_string() + &p.to_string() + &".seq");
    let path = Path::new(s);
    let mut f = File::create(path).expect("Invalid file ?");

    let mut result_string = "".to_string();
    for seq in result {
        assert!(&seq.to_qs().is_perfect());
        result_string += &seq.to_qs().to_string_raw();
        result_string += &"\n";
    }
    f.write(result_string.as_bytes()).expect("Error when writing in the file");
    
}


fn main() {
    let args : Vec<String> = std::env::args().collect();
    let count = args.len();

    if count == 1 {
        //find_pqs(None);
        //find_williamson();
        //find_pqs(Some(Symmetry::I));
        //convert_qs_to_matrices();
        //find_unique_williamson_type_of_size(9);
        //find_q24(8, None);
        //find_write::write_pairs(7);
        hadamard_equivalence_from_file("results/pairs/wts/find_13/result.seq".to_string());
    }
    else if count == 2 {
        match &args[1] {
            s if s == "pqs" => {find_pqs(None)}
            s if s == "ws" => {find_williamson()}
            s if s == "wts" => {find_williamson_type()}
            _ => {}
        }
    }
    else if count == 3 {
        let i = match str::parse::<usize>(&args[2]){
                    Ok(a) => {a},
                    Err(_) => {panic!("argument isn't an integer !")}};

        match &args[1] {
            s if s == "pqs" => {find_pqs_of_type(i, &None)}
            s if s == "ws" => {find_williamson_of_size(i)}
            s if s == "wts" => {find_williamson_type_of_size(i)}
            s if s == "unique" => {find_unique_williamson_type_of_size(i)}
            s if s == "equation" => {find_with_rowsum::find(i, SequenceType::WilliamsonType)}
            s if s == "matching" => {find_matching_algorithm(i)}
            s if s == "pairs" => {find_write::write_pairs(i)}
            s if s == "join" => {find_write_wts(i)}
            s if s == "convert" => {hadamard_equivalence_from_file("results/pairs/wts/find_".to_string() + &i.to_string() + &"/result.seq".to_string())}
            _ => {}
        }
    }



}
