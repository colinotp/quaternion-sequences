use std::usize::MAX;


use crate::sequences::{rowsum::{generate_rowsums, Quad, generate_sequences_with_rowsum}, fourier::{iter_over_filtered_dft}};

use super::{williamson::SequenceTag, symmetries::SequenceType};



enum OpType {
    PLUS, MINUS
}


pub fn generate_equations(seq1 : Vec<i8>, tag1 : SequenceTag, seq2 : Vec<i8>, tag2 : SequenceTag, seqtype : SequenceType) {

    match seqtype {
        SequenceType::WilliamsonType => {

        }
        SequenceType::ItoType => {

        }
        _ => {/* TODO */}
    }
}


fn equations_williamson_type(seq1 : Vec<i8>, tag1 : SequenceTag, seq2 : Vec<i8>, tag2 : SequenceTag) {

    match (tag1, tag2) {
        (SequenceTag::A, SequenceTag::B) => {
        }
        _ => {/* TODO */}
    }


}

fn equations_crosscorrelation(seq1 : Vec<i8>, op_type1 : OpType, seq2 : Vec<i8>, op_type2 : OpType) {

    let op1 = match op_type1 {
        OpType::PLUS => { |x : i8, y : i8| x + y}
        OpType::MINUS => { |x : i8, y : i8| x - y}
    };
    let op2 = match op_type2 {
        OpType::PLUS => { |x : i8, y : i8| x + y}
        OpType::MINUS => { |x : i8, y : i8| x - y}
    };

    let n = seq1.len();
    let mut values = vec![0 ; 2*n];

    for i in 1..=n {

        for k in 0..n {
            values[k] = op1(seq1[k+1-i], seq1[k+i-1]); 
        }
        
        for k in 0..n {
            values[k+n] = - op2(seq2[k+1-i], seq2[k+i-1]); 
        }

        generate_equation_from(&values);
    }
}


fn generate_equation_from(values : &Vec<i8>) {

    let x = &" x".to_string();
    let space = &" ".to_string();

    let mut result = "".to_string();
    let mut variable_counter = 1;

    for &elm in values.iter() {
        if elm > 0 {
            result += &(elm.to_string() + x + &variable_counter.to_string() + space);
            variable_counter += 1;
        }
    }

    result += ";";
    println!("{}", result);
}