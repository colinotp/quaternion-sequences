#![allow(dead_code)]
#[macro_use]
extern crate lazy_static;

use std::{fs::File, env};
use std::io::{self, BufRead, Write};
use std::path::Path;

use sequences::matrices::QHM;
use time::*;

mod sequences;
mod tests;
mod find;
use crate::find::*;
use crate::find::find_unique::reduce_to_equivalence;
use crate::sequences::{williamson::*, sequence::*, symmetries::*};
use sequences::matrix_equivalence::hadamard_equivalence_from_file;

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
        find_williamson_type_of_size(i);
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

fn find_write_quad_seq(i : usize, seqtype : SequenceType){

    let result = find_write::join_pairs(i, seqtype.clone());

    if matches!(seqtype, SequenceType::QuaternionType) {
        // Check to see if also valid WTS
        for qts in &result {
            assert!(qts.verify_qts(), "Sequence failed auto/cross correlation conditions: {}", qts.to_string());
            if !qts.is_amicable() {
                print!("Seq is valid QTS, but not amicable (not WTS): {}", qts.to_string())
            }
        }
    }
    
    let folder = seqtype.to_string();
    
    let s = &("./results/pairs/".to_string() + &folder + &"/find_".to_string() + &i.to_string() + &"/result.seq");
    let qs = &("./results/pairs/".to_string() + &folder + &"/find_".to_string() + &i.to_string() + &"/result.qseq");
    
    let path_seq = Path::new(s);
    let path_qseq = Path::new(qs);

    let mut f_seq = File::create(path_seq).expect("Invalid file ?");
    let mut f_qseq = File::create(path_qseq).expect("Invalid file ?");
    
    let seq_res_string = result.iter().map(|w| w.to_qs().to_string_raw() + &"\n").fold("".to_string(), |s, t| s + &t);
    let qseq_res_string = result.iter().map(|w| w.to_string() + &"\n").fold("".to_string(), |s, t| s + &t);

    f_seq.write(seq_res_string.as_bytes()).expect("Error when writing in the file");
    f_qseq.write(qseq_res_string.as_bytes()).expect("Error when writing in the file");
}


fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, { // compact code to read a file
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

fn convert_qs_to_matrices() {
    for i in 1.. {
        println!("{}", &("./results/pairs/qts/find_".to_string() + &i.to_string() + &"/result.seq"));
        if let Ok(lines) = read_lines(&("./results/pairs/qts/find_".to_string() + &i.to_string() + &"/result.seq")) {
            // Consumes the iterator, returns an (Optional) String
            let s = &("./results/pairs/qts/find_".to_string() + &i.to_string() + &"/result.qhm");
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

// Verify that all qts of a given length are also WTS (amicable)
fn verify_qts_eq_wts(p : usize) {
    let mut seqs = vec![];

    let pathname = "results/pairs/qts/find_".to_string() + &p.to_string() + &"/result.seq".to_string();

    println!("{:?}",env::current_dir());
    println!("{pathname}");
    for line_res in read_lines(&pathname).expect("error reading the file") {
        let line = line_res.expect("Error reading line");
        println!("{}", &line);
        seqs.push(QS::from_str(&line.to_string()));
    }

    let qts_list : Vec<QuadSeq> = seqs.iter().map(|s| QuadSeq::from_pqs(s)).collect();

    for qts in &qts_list {
        assert!(qts.verify_qts(), "Sequence fails auto/cross correlation condition: {}", qts.to_string());
        assert!(qts.is_amicable(), "Valid QTS fails amicabililty condition: {}", qts.to_string());
    }

    println!("Length {} checked, all QTS = WTS", p);
}

fn str_to_rowsum_pairing(n : &String) -> Option<RowsumPairing> {
    match n.as_str() {
        "XY" => Some(RowsumPairing::XY),
        "XZ" => Some(RowsumPairing::XZ),
        "XW" => Some(RowsumPairing::XW),
        _ => None
    }
}

fn str_to_seqtype(n : &str) -> SequenceType {
    match n {
        "qts" => SequenceType::QuaternionType,
        "wts" => SequenceType::WilliamsonType,
        "ws" => SequenceType::Williamson,
        _ => {panic!("Invalid sequence type passed")}
    }
}

fn str_to_usize(source : &str) -> usize {
    match str::parse::<usize>(source) {
        Ok(a) => {a},
        Err(_) => {panic!("argument isn't an integer !")}
    }
}

fn str_to_isize(source : &str) -> isize {
    match str::parse::<isize>(source) {
        Ok(a) => {a},
        Err(_) => {panic!("argument isn't an integer !")}
    }
}

fn main() {
    let args : Vec<String> = std::env::args().collect();

    match args.len() {
        // For running individual functions
        1 => {
            //find_pqs(None);
            //find_williamson();
            //find_pqs(Some(Symmetry::I));
            convert_qs_to_matrices();
            //find_unique_williamson_type_of_size(9);
            //find_q24(8, None);
            //find_write::write_pairs(7);
        },
        // General operations on sequences (mostly deprecated)
        2 => match args[1].as_str() {
            "pqs" => {find_pqs(None);},
            "ws" => {find_williamson();},
            "wts" => {find_williamson_type();},
            _ => {panic!("Invalid argument");}
        },
        // Operations on sequences of a given length p
        3 => {
            let p = str_to_usize(&args[2]);
            match args[1].as_str() {
                "pqs" => {find_pqs_of_type(p, &None);}
                "ws" => {find_williamson_of_size(p);}
                "wts" => {find_williamson_type_of_size(p);}
                "unique" => {find_unique_williamson_type_of_size(p);}
                "equation" => {find_with_rowsum::find(p, SequenceType::QuaternionType);}
                "matching" => {find_matching_algorithm(p);}
                "amicable" => {verify_qts_eq_wts(p);}
                _ => {}
            };
        },
        // Matching step/rowsum generation/Hadamard equivalence reduction
        4 => {
            let seqtype = str_to_seqtype(&args[2]);
            let p = str_to_usize(&args[3]);
            match args[1].as_str() {
                "join" => {find_write_quad_seq(p, seqtype);},
                "rowsums" => {find_write::write_rowsums(p, seqtype);},
                "convert" => {hadamard_equivalence_from_file("results/pairs/".to_string() + &seqtype.to_string() + &"/find_".to_string() + &p.to_string() + &"/result.seq".to_string());}
                _ => {}
            }
        }
        // Pair generation
        5 => {
            let seqtype = str_to_seqtype(&args[2]);
            let p = str_to_usize(&args[3]);
            let pairing = str_to_rowsum_pairing(&args[4]);
            match args[1].as_str() {
                "pairs" => {find_write::write_pairs(p, seqtype, pairing);}
                _ => {}
            }
        },
        // Pair generation for a single set of rowsums (or file generation)
        9 => {
            let folder = str_to_seqtype(&args[2]).to_string();  // verifies seqtype input is correct
            let p = str_to_usize(&args[3]);     // length
            let a = str_to_isize(&args[4]);     // rowsum 1
            let b = str_to_isize(&args[5]);     // rowsum 2
            let c = str_to_isize(&args[6]);     // rowsum 3
            let d = str_to_isize(&args[7]);     // rowsum 4

            let pairing = str_to_rowsum_pairing(&args[8]);      // Rowsum pairing

            match args[1].as_str() {
                "pairs_rowsum" => {find_write::write_pairs_rowsum(&folder, (a,b,c,d), p, pairing)}
                "create" => {find_write::create_rowsum_dirs(folder, p, (a,b,c,d), pairing);}
                _ => {}
            }
        },
        // Pair generation for a single pair from a single set of rowsums
        10 => {
            let folder = str_to_seqtype(&args[2]).to_string();  // verifies seqtype input is correct
            let p = str_to_usize(&args[3]);     // length
            let a = str_to_isize(&args[4]);     // rowsum 1
            let b = str_to_isize(&args[5]);     // rowsum 2
            let c = str_to_isize(&args[6]);     // rowsum 3
            let d = str_to_isize(&args[7]);     // rowsum 4

            let pairing = str_to_rowsum_pairing(&args[8]);      // Rowsum pairing
            
            // First or second pair (1 or 2)?
            let pair = match str::parse::<u8>(&args[9]) {
                Ok(a) => {a},
                Err(_) => {panic!("argument isn't an integer !")}
            };

            match args[1].as_str() {
                "pair_single" => {find_write::write_pair_single_rowsum(folder, (a,b,c,d), p, pairing, pair);},
                _ => {}
            }
        },
        _ => {}
    };
}