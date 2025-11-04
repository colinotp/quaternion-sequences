use std::{collections::{HashMap, HashSet}, fs::File, io::Write, path::Path};

use itertools::Itertools;
use petgraph::{graph::NodeIndex, Graph, Undirected};

use crate::{read_lines, sequences::{equivalence::ns_canonical, symmetries::SequenceType, williamson::QuadSeq}};

use super::{matrices::HM, sequence::QS};

use rayon::{iter::*};

use graph_canon::{self, CanonLabeling};


pub fn graph_from_hm(mat : &HM) -> Graph<i32,i32,Undirected> {
    let mut g = Graph::new_undirected();
    
    let size = mat.size();

    for i in 0..4*size{
        let node = g.add_node(0);
        if i < 2*size {
            g.add_edge(node, node, 0);
        }
    }

    for row in 0..size {
        for col in 0..size {
            if mat.get(row,col) == 1 {
                g.add_edge(NodeIndex::new(row), NodeIndex::new(2*size+col), 0);
                g.add_edge(NodeIndex::new(size + row), NodeIndex::new(3*size+col), 0);
            }
            else{
                g.add_edge(NodeIndex::new(row), NodeIndex::new(3*size+col), 0);
                g.add_edge(NodeIndex::new(size + row), NodeIndex::new(2*size+col), 0);
            }
        }
    }

    g
}


pub fn are_isomorphic(g1 : Graph<i32,i32,Undirected>, g2 : Graph<i32,i32,Undirected>) -> bool {

    let canon1 = graph_canon::CanonLabeling::new(&g1);
    let canon2 = graph_canon::CanonLabeling::new(&g2);

    canon1 == canon2    
}

fn canon_hm(mat : &HM) -> CanonLabeling {
    graph_canon::CanonLabeling::new(&graph_from_hm(mat))
}


pub fn reduce_to_hadamard_equivalence( mats : &Vec<HM>) -> Vec<&HM> {
    mats.iter().unique_by(|mat|{canon_hm(mat)}).collect_vec()
}

pub fn reduce_to_ns_equivalence(sequences : &Vec<QuadSeq>) -> Vec<QuadSeq> {
    let set : HashSet<QuadSeq> = sequences.iter().map(|seq| ns_canonical(seq)).collect();

    set.into_iter().collect()
}

fn is_ns_canonical(seq : &QuadSeq) -> bool {
    ns_canonical(seq) == *seq
}

pub fn hadamard_equivalence_from_file(pathname : String, seqtype : SequenceType) {

    let mut seqs = vec![];

    // Taking input from the list of filtered sequences means the sequences have already been reduced via QT equivalence operations
    println!("Converting sequences found in {pathname} to Hadamard matrices up to Hadamard equivalence ...");
    for line_res in read_lines(&pathname).expect(&format!("Error reading file. Make sure sequences have already been generated for this length (e.g., {} should exist and not be empty)", pathname)) {
        let line = line_res.expect("Error reading line");

        let qs = QS::from_str(&line.to_string());
        debug_assert!(is_ns_canonical(&QuadSeq::from_pqs(&qs)));
        
        seqs.push(qs);
    }

    let quad_seq_list : Vec<QuadSeq> = seqs.into_iter().map(|s| QuadSeq::from_pqs(&s)).collect();

    for quad_seq in &quad_seq_list {
        quad_seq.verify(seqtype);
    }

    // Reduce via graph isomorphism checking
    println!("Reducing matrices to equivalence via graph isomorphism...");
    let canon_reps : HashMap<CanonLabeling, HM> = quad_seq_list.par_iter().map(|seq| {
        let hmat = HM::from_williamson(seq, SequenceType::QuaternionType);
        (canon_hm(&hmat), hmat)
    }).collect();
    
    let equ = canon_reps.into_iter().map(|(_,v)| v).collect::<Vec<_>>();

    let count = equ.len();

    println!("Number of matrices up to equivalence : {}", /*equ.len()*/ count);

    //println!("From sequences :\n");
    //for mat in &equ {
    //    println!("{}", mat.get_qts().to_string());
    //}

    let input_file = Path::new(&pathname);
    let result_path = input_file.parent().expect("Invalid file").join("result.mat");
    let mut result_file = File::create(result_path).expect("Invalid file ?");

    let mut result_string = "".to_string();
    for mat in equ {
        result_string += &mat.to_string_magma();
        result_string += &"\n";
    }

    result_file.write(result_string.as_bytes()).expect("Error when writing in the file");

}
