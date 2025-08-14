use std::{collections::HashMap, env, fs::File, io::Write, path::Path, time::Instant};

use itertools::Itertools;
use petgraph::{graph::NodeIndex, Graph, Undirected};

use crate::{find::find_unique::reduce_to_equivalence, read_lines, sequences::{equivalence::{equivalent_negate_swap, filter_by_rowsums, generate_equivalent_quad_seqs}, symmetries::SequenceType, williamson::QuadSeq}};

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



pub fn hadamard_equivalence_from_file(pathname : String, seqtype : SequenceType) {

    let mut seqs = vec![];

    println!("{:?}",env::current_dir());
    println!("{pathname}");
    for line_res in read_lines(&pathname).expect("error reading the file") {
        let line = line_res.expect("Error reading line");
        println!("{}", &line);
        seqs.push(QS::from_str(&line.to_string()));
    }

    let quad_seq_list : Vec<QuadSeq> = seqs.into_iter().map(|s| QuadSeq::from_pqs(&s)).collect();

    for quad_seq in &quad_seq_list {
        quad_seq.verify(seqtype);
    }

    let all = generate_equivalent_quad_seqs(&quad_seq_list, seqtype);
    let all_size = all.len();

    for elm in &all {
        assert!(elm.to_qs().is_perfect());
    }

    println!("Generated all {} {} including equivalent sequences via {} equivalence operations\n", all_size, seqtype.to_string(), seqtype.to_string());

    println!("Filtering sequences by rowsums...");
    let time = Instant::now();
    let filtered_rowsum = filter_by_rowsums(&all);
    let elapsed = time.elapsed().as_secs();
    println!("Filtering sequences by rowsums took {} seconds. Reduced to {} sequences.\n", elapsed, filtered_rowsum.len());

    // Filter via equivalence operations
    println!("Filtering sequences via Hadamard equivalence operations...");
    let time = Instant::now();
    let reduced = reduce_to_equivalence(&all, seqtype, &vec![equivalent_negate_swap]);
    let elapsed = time.elapsed().as_secs();
    println!("Filtering sequences via equivalence operations took {} seconds. Reduced to {} sequences.\n", elapsed, reduced.len());

    // Fully reduce via graph isomorphism checking
    println!("Reducing matrices to equivalence via graph isomorphism...");
    let canon_reps : HashMap<CanonLabeling, HM> = reduced.par_iter().map(|seq| {
        let hmat = HM::from_williamson(seq, SequenceType::QuaternionType);
        (canon_hm(&hmat), hmat)
    }).collect();
    
    let equ = canon_reps.into_iter().map(|(_,v)| v).collect::<Vec<_>>();

    let count = equ.len();

    println!("Number of matrices up to equivalence : {}", /*equ.len()*/ count);

    let result_name = pathname.split(".").next().expect("No first element").to_string() + &".mat";

    println!("{result_name}");
    let path = Path::new(&result_name);
    let mut result_file = File::create(path).expect("Invalid file ?");

    let mut result_string = "".to_string();
    for mat in equ {
        result_string += &mat.to_string_magma();
        result_string += &"\n";
    }

    result_file.write(result_string.as_bytes()).expect("Error when writing in the file");

}