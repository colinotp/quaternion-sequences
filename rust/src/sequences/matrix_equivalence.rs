use std::{path::Path, fs::File, io::Write, env, collections::{HashMap}};

use itertools::Itertools;
use petgraph::{graph::NodeIndex, visit::{EdgeIndexable, NodeIndexable}, Graph, Undirected};

use crate::{sequences::{equivalence::generate_equivalence_classes, williamson::Williamson, symmetries::SequenceType}, read_lines};

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



pub fn hadamard_equivalence_from_file(pathname : String) {

    let mut seqs = vec![];

    println!("{:?}",env::current_dir());
    println!("{pathname}");
    for line_res in read_lines(&pathname).expect("error reading the file") {
        let line = line_res.expect("Error reading line");
        println!("{}", &line);
        seqs.push(QS::from_str(&line.to_string()));
    }

    let wills = seqs.iter().map(|s| Williamson::from_pqs(s)).collect();
    let all = generate_equivalence_classes(&wills);

    println!("generated all sequences");

    for elm in &all {
        assert!(elm.to_qs().is_perfect());
    }

    let liste: Vec<HM> = all.iter().map(|w| HM::from_williamson(w,SequenceType::WilliamsonType)).collect();

    println!("generated all matrices");

    let equ = liste.into_par_iter().map(|mat| (canon_hm(&mat), mat)).collect::<HashMap<_,_>>().into_iter().map(|(_,v)| v).collect::<Vec<_>>();

    let count = equ.len();

    println!("number of matrices up to equivalence : {}", /*equ.len()*/ count);

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