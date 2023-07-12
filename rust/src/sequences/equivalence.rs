use std::collections::HashSet;


use itertools::iproduct;

use super::williamson::{Williamson, SequenceTag};



// * Computes a table of coprimes at compile time

pub fn coprime(mut a : usize, mut b : usize) -> bool {
    // tests if two numbers are coprime or not
    if b < a {
        (a,b) = (b,a);
    }
    while b != 0 {
        (a,b) = (b, a % b)
    }
    a == 1
}

const N : usize = 101;

lazy_static!(
    // calculates at compile-time the table of coprimes up to N
    pub static ref COPRIMES : Vec<Vec<usize>> = {
        let mut list = vec![];

        for i in 0..N {

            let mut coprimes = vec![];

            for j in 1..=i {
                if coprime(i,j) {
                    coprimes.push(j);
                }
            }
            list.push(coprimes);
        }

        list
    };
);




// * defining the lexical order on williamson sequences

pub fn seq_less_than(seq1 : &Vec<i8>, seq2 : &Vec<i8>) -> bool {

    let mut index = 0;
    while index < seq1.len() {
        if seq1[index] < seq2[index] {
            return true;
        }
        else if seq1[index] > seq2[index] {
            return false;
        }
        else {
            index += 1;
        }
    }

    false
}


enum Comp {
    LT, EQUAL, GT
}

fn will_less_than_aux(seq1 : &Vec<i8>, seq2 : &Vec<i8>) -> Comp {

    let mut index = 0;
    while index < seq1.len() {
        if seq1[index] < seq2[index] {
            return Comp::LT;
        }
        else if seq1[index] > seq2[index] {
            return Comp::GT;
        }
        else {
            index += 1;
        }
    }

    Comp::EQUAL
}

pub fn will_less_than(will1 : &Williamson, will2 : &Williamson) -> bool {

    let (a1,b1,c1,d1) = will1.sequences();
    let (a2,b2,c2,d2) = will2.sequences();
    match will_less_than_aux(&a1, &a2) {
        Comp::LT => {true}
        Comp::GT => {false}
        Comp::EQUAL => {
            match will_less_than_aux(&b1, &b2) {
                Comp::LT => {true}
                Comp::GT => {false}
                Comp::EQUAL => {
                    match will_less_than_aux(&c1, &c2) {
                        Comp::LT => {true}
                        Comp::GT => {false}
                        Comp::EQUAL => {
                            seq_less_than(&d1, &d2)
                        }
    }}}}}
}



// * Functions to treat the equivalences

pub fn generate_canonical_representative(seq : &Williamson) -> Williamson{
    let set = generate_equivalence_class(seq);
    let mut mini = seq.clone();
    for elm in set {
        if will_less_than(&elm, &mini) {
            mini = elm;
        }
    }
    mini
}




pub fn generate_equivalence_class(seq : &Williamson) -> HashSet<Williamson> {
    // This function generates the representative of the equivalence class that seq belongs to
    
    let mut class = HashSet::new();
    class.insert(seq.clone());

    loop {
        let mut new = HashSet::new();

        for seq in &class {
            for equivalence in [equivalent_reorder, equivalent_negate, equivalent_uniform_shift, equivalent_reorder, equivalent_alternated_negation, equivalent_automorphism] {
                for equ in equivalence(&seq){
                    if !class.contains(&equ) {
                        new.insert(equ);
                    }
                }
            }

        }


        if new.len() == 0 {break;}
        else {
            for seq in new {
                class.insert(seq);
            }
        }
    }

    class
}




fn swap(will : &mut Williamson, seqtag1 : SequenceTag, seqtag2 : SequenceTag) {
    
    let (a,b,c,d) = will.sequences();

    let (seq1, seq2) = match (&seqtag1, &seqtag2) {
        (SequenceTag::A, SequenceTag::B) => {(a, b)},
        (SequenceTag::A, SequenceTag::C) => {(a, c)},
        (SequenceTag::A, SequenceTag::D) => {(a, d)},
        (SequenceTag::B, SequenceTag::C) => {(b, c)},
        (SequenceTag::B, SequenceTag::D) => {(b, d)},
        (SequenceTag::C, SequenceTag::D) => {(c, d)},
        _ => {panic!("Incorrect tags entered !")}
    };

    will.set_sequence(seq2, &seqtag1);
    will.set_sequence(seq1, &seqtag2);
}

pub fn equivalent_reorder(seq : &Williamson) -> Vec<Williamson> {
    // computes all equivalent sequences by reorder

    let mut res = vec![seq.clone()];

    let couples = [(SequenceTag::A, SequenceTag::B), (SequenceTag::A, SequenceTag::C), (SequenceTag::A, SequenceTag::D), (SequenceTag::B, SequenceTag::C), (SequenceTag::B, SequenceTag::D), (SequenceTag::C, SequenceTag::D)];

    for (couple1, couple2) in iproduct!(couples.clone(), couples) {
        let mut new_seq = seq.clone();
        
        let (seq11, seq12) = couple1;
        let (seq21, seq22) = couple2;
        if !(seq11 == seq21 && seq12 == seq22){
            swap(&mut new_seq, seq11, seq12);
            swap(&mut new_seq, seq21, seq22);
        }

        res.push(new_seq);
    
    }

    res
}




pub fn equivalent_uniform_shift(seq : &Williamson) -> Vec<Williamson> {
    // computes all equivalent sequences by shift

    let mut res = vec![seq.clone()];

    for offset in 1..seq.size() {
        let mut s = Williamson::new(seq.size());
        for index in 0..seq.size() {
            s.set_sequence_value(&seq.values((index + offset) % seq.size()), index)
        }
        res.push(s);
    }

    res
}



pub fn negated(seq : &Vec<i8>) -> Vec<i8> {
    let mut s = vec![];
    for i in 0..seq.len() {
        s.push(-seq[i]);
    }
    s
}

pub fn alt_negated(seq : &Vec<i8>, frequency : usize) -> Vec<i8> {
    let mut s = vec![];
    let mut count = 0;
    for i in 0..seq.len() {
        if count % frequency == frequency - 1 {
            s.push(-seq[i]);
        }
        else {
            s.push(seq[i]);
        }
        count += 1;
    }
    s
}



pub fn equivalent_negate(seq : &Williamson) -> Vec<Williamson> {
    // computes all equivalent sequences by negation

    let (a,b,c,d) = seq.sequences();

    let mut res = vec![];

    for tag_couple in [(SequenceTag::A, SequenceTag::B), (SequenceTag::A, SequenceTag::C), (SequenceTag::A, SequenceTag::D), (SequenceTag::B, SequenceTag::C), (SequenceTag::B, SequenceTag::D), (SequenceTag::C, SequenceTag::D)] {
        // this loops through all the couples of a,b,c,d (ordered couples)
        let quad = match tag_couple {
            (SequenceTag::A, SequenceTag::B) => {(negated(&a.clone()), negated(&b.clone()), c.clone(), d.clone())},
            (SequenceTag::A, SequenceTag::C) => {(negated(&a.clone()), b.clone(), negated(&c.clone()), d.clone())},
            (SequenceTag::A, SequenceTag::D) => {(negated(&a.clone()), b.clone(), c.clone(), negated(&d.clone()))},
            (SequenceTag::B, SequenceTag::C) => {(a.clone(), negated(&b.clone()), negated(&c.clone()), d.clone())},
            (SequenceTag::B, SequenceTag::D) => {(a.clone(), negated(&b.clone()), c.clone(), negated(&d.clone()))},
            (SequenceTag::C, SequenceTag::D) => {(a.clone(), b.clone(), negated(&c.clone()), negated(&d.clone()))},
            _ => {panic!("Incorrect tags entered !")}
        };
        let mut s = Williamson::new(seq.size());
        s.set_all_values(&quad);

        res.push(s);
    }

    res
}


pub fn equivalent_alternated_negation(seq : &Williamson) -> Vec<Williamson> {

    if seq.size() % 2 == 1 {
        return vec![seq.clone()];
    }

    let frequency = 2;

    let (a,b,c,d) = seq.sequences();

    let mut res = vec![];

    let quads = (alt_negated(&a, frequency), alt_negated(&b, frequency), alt_negated(&c, frequency), alt_negated(&d, frequency));

    let mut s = Williamson::new(seq.size());
    s.set_all_values(&quads);

    res.push(seq.clone());
    res.push(s);

    res
}




fn permute(seq : &Vec<i8>, coprime : usize) -> Vec<i8> {
    // permutes all elements of seq using
    // the automorphism of the cyclic group defined by a coprime

    let n = seq.len();

    let mut result = vec![];

    for index in 0..n {
        result.push(seq[index * coprime % n]);
    }

    result
}





pub fn equivalent_automorphism(seq : &Williamson) -> Vec<Williamson> {
    // computes all equivalent sequences by permutation

    let mut result = vec![];
    let size = seq.size();

    let mut identity = vec![];
    for i in 0..size {
        identity.push(i);
    }

    for k in COPRIMES[size].iter() {

        let mut will = Williamson::new(size);

        let (a,b,c,d) = seq.sequences();

        let quad = (permute(&a, *k), permute(&b, *k), permute(&c, *k), permute(&d, *k));

        will.set_all_values(&quad);

        result.push(will);

    }


    result
}