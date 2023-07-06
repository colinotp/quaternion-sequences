use itertools::{iproduct, Itertools};

use super::{williamson::{Williamson}};



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

pub fn generate_equivalent(seq : &Williamson) -> Williamson {
    // This function generates the representative of the equivalence class that seq belongs to
    
    let mut lexical_minimum = seq.clone();

    minimize_wrt_equivalence(&mut lexical_minimum, Box::new(equivalent_reorder)); // E1: reorder
    minimize_wrt_equivalence(&mut lexical_minimum, Box::new(equivalent_negate)); // E2 : negate
    minimize_wrt_equivalence(&mut lexical_minimum, Box::new(equivalent_shift)); // E3 : shift
    minimize_wrt_equivalence(&mut lexical_minimum, Box::new(equivalent_automorphism)); // E4 : automorphisms
    minimize_wrt_equivalence(&mut lexical_minimum, Box::new(equivalent_alternated_negation)); // E5 : alternate negation
    
    lexical_minimum
}


pub fn minimize_wrt_equivalence(will : &mut Williamson, equivalence : Box<dyn Fn (&Williamson) -> Vec<Williamson>>) {
    
    let mut best = will.clone();

    for current in equivalence(will) {
        if will_less_than(&current, &best) {
            best = current;
        }
    }

    will.set_all_values(&best.sequences());
}



pub fn equivalent_reorder(seq : &Williamson) -> Vec<Williamson> {
    // computes all equivalent sequences by reorder

    let mut res = vec![];

    let (a,b,c,d) = seq.sequences();

    for quad in [a,b,c,d].iter().permutations(4) {
        let mut s = Williamson::new(seq.size());
        let values = (quad[0].clone(), quad[1].clone(), quad[2].clone(), quad[3].clone());
        s.set_all_values(&values);
        res.push(s);
    }

    res
}




pub fn equivalent_shift(seq : &Williamson) -> Vec<Williamson> {
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

    for quads in iproduct!([a.clone(), negated(&a)], [b.clone(), negated(&b)], [c.clone(), negated(&c)], [d.clone(), negated(&d)]) {
        // this macro does the cartesian product of whatever sequences is inside
        // The result is looping over all possible quadruplets or their negated counterpart
        let mut s = Williamson::new(seq.size());
        s.set_all_values(&quads);

        res.push(s);
    }

    res
}


pub fn equivalent_alternated_negation(seq : &Williamson) -> Vec<Williamson> {

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