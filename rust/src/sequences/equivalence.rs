use std::collections::HashSet;


use itertools::iproduct;

use crate::sequences::{symmetries::SequenceType, williamson::tag_to_string};

use super::williamson::{QuadSeq, SequenceTag};



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

        // Fixes indexing
        list.push(vec![]);

        for i in 1..N {

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

pub fn will_less_than(will1 : &QuadSeq, will2 : &QuadSeq) -> bool {

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

pub fn generate_canonical_representative(seq : &QuadSeq, seqtype : SequenceType) -> QuadSeq{
    let set = generate_equivalence_class(seq, seqtype);
    let mut mini = seq.clone();
    for elm in set {
        if will_less_than(&elm, &mini) {
            mini = elm;
        }
    }
    mini
}




pub fn generate_equivalence_class(seq : &QuadSeq, seqtype : SequenceType) -> HashSet<QuadSeq> {
    // This function generates the equivalence class that seq belongs to
    
    let mut class = HashSet::new();
    class.insert(seq.clone());

    loop {
        let mut new = HashSet::new();

        for seq in &class {
            for equivalence in seqtype.equivalences() {
                for equ in equivalence(&seq, seqtype.clone()).into_iter() {
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


pub fn generate_equivalent_quad_seqs(quad_seq_list : &Vec<QuadSeq>, seqtype : SequenceType) -> Vec<QuadSeq> {

    let mut result = HashSet::new();

    for quad_seq in quad_seq_list {
        if result.contains(quad_seq) {
            continue;
        }

        let class = generate_equivalence_class(quad_seq, seqtype.clone());
        for elm in class {
            result.insert(elm);
        }
    }

    result.into_iter().collect()
}



fn swap(will : &mut QuadSeq, seqtag1 : SequenceTag, seqtag2 : SequenceTag) {
    
    let (a,b,c,d) = will.sequences();

    let (seq1, seq2) = match (&seqtag1, &seqtag2) {
        (SequenceTag::W, SequenceTag::X) => {(a, b)},
        (SequenceTag::W, SequenceTag::Y) => {(a, c)},
        (SequenceTag::W, SequenceTag::Z) => {(a, d)},
        (SequenceTag::X, SequenceTag::Y) => {(b, c)},
        (SequenceTag::X, SequenceTag::Z) => {(b, d)},
        (SequenceTag::Y, SequenceTag::Z) => {(c, d)},
        _ => {panic!("Incorrect tags entered !")}
    };

    will.set_sequence(&seq2, &seqtag1);
    will.set_sequence(&seq1, &seqtag2);
}

// Applies a single swap and a single negation
pub fn equivalent_negate_swap(seq : &QuadSeq, seqtype : SequenceType) -> HashSet<QuadSeq> {
    let mut res : HashSet<QuadSeq> = HashSet::new();
    res.insert(seq.clone());
    let (a,b,c,d) = seq.sequences();
    let (nega_a, nega_b, nega_c, nega_d) = (negated(&a), negated(&b), negated(&c), negated(&d));

    let couples = [(SequenceTag::W, SequenceTag::X), (SequenceTag::W, SequenceTag::Y), (SequenceTag::W, SequenceTag::Z), (SequenceTag::X, SequenceTag::Y), (SequenceTag::X, SequenceTag::Z), (SequenceTag::Y, SequenceTag::Z)];

    for couple in couples {
        for tag in [SequenceTag::W, SequenceTag::X, SequenceTag::Y, SequenceTag::Z] {
            let mut new_seq = seq.clone();
            
            match tag {
                SequenceTag::W => new_seq.set_sequence(&nega_a, &tag),
                SequenceTag::X => new_seq.set_sequence(&nega_b, &tag),
                SequenceTag::Y => new_seq.set_sequence(&nega_c, &tag),
                SequenceTag::Z => new_seq.set_sequence(&nega_d, &tag)
            }

            swap(&mut new_seq, couple.0.clone(), couple.1.clone());
            debug_assert!(new_seq.verify(seqtype.clone()), "equivalent_negate_swap produced bad seq. Init: {}\nRes: {}\nSwap: {}\n Negation: {}", seq.to_string(), new_seq.to_string(), tag_to_string(&couple.0) + &tag_to_string(&couple.1), tag_to_string(&tag));
            res.insert(new_seq);
        }
    }

    res
}

// Reorder the sequences A, B, C, D in any way
pub fn equivalent_reorder(seq : &QuadSeq, seqtype : SequenceType) -> HashSet<QuadSeq> {
    // computes all equivalent sequences by reordering, one swap at a time
    let mut res : HashSet<QuadSeq> = HashSet::new();
    res.insert(seq.clone());

    let couples = [(SequenceTag::W, SequenceTag::X), (SequenceTag::W, SequenceTag::Y), (SequenceTag::W, SequenceTag::Z), (SequenceTag::X, SequenceTag::Y), (SequenceTag::X, SequenceTag::Z), (SequenceTag::Y, SequenceTag::Z)];

    for couple in couples {
        let mut new_seq = seq.clone();
        
        swap(&mut new_seq, couple.0, couple.1);

        debug_assert!(new_seq.verify(seqtype.clone()), "equivalent_reorder function produced invalid {}", seqtype.to_string());
        res.insert(new_seq);
    }

    res
}

// Swap any two pairs of A, B, C, D
pub fn equivalent_double_reorder(seq : &QuadSeq, seqtype : SequenceType) -> HashSet<QuadSeq> {
    // computes all equivalent sequences by reordering, two swaps at a time

    let mut res : HashSet<QuadSeq> = HashSet::new();
    res.insert(seq.clone());

    let couples = [(SequenceTag::W, SequenceTag::X), (SequenceTag::W, SequenceTag::Y), (SequenceTag::W, SequenceTag::Z), (SequenceTag::X, SequenceTag::Y), (SequenceTag::X, SequenceTag::Z), (SequenceTag::Y, SequenceTag::Z)];

    for (couple1, couple2) in iproduct!(couples.clone(), couples) {
        let mut new_seq = seq.clone();
        
        let (seq11, seq12) = couple1;
        let (seq21, seq22) = couple2;
        if !(seq11 == seq21 && seq12 == seq22){
            swap(&mut new_seq, seq11, seq12);
            swap(&mut new_seq, seq21, seq22);
        }

        debug_assert!(new_seq.verify(seqtype.clone()), "equivalent_double_reorder function produced invalid {}", seqtype.to_string());

        res.insert(new_seq);
    
    }

    res
}

// If n is even, cyclically shift all the entries in any of A, B, C, or D by an offset of n/2
pub fn equivalent_uniform_half_shift(seq : &QuadSeq, seqtype : SequenceType) -> HashSet<QuadSeq> {
    // computes equivalent sequences by shifting a single sequence by half length (assuming sequence is even length)
    let mut res : HashSet<QuadSeq> = HashSet::new();
    res.insert(seq.clone());

    if seq.size() % 2 == 1 {
        return res;
    }

    let offset = seq.size() / 2;
    for tag in [SequenceTag::W, SequenceTag::X, SequenceTag::Y, SequenceTag::Z] {
        let mut s = seq.clone();
        let tag_seq = seq.sequence(tag.clone());

        for index in 0..seq.size() {
            s.set_single_value(tag_seq[(index + offset) % seq.size()], &tag, index);
        }

        debug_assert!(s.verify(seqtype.clone()), "equivalent_uniform_half_shift function produced invalid {}", seqtype.to_string());
        res.insert(s);
    }


    res
}

// Cyclically shift all the entries of A, B, C, and D simultaneously by any amount
pub fn equivalent_uniform_shift(seq : &QuadSeq, seqtype : SequenceType) -> HashSet<QuadSeq> {
    // computes all equivalent sequences by shift

    let mut res : HashSet<QuadSeq> = HashSet::new();
    res.insert(seq.clone());

    for offset in 1..seq.size() {
        let mut s = QuadSeq::new(seq.size());
        for index in 0..seq.size() {
            s.set_sequence_value(&seq.values((index + offset) % seq.size()), index)
        }

        debug_assert!(s.verify(seqtype.clone()), "equivalent_uniform_shift function produced invalid {}", seqtype.to_string());
        res.insert(s);
    }

    res
}

// Reverse every sequence simultaneously
pub fn equivalent_reverse(seq : &QuadSeq, seqtype : SequenceType) -> HashSet<QuadSeq> {
    // computes all equivalent sequences by shift
    let mut res : HashSet<QuadSeq> = HashSet::new();
    res.insert(seq.clone());

    let mut s = QuadSeq::new(seq.size());
    for index in 0..seq.size() {
        s.set_sequence_value(&seq.values(seq.size() - 1 - index), index)
    }
    debug_assert!(s.verify(seqtype.clone()), "equivalent_reverse function produced invalid {}", seqtype.to_string());
    res.insert(s);

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

// Negate all the entries of any of A, B, C, or D
pub fn equivalent_negate(seq : &QuadSeq, seqtype : SequenceType) -> HashSet<QuadSeq> {
    // computes all equivalent sequences by negation of any sequence
    let (a,b,c,d) = seq.sequences();
    let (nega_a,nega_b,nega_c,nega_d) = (negated(&a), negated(&b), negated(&c), negated(&d));

    let mut res : HashSet<QuadSeq> = HashSet::new();
    res.insert(seq.clone());

    for tag in [SequenceTag::W, SequenceTag::X, SequenceTag::Y, SequenceTag::Z] {
        let quad = match tag {
            SequenceTag::W => (&nega_a, &b, &c, &d),
            SequenceTag::X => (&a, &nega_b, &c, &d),
            SequenceTag::Y => (&a, &b, &nega_c, &d),
            SequenceTag::Z => (&a, &b, &c, &nega_d)
        };
        let mut new_seq = QuadSeq::new(seq.size());
        new_seq.set_all_values(quad);

        debug_assert!(new_seq.verify(seqtype.clone()), "equivalent_negate function produced invalid {}", seqtype.to_string());
        res.insert(new_seq);
    }

    res
}

// Negate any two sequences of A, B, C, or D
pub fn equivalent_double_negate(seq : &QuadSeq, seqtype : SequenceType) -> HashSet<QuadSeq> {
    // computes all equivalent sequences by negation of two sequences

    let (a,b,c,d) = seq.sequences();
    let (nega_a,nega_b,nega_c,nega_d) = (negated(&a), negated(&b), negated(&c), negated(&d));

    let mut res : HashSet<QuadSeq> = HashSet::new();
    res.insert(seq.clone());

    for tag_couple in [(SequenceTag::W, SequenceTag::X), (SequenceTag::W, SequenceTag::Y), (SequenceTag::W, SequenceTag::Z), (SequenceTag::X, SequenceTag::Y), (SequenceTag::X, SequenceTag::Z), (SequenceTag::Y, SequenceTag::Z)] {
        // this loops through all the couples of a,b,c,d (ordered couples)
        let quad = match tag_couple {
            (SequenceTag::W, SequenceTag::X) => {(&nega_a, &nega_b, &c, &d)},
            (SequenceTag::W, SequenceTag::Y) => {(&nega_a, &b, &nega_c, &d)},
            (SequenceTag::W, SequenceTag::Z) => {(&nega_a, &b, &c, &nega_d)},
            (SequenceTag::X, SequenceTag::Y) => {(&a, &nega_b, &nega_c, &d)},
            (SequenceTag::X, SequenceTag::Z) => {(&a, &nega_b, &c, &nega_d)},
            (SequenceTag::Y, SequenceTag::Z) => {(&a, &b, &nega_c, &nega_d)},
            _ => {panic!("Incorrect tags entered !")}
        };
        let mut s = QuadSeq::new(seq.size());
        s.set_all_values(quad);

        debug_assert!(s.verify(seqtype.clone()), "equivalent_double_negate function produced invalid {}. Original: {}\nNew: {}", seqtype.to_string(), seq.to_string(), s.to_string());
        res.insert(s);
    }

    res
}

// Multiply every other element of A, B, C, and D by -1 simultaneously
pub fn equivalent_alternated_negation(seq : &QuadSeq, seqtype : SequenceType) -> HashSet<QuadSeq> {
    let frequency = 2;

    let (a,b,c,d) = seq.sequences();

    let mut res : HashSet<QuadSeq> = HashSet::new();
    res.insert(seq.clone());

    let quads = (&alt_negated(&a, frequency), &alt_negated(&b, frequency), &alt_negated(&c, frequency), &alt_negated(&d, frequency));

    let mut s = QuadSeq::new(seq.size());
    s.set_all_values(quads);

    debug_assert!(s.verify(seqtype.clone()), "equivalent_alternated_negation function produced invalid {}", seqtype.to_string());
    res.insert(seq.clone());
    res.insert(s);

    res
}

// If n is even, multiply every other element of A, B, C, and D by -1 simultaneously
pub fn equivalent_even_alternated_negation(seq : &QuadSeq, seqtype : SequenceType) -> HashSet<QuadSeq> {
    let mut res : HashSet<QuadSeq> = HashSet::new();
    res.insert(seq.clone());

    if seq.size() % 2 == 1 {
        return res;
    }

    let frequency = 2;

    let (a,b,c,d) = seq.sequences();

    let quads = (&alt_negated(&a, frequency), &alt_negated(&b, frequency), &alt_negated(&c, frequency), &alt_negated(&d, frequency));

    let mut s = QuadSeq::new(seq.size());
    s.set_all_values(quads);

    debug_assert!(s.verify(seqtype.clone()), "equivalent_even_alternated_negation function produced invalid {}", seqtype.to_string());
    res.insert(seq.clone());
    res.insert(s);

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




// Apply an automorphism of the cyclic group C_n to all the indices of the entries of each of A, B, C, and D simultaneously
pub fn equivalent_automorphism(seq : &QuadSeq, seqtype : SequenceType) -> HashSet<QuadSeq> {
    // computes all equivalent sequences by permutation

    let mut res : HashSet<QuadSeq> = HashSet::new();
    res.insert(seq.clone());
    let size = seq.size();

    for k in COPRIMES[size].iter() {

        let mut will = QuadSeq::new(size);

        let (a,b,c,d) = seq.sequences();

        let quad = (&permute(&a, *k), &permute(&b, *k), &permute(&c, *k), &permute(&d, *k));

        will.set_all_values(quad);

        debug_assert!(will.verify(seqtype.clone()), "equivalent_automorphism function produced invalid {}", seqtype.to_string());
        res.insert(will);

    }


    res
}
