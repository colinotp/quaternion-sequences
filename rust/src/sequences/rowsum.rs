use itertools::*;

use crate::sequences::symmetries::SequenceType;


// * Sequence generation with specific rowsum

pub type Quad = (isize,isize,isize,isize);

pub fn quad_contains_zero(quad : &Quad) -> bool {
    (quad.0 == 0) || (quad.1 == 0) || (quad.2 == 0) || (quad.3 == 0)
}

pub fn quad_contains_dup(quad : &Quad) -> bool {
    (quad.0 == quad.1) || (quad.0 == quad.2) || (quad.0 == quad.3) || (quad.1 == quad.2) || (quad.1 == quad.3) || (quad.2 == quad.3)
}

pub fn sequence_to_string(seq: &Vec<i8>) -> String {
    let mut res_str: String = "".to_owned();

    for i in seq.iter(){
        let str = match i {
            1 => "+",
            -1 => "-",
            _ => {panic!("not a +-1 sequence !")}
        };
        res_str.push_str(&str);
    }

    res_str.to_string()
}



pub fn rowsum(seq : Vec<i8>) -> isize {
    // computes the rowsum of the sequence seq
    seq.iter().map(|i| *i as isize).sum()
}



pub fn generate_sequences_with_rowsum(rowsum: isize, size : usize) -> Vec<Vec<i8>> {
    // generates all sequences of length size and whose sum equals rowsum

    if (rowsum.abs() % 2) as usize != size % 2 {
        // no combination will work
        return vec![];
    }

    // we get the number of ones in our sequence
    let nb_ones = (size as isize + rowsum)/2;

    // and we call the recursive function
    let mut seq : Vec<i8> = vec![-1;size];
    gen_seq_rec(&mut seq, nb_ones as usize, 0)
}


pub fn gen_seq_rec(seq : &mut Vec<i8>, remaining_ones : usize, current_pos : usize) -> Vec<Vec<i8>> {

    if remaining_ones == 0 {
        // We're done ! There are no more ones to place
        return vec![seq.clone()];
    }
    if current_pos + remaining_ones > seq.len() {
        // We can't possibly fit the remaining ones in the rest of the sequence, so this is impossible
        return vec![];
    }

    // Either there's a -1 in position current_pos...
    let mut results2 = gen_seq_rec(seq, remaining_ones, current_pos + 1);
    
    // or there's not.
    seq[current_pos] = 1;
    let mut results1 = gen_seq_rec(seq, remaining_ones - 1, current_pos + 1);

    seq[current_pos] = -1; // We reset the sequence to it's original state

    // We concatenate both results
    results1.append(&mut results2);

    results1
}



// * Rowsum generation


fn square_sum(s : &Quad) -> usize {
    (s.0*s.0 + s.1*s.1 + s.2*s.2 + s.3*s.3) as usize
}

fn increment_squares(s : &mut Quad, bound : usize) {
    // auxiliary function of sum_of_four_squares.
    // calling this function transforms s to the next quadruplet to test

    s.3 += 1;
    if square_sum(&s) <= bound {
        return;
    }

    s.2 += 1;
    s.3 = s.2;
    if square_sum(&s) <= bound {
        return;
    }
    
    s.1 += 1;
    s.2 = s.1;
    s.3 = s.1;
    if square_sum(&s) <= bound {
        return;
    }

    s.0 += 1;
    s.1 = s.0;
    s.2 = s.0;
    s.3 = s.0;
}



pub fn sum_of_four_squares(p : usize) -> Vec<Quad> {
    // returns all possible combination of quadruplets of integers whose square sum is p
    // each quadruplet is sorted in increasing order, so there are no permutations.

    let mut squares_list : Vec<Quad> = vec![];
    let mut squares : Quad = (0,0,0,0);

    while (squares.0*squares.0) as usize <= p {
        if square_sum(&squares) == p {
            squares_list.push(squares.clone())
        }

        increment_squares(&mut squares, p);
    }

    squares_list
}




pub fn generate_other_quadruplets(quad : &Quad, seqtype : SequenceType) -> Vec<Quad> {
    // returns all the permutations of the quadruplet up to equivalence
    
    // generate all permutations and filter them
    let quad_vec = vec![quad.0,quad.1,quad.2,quad.3];

    let mut isomorphism = vec![];

    for perm in quad_vec.iter().permutations(4) {
        let quad = (*perm[0], *perm[1], *perm[2], *perm[3]);
        isomorphism.push(quad);
    }

    // filtering and adding negated elements
    let mut result = vec![];
    for elm in isomorphism.into_iter().unique_by(|q| equivalent(&q)).collect::<Vec<Quad>>() {
        result.push(elm.clone());
        if matches!(seqtype, SequenceType::QuaternionType) && !quad_contains_zero(&elm) && !quad_contains_dup(&elm) {
            let mut nega_elm = elm.clone();
            nega_elm.0 = - nega_elm.0;
            result.push(nega_elm)
        }
    }

    result.into_iter().unique().collect()
}




fn swap(quad : &mut Quad, i1 : usize, i2 : usize) {
    // swaps two indices of a Quad
    let val1 = match i1 {
        0 => {quad.0},
        1 => {quad.1},
        2 => {quad.2},
        3 => {quad.3},
        _ => {panic!()}
    };
    let val2 = match i2 {
        0 => {quad.0},
        1 => {quad.1},
        2 => {quad.2},
        3 => {quad.3},
        _ => {panic!()}
    };
    match i1 {
        0 => {quad.0 = val2},
        1 => {quad.1 = val2},
        2 => {quad.2 = val2},
        3 => {quad.3 = val2},
        _ => {panic!()}
    };
    match i2 {
        0 => {quad.0 = val1},
        1 => {quad.1 = val1},
        2 => {quad.2 = val1},
        3 => {quad.3 = val1},
        _ => {panic!()}
    };

}

fn better_than(q1 : &Quad, q2 : &Quad) -> bool {
    // compares two quads
    if q1.0 < q2.0 {true}
    else if q1.0 > q2.0 {false}
    else{
        if q1.1 < q2.1 {true}
        else if q1.1 > q2.1 {false}
        else{
            if q1.2 < q2.2 {true}
            else if q1.2 > q2.2 {false}
            else{
                if q1.3 < q2.3 {true}
                else{false}
            }
        }
    }
}

fn equivalent(quad : &Quad) -> Quad {
    // finds the representative of the equivalence class that quad belongs to
    let mut final_quad = quad.clone();

    for swaps in [(0,1,2,3), (0,2,1,3), (0,3,1,2), (0,1,1,2), (0,1,1,3), (0,2,1,2), (0,2,2,3), (0,3,1,3), (0,3,2,3), (1,2,2,3), (1,3,2,3)] {
        let mut test_quad = quad.clone(); // at the start, this is 1 2 3 4
        swap(&mut test_quad, swaps.0, swaps.1); 
        swap(&mut test_quad, swaps.2, swaps.3); // 2 1 4 3
        if better_than(&test_quad, &final_quad) {final_quad = test_quad.clone()}
    }

    final_quad
}




pub fn generate_rowsums(p : usize, seqtype : SequenceType) -> Vec<Quad>{
    // generates all quadruplets of integers such that the sum of their squares equal 4*p
    // it also generates their permutation, but only up to equivalence
    let quads = sum_of_four_squares(4*p);

    let mut total_quadruplets = vec![];

    let parity = (p % 2) as isize;

    for elm in quads {
        if parity == elm.0 % 2 && parity == elm.1 % 2 && parity == elm.2 % 2 && parity == elm.3 % 2 {
            total_quadruplets.append(&mut generate_other_quadruplets(&elm, seqtype.clone()));
        }
    }

    total_quadruplets
}