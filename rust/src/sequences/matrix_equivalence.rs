use std::collections::HashSet;

use itertools::Itertools;

use super::matrices::HM;



pub fn generate_equivalence_class(hm : &HM) -> HashSet<HM> {
    // This function generates the representative of the equivalence class that seq belongs to
    
    let mut class = HashSet::new();
    class.insert(hm.clone());

    loop {
        let mut new = HashSet::new();

        for seq in &class {
            for equivalence in [equivalent_negate_col, equivalent_negate_row, equivalent_swapping_col, equivalent_swapping_row] {
                for equ in equivalence(&seq){
                    if !class.contains(&equ) {
                        new.insert(equ);
                    }
                }
            }

        }

        if new.len() == 0 {break;}
        else {
            println!("{} matrices added", new.len());
            for seq in new {
                class.insert(seq);
            }
        }
    }

    class
}




pub fn equivalent_swapping_row(hm : &HM) -> Vec<HM> {

    let mut res = vec![];

    for perm in (0..hm.size()).permutations(2) {
        let (i,j) = (perm[0], perm[1]);
        let mut new_hm = hm.clone();
        for index in 0..hm.size() {
            new_hm.set_value(i, index, hm.get(j, index));
            new_hm.set_value(j, index, hm.get(i, index));
        }
        res.push(new_hm);
    }

    res
}

pub fn equivalent_swapping_col(hm : &HM) -> Vec<HM> {
    
    let mut res = vec![];

    for perm in (0..hm.size()).permutations(2) {
        let (i,j) = (perm[0], perm[1]);
        let mut new_hm = hm.clone();
        for index in 0..hm.size() {
            new_hm.set_value(index, i, hm.get(j, index));
            new_hm.set_value(index, j, hm.get(i, index));
        }
        res.push(new_hm);
    }

    res
}

pub fn equivalent_negate_row(hm : &HM) -> Vec<HM> {

    let mut res = vec![];

    for row in 0..hm.size() {
        let mut new_hm = hm.clone();
        for index in 0..hm.size() {
            new_hm.set_value(row, index, -hm.get(row, index));
        }
        res.push(new_hm);
    }

    res
}

pub fn equivalent_negate_col(hm : &HM) -> Vec<HM> {

    let mut res = vec![];

    for col in 0..hm.size() {
        let mut new_hm = hm.clone();
        for index in 0..hm.size() {
            new_hm.set_value(index, col, -hm.get(index, col));
        }
        res.push(new_hm);
    }

    res
}