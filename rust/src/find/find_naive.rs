
use crate::sequences::sequence::*;
use crate::sequences::symmetries::Symmetry;

#[allow(dead_code)]
pub fn find(size : usize, symmetry : Option<Symmetry>) -> usize{
    let mut count = 0;
    let mut pqs = QS::new(size, symmetry);

    find_recursive(&mut pqs, size, 1, &mut count);

    return count;
}

#[allow(dead_code)]
fn find_recursive(pqs : &mut QS, size : usize, index : usize, count: &mut usize){

    if index >= pqs.search_size(){
        //println!("{}", pqs.to_string());
        if pqs.is_perfect(){
            *count+=1;
        }
        return;
    }

    for value_to_test in QPLUS.iter(){ // tries every element possible recursively
        pqs.set_value(value_to_test.clone(), index);

        find_recursive(pqs, size, index+1, count);
    }
}


