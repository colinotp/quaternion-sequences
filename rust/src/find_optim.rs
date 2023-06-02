
use rayon;
use rayon::prelude::{IntoParallelRefIterator, ParallelIterator};
use crate::sequence::*;
use crate::symmetries::Symmetry;


pub fn find(size : usize, symmetry : Option<Symmetry>) -> usize{
    let mut pqs = QS::new(size, symmetry);

    return find_recursive(&mut pqs, size, 1);
}

fn find_recursive(pqs : &mut QS, size : usize, index : usize) -> usize{

    if index >= pqs.search_size(){
        if pqs.is_perfect(){
            //println!("{}", pqs.to_string());
            return 1;
        }
        return 0;
    }

    return QPLUS.par_iter() // we use parallel iterators
        .map(|q| {
            let mut new_pqs = pqs.clone(); // we take a new sequence, modify it, and repeat recursively
            new_pqs.set_value(q.clone(), index);
            find_recursive(&mut new_pqs, size, index+1)
        })
        .fold(|| 0, |a,b| a+b)
        .reduce(|| 0, |a,b| a+b);
}


