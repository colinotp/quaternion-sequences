
use crate::williamson::{Williamson, QUADRUPLETS};




pub fn find(size : usize) -> usize{
    let mut will = Williamson::new(size);

    return find_recursive(&mut will, size, 1);
}

fn find_recursive(will : &mut Williamson, size : usize, index : usize) -> usize{

    if index >= will.search_size(){
        if will.is_perfect_complementary(){
            if !will.to_qs().is_pqs(){
                println!("{}", will.to_string());
                println!("{}", will.to_qs().to_string());
            }
            return 1;
        }
        return 0;
    }

    let mut count = 0;
    for value_to_test in QUADRUPLETS.iter(){
        let mut will1 = will.clone();
        will1.set_sequence_value(value_to_test, index);
        count += find_recursive(&mut will1, size, index+1);
    }

    count
}


