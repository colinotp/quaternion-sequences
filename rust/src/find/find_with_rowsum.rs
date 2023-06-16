use std::usize::MAX;

use itertools::Itertools;

use crate::sequences::{rowsum::{generate_rowsums, Quad, generate_sequences_with_rowsum}, williamson::SequenceTag, fourier::{iter_over_filtered_dft}, self};


fn get_two_best(quad: &Quad) -> ((usize, usize),(usize, usize)){
    let (mut maxi, mut index) = (quad.0, 0);
    if quad.1 < maxi {(maxi, index) = (quad.1, 1)}
    if quad.2 < maxi {(maxi, index) = (quad.2, 2)}
    if quad.3 < maxi {(maxi, index) = (quad.3, 3)}

    let (mut maxi2, mut index2) = (MAX, 4);
    if quad.0 < maxi2 && index != 0 {(maxi2, index2) = (quad.1, 1)}
    if quad.1 < maxi2 && index != 1 {(maxi2, index2) = (quad.1, 1)}
    if quad.2 < maxi2 && index != 2 {(maxi2, index2) = (quad.2, 2)}
    if quad.3 < maxi2 && index != 3 {(maxi2, index2) = (quad.3, 3)}
    ((maxi,index),(maxi2,index2))
}


pub fn find(p : usize) {

    let rowsums = generate_rowsums(4*p);

    for rs in rowsums {
        let ((maxi,index),(maxi2,index2)) = get_two_best(&rs);

        let sequences_1 = generate_sequences_with_rowsum(maxi, p);
        let sequences_2 = generate_sequences_with_rowsum(maxi2, p);

        for pair1 in iter_over_filtered_dft(&sequences_1, &(4.*p as f64)){
            for pair2 in iter_over_filtered_dft(&sequences_2, &(4.*p as f64)){
                // convert to equations
    
                // solve using libexact
            }
        }
    }

}