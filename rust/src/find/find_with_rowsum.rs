use std::{usize::MAX, fs::{*, self}, path::Path, io::Write, env};


use crate::sequences::{rowsum::{generate_rowsums, Quad, generate_sequences_with_rowsum}, fourier::{iter_over_filtered_dft}, equations::generate_equations, williamson::SequenceTag, symmetries::SequenceType};


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


fn index_to_tag(index: usize) -> SequenceTag {
    match index {
        0 => {SequenceTag::A}
        1 => {SequenceTag::B}
        2 => {SequenceTag::C}
        3 => {SequenceTag::D}
        _ => {panic!("incorrect index")}
    }
}




pub fn find(p : usize, seqtype : SequenceType) {

    let rowsums = generate_rowsums(4*p);
    println!("generated {} different rowsums", rowsums.len());
    println!("Current directory: {:?}", env::current_dir().ok());

    for rs in rowsums {
        let ((maxi,index),(maxi2,index2)) = get_two_best(&rs); // should we get the two best or just the same two each time ?
        let (tag1, tag2) = (index_to_tag(index), index_to_tag(index2));

        let sequences_1 = generate_sequences_with_rowsum(maxi, p);
        let sequences_2 = generate_sequences_with_rowsum(maxi2, p);
        let string_path = "results/equations/rowsum_".to_string() + &maxi.to_string() + &"-" + &maxi2.to_string() + &"_at_" + &index.to_string() + &"-" + &index2.to_string();

        println!("{}",string_path);
        fs::create_dir_all(&string_path).expect("Error when creating the dir");

        let mut count = 0;
        for seq1 in iter_over_filtered_dft(&sequences_1, &(4.*p as f64)){
            for seq2 in iter_over_filtered_dft(&sequences_2, &(4.*p as f64)){
                let s = string_path.clone() + &"/num_" + &count.to_string() + &".opb";
                println!("{}",s);
                let path = Path::new(&s);
                let mut f = File::create(path).expect("Invalid file ?");

                // convert to equations
                let equations = generate_equations(&seq1, &tag1, &seq2, &tag2, &seqtype);
                f.write(equations.as_bytes()).expect("Error when writing in the file");

                // solve using libexact


                count +=1;
            }
        }
    }

}