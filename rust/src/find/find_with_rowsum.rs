use std::{isize::MIN, fs::{*, self}, path::Path, io::Write, env};


use crate::sequences::{rowsum::{generate_rowsums, Quad, generate_sequences_with_rowsum, sequence_to_string}, fourier::{iter_over_filtered_dft}, equations::generate_equations, williamson::{SequenceTag, tag_to_string}, symmetries::SequenceType};


fn get_two_best(quad: &Quad) -> ((isize, usize),(isize, usize)){
    let (mut maxi, mut index) = (quad.0, 0);
    if quad.1 > maxi {(maxi, index) = (quad.1, 1)}
    if quad.2 > maxi {(maxi, index) = (quad.2, 2)}
    if quad.3 > maxi {(maxi, index) = (quad.3, 3)}

    let (mut maxi2, mut index2) = (MIN, 4);
    if quad.0 >= maxi2 && index != 0 {(maxi2, index2) = (quad.0, 0)}
    if quad.1 >= maxi2 && index != 1 {(maxi2, index2) = (quad.1, 1)}
    if quad.2 >= maxi2 && index != 2 {(maxi2, index2) = (quad.2, 2)}
    if quad.3 >= maxi2 && index != 3 {(maxi2, index2) = (quad.3, 3)}

    if index < index2 {
        ((maxi,index),(maxi2,index2))
    }
    else{
        ((maxi2,index2),(maxi,index))
    }
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

fn generate_comment(seq1 : &Vec<i8>, tag1 : &SequenceTag, seq2 : &Vec<i8>, tag2 : &SequenceTag, rowsum : &Quad) -> String {
    let mut comment = "* ".to_string() + &tag_to_string(tag1) + &": " + &sequence_to_string(seq1) + &"\n* " + &tag_to_string(tag2) + &": " + &sequence_to_string(seq2) + &"\n";
    comment += &("* rowsum: (".to_string() + &rowsum.0.to_string() + &"," + &rowsum.1.to_string() + &"," + &rowsum.2.to_string() + &"," + &rowsum.3.to_string() + &")");
    comment
}




pub fn find(p : usize, seqtype : SequenceType) {

    let rowsums = generate_rowsums(p);
    println!("generated {} different rowsums", rowsums.len());
    println!("Current directory: {:?}", env::current_dir().ok().unwrap());

    let folder = match seqtype {
        SequenceType::WilliamsonType => {"wts"}
        _ => {panic!("not implemented yet")} // TODO
    };

    for rs in rowsums {
        let ((maxi,index),(maxi2,index2)) = get_two_best(&rs); // should we get the two best or just the same two each time ?
        let (tag1, tag2) = (index_to_tag(index), index_to_tag(index2));

        let sequences_1 = generate_sequences_with_rowsum(maxi, p);
        let sequences_2 = generate_sequences_with_rowsum(maxi2, p);
        let string_path = "results/equations/".to_string()+ &folder + &"/find_" + &p.to_string() + &"/rowsum_" + &maxi.to_string() + &"-" + &maxi2.to_string() + &"_at_" + &index.to_string() + &"-" + &index2.to_string();

        println!("{}",string_path);
        fs::create_dir_all(&string_path).expect("Error when creating the dir");

        let mut count = 0;
        for seq1 in iter_over_filtered_dft(&sequences_1, 4.*p as f64){
            println!("loop");
            for seq2 in iter_over_filtered_dft(&sequences_2, 4.*p as f64){
                let s = string_path.clone() + &"/num_" + &count.to_string() + &".opb";
                println!("{}",s);
                let path = Path::new(&s);

                // convert to equations
                let equations = generate_equations(&seq1, &tag1, &seq2, &tag2, &seqtype, &rs);
                // TODO add rowsum constraints

                //if there are equations to write, create a new file
                if equations != "" {
                    let mut f = File::create(path).expect("Invalid file ?");
                    let comment = generate_comment(seq1, &tag1, seq2, &tag2, &rs);
                    f.write(comment.as_bytes()).expect("Error when writing in the file");
                    f.write(equations.as_bytes()).expect("Error when writing in the file");
                    count +=1;
                }

                // solve using libexact


            }
        }
    }

}