use std::{isize::MIN, fs::{*, self}, path::Path, io::Write, env};


use crate::sequences::{rowsum::{generate_rowsums, Quad, generate_sequences_with_rowsum, sequence_to_string}, fourier::{iter_over_filtered_dft, iter_over_filtered_couples}, equations::generate_equations, williamson::{SequenceTag, tag_to_string, Williamson}, symmetries::SequenceType, matching::{generate_matching_table, MatchData, compute_complementary_auto_correlations, compute_complementary_cross_correlations, verify_cross_correlation}};


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
        0 => {SequenceTag::X}
        1 => {SequenceTag::Y}
        2 => {SequenceTag::Z}
        3 => {SequenceTag::W}
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



pub fn sort(quad : &Quad) -> (Vec<isize>, Vec<usize>){

    let mut tab = vec![quad.0, quad.1, quad.2, quad.3];
    let mut indices = vec![0,1,2,3];

    for i in 0..4 {
        
        let (mut maxi, mut index) = (tab[i], i);
        for j in i+1..4 {
            if tab[j] > maxi {(maxi, index) = (tab[j], j)}
        }
        (tab[i], tab[index]) = (tab[index], tab[i]);
        (indices[i], indices[index]) = (indices[index], indices[i]);
    }

    (tab, indices)
}


pub fn find_matching(p : usize) -> Vec<Williamson>{

    // The resulting list of sequences
    let mut matches = vec![];

    // all the possible rowsums of p
    let rowsums = generate_rowsums(p);
    for rs in &rowsums {
        println!("{:?}", rs);
    }
    println!("generated {} different rowsums", rowsums.len());

    for rs in rowsums {
        let (rowsums, indices) = sort(&rs); // we sort the rowsum in decreasing order, and we keep track of their original indices
        let tags : Vec<SequenceTag> = indices.iter().map(|i| index_to_tag(*i)).collect(); // we convert the indices to their respective tags

        // We generate all the sequences possible for each rowsums
        let sequences_1 = generate_sequences_with_rowsum(rowsums[0], p);
        let sequences_2 = generate_sequences_with_rowsum(rowsums[1], p);
        let sequences_3 = generate_sequences_with_rowsum(rowsums[2], p);
        let sequences_4 = generate_sequences_with_rowsum(rowsums[3], p);

        // We create the matching table : all the sequences that give the same auto-correlation and cross-correlation values are put in the same entry, in a list
        let match_table = generate_matching_table(&sequences_3, &sequences_4, &(tags[2].clone(), tags[3].clone()), p);

        for (seq1, seq2) in iter_over_filtered_couples(&sequences_1, &sequences_2, 4.*p as f64) {
            // We iterate over the couples of sequences, but we filter out some with the dft checks

            // We compute the auto and cross correlation values when considered on the other side of the equation
            let autoc_values = compute_complementary_auto_correlations(seq1, seq2, p);
            let crossc_values = compute_complementary_cross_correlations(seq1, seq2, &(tags[0].clone(), tags[1].clone()));

            // We construct a new MatchData to find all the sequences that match using the Match Table
            let match_data = MatchData::new(autoc_values, crossc_values);
            let possible_matches = match_table.get(&match_data);

            if possible_matches == None {
                continue; // There are no matches found
            }

            for (ref_seq3, ref_seq4) in possible_matches.unwrap() {
                // If we found matches, we test them all
                let (seq3, seq4) = (*ref_seq3, *ref_seq4); // The sequences are behind two references

                if verify_cross_correlation(&[&seq1, &seq2, seq3, seq4], &tags) {

                    let mut will = Williamson::new(p);
                    will.set_all_values(&(seq1.clone(), seq2.clone(), seq3.clone(), seq4.clone()));

                    matches.push(will)
                }
            }

        }
    }

    matches
}