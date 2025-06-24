use std::{time::Instant, fs::{self, File, DirEntry}, io::Write, io::Error};
use itertools::iproduct;
use memory_stats::memory_stats;

use crate::{find::find_unique::reduce_to_equivalence, read_lines, sequences::{fourier::iter_over_enumerate_filtered_couples, matching::{compute_auto_correlations, compute_cross_correlations}, rowsum::{generate_rowsums, generate_sequences_with_rowsum, rowsum, sequence_to_string, Quad}, symmetries::*, williamson::{tag_to_string, QuadSeq, SequenceTag}}, str_to_seqtype};



pub fn quad_to_string(q : (isize, isize, isize, isize)) -> String {
    // Transforms a quadruplet of integers into a string

    let (a,b,c,d) = q;
    a.to_string() + &" " + &b.to_string() + &" " + &c.to_string() + &" " + &d.to_string() + &"\n"
}

pub fn write_rowsums(p : usize, seqtype : SequenceType) {
    // Stores the possible rowsums for qts sequences of length p
    let folder = seqtype.to_string();

    let folder_path = "results/pairs/".to_string()+ &folder + &"/find_" + &p.to_string();
    let path = folder_path.clone() + &"/rowsums.quad";
    let mut f = File::create(path).expect("Invalid file ?");

    let rs = generate_rowsums(p);

    let s = rs.iter().map(|e| quad_to_string(*e)).fold("".to_string(), |a,b| a + &b);

    f.write(s.as_bytes()).expect("Error writing file");
}




pub enum EquationSide {
    LEFT, RIGHT
}


pub fn sort(quad : &Quad) -> (Vec<isize>, Vec<usize>){
    // Sorts a quadruplet of integers

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

pub fn print_memory_usage(message : &str) {
    eprintln!("{}",message);
    if let Some(usage) = memory_stats() {
        eprintln!("physical memory usage: {} bytes", usage.physical_mem);
        eprintln!("virtual memory usage: {} bytes ", usage.virtual_mem);
    } else {
        eprintln!("Couldn't get the current memory usage :(");
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

pub fn write_sequences(sequences : &Vec<Vec<i8>>, tag : &SequenceTag, folder_path : &String) {
    // stores the sequences

    let path = folder_path.clone() + &"/seq_" + &tag_to_string(tag) + ".seq";
    let mut f = File::create(path).expect("Invalid file ?");

    for seq in sequences {
        f.write((sequence_to_string(seq) + &"\n").as_bytes()).expect("Error when writing in the file");
    }
}

pub fn verify_rowsums(sequences : (&Vec<Vec<i8>>, &Vec<Vec<i8>>), tags : (&SequenceTag, &SequenceTag), rs : (isize, isize, isize, isize)) -> bool {
    let rowsum_0 : isize = match tags.0 {
            SequenceTag::X => rs.0,
            SequenceTag::Y => rs.1,
            SequenceTag::Z => rs.2,
            SequenceTag::W => rs.3
        };
    

    for seq in sequences.0 {
        if rowsum(seq.clone()) != rowsum_0 {
            return false;
        }
    }

    let rowsum_1 : isize = match tags.1 {
            SequenceTag::X => rs.0,
            SequenceTag::Y => rs.1,
            SequenceTag::Z => rs.2,
            SequenceTag::W => rs.3
        };
    

    for seq in sequences.1 {
        if rowsum(seq.clone()) != rowsum_1 {
            return false;
        }
    }
    
    true
}

pub fn write_seq_pairs(sequences : (&Vec<Vec<i8>>, &Vec<Vec<i8>>), tags : (&SequenceTag, &SequenceTag), seqtype : SequenceType, rs : (isize, isize, isize, isize), p : usize, folder_path : &String, side : EquationSide) {
    // This function generates the files that end in .pair used for the algorithm

    assert!(verify_rowsums(sequences, tags, rs));

    let path = folder_path.clone() + &"/pair_" + &tag_to_string(&tags.0) + &tag_to_string(&tags.1) + ".pair";
    let mut f = File::create(path).expect("Invalid file ?");

    let op = match side {
        EquationSide::LEFT => {|x : isize| x}
        EquationSide::RIGHT => {|x : isize| -x}
    };

    // Instead of writing each line one by one n the file, we use a buffer to write them by chunks of 1000 lines
    let mut buffer = "".to_string();
    let mut buffer_counter = 0;

    // We iterate over the couples of sequences, but we filter out some with the dft checks
    for ((index0, seq0), (index1, seq1)) in iter_over_enumerate_filtered_couples(sequences.0, sequences.1, 4.*p as f64){
        let mut result = "".to_string();

        // We compute the auto and cross correlation values when considered on the other side of the equation
        let autoc_values = compute_auto_correlations(seq0, seq1);
        let crossc_values = compute_cross_correlations(seq0, seq1, &(tags.0.clone(), tags.1.clone()));
        
        // We add these values to the current line
        for a in autoc_values {
            result += &(op(a).to_string() + &"_");
        }
        
        // In the case of WTS, we want to filter out pairs with nonzero crosscorrelation values
        if matches!(seqtype, SequenceType::WilliamsonType) {
            let mut crossc_result = "".to_string();
            for c in crossc_values {
                if c == 0 {
                    crossc_result = "".to_string();
                    break;
                }
                crossc_result += &(op(c).to_string() + &"_");
                
            }
            result += &crossc_result;
        }
        else {
            for c in crossc_values {
                result += &(op(c).to_string() + &"_");
            }
        }
        

        result += &(":_".to_string() + &index0.to_string() + "_" + &index1.to_string() + &"\n");

        buffer += &result;
        buffer_counter += 1;

        if buffer_counter > 1000 {
            f.write(buffer.as_bytes()).expect("Error when writing in the file");
            buffer = "".to_string();
            buffer_counter = 0;
        }

    }
    
    f.write(buffer.as_bytes()).expect("Error when writing in the file");

}

pub fn get_indices(pairing: Option<RowsumPairing>, pair: u8) -> Option<(usize, usize)> {
    match (pairing, pair) {
        (Some(RowsumPairing::XY), 1) => Some((0, 1)),     // XY
        (Some(RowsumPairing::XY), 2) => Some((2, 3)),     // ZW
        (Some(RowsumPairing::XZ), 1) => Some((0, 2)),     // XZ
        (Some(RowsumPairing::XZ), 2) => Some((1, 3)),     // YW
        (Some(RowsumPairing::XW), 1) => Some((0, 3)),     // XW
        (Some(RowsumPairing::XW), 2) => Some((1, 2)),     // YZ
        _ => None
    }
}

pub fn write_pair_single(seqtype : SequenceType, p: usize, pairing: Option<RowsumPairing>, pair: u8) {
    // This function is identical to write_pairs(), except for the purpose of running pairs individually on separate processors
    // `pair` should be either a 1 or a 2, which decides whether to look at the first or second pair given by the chosen pairing

    // all the possible rowsums of p
    let rowsums = generate_rowsums(p);
    for rs in &rowsums {
        eprintln!("{:?}", rs);
    }
    eprintln!("generated {} different rowsums", rowsums.len());

    let folder = seqtype.to_string();

    for rs in rowsums {
        write_pair_single_rowsum(folder.clone(), rs, p, pairing.clone(), pair);
    }

}

pub fn write_pair_single_rowsum(folder : String, rs : (isize, isize, isize, isize), p : usize, pairing: Option<RowsumPairing>, pair: u8) {
    let (rowsums, indices) = sort(&rs); // we sort the rowsum in decreasing order, and we keep track of their original indices
    let tags : Vec<SequenceTag> = indices.iter().map(|i| index_to_tag(*i)).collect(); // we convert the indices to their respective tags

    let folder_path = "results/pairs/".to_string()+ &folder + &"/find_" + &p.to_string() + &"/rowsum_" + &(rs.0).to_string() + &"_" + &(rs.1).to_string() + &"_" + &(rs.2).to_string() + &"_" + &(rs.3).to_string();
    println!("{}",folder_path);
    fs::create_dir_all(&folder_path).expect("Error when creating the dir");     // This is safe to do concurrently across multiple processes according to the documentation

    let sequences_0: Vec<Vec<i8>>;
    let sequences_1: Vec<Vec<i8>>;

    let pair_indices;
    match get_indices(pairing, pair) {
        Some(s) => {pair_indices = s},
        None => {
            println!("ERROR: get_indices() returned None");
            return;
        }
    }

    let now = Instant::now();
    sequences_0 = generate_sequences_with_rowsum(rowsums[pair_indices.0], p);
    sequences_1 = generate_sequences_with_rowsum(rowsums[pair_indices.1], p);

    write_sequences(&sequences_0, &tags[pair_indices.0], &folder_path);
    write_sequences(&sequences_1, &tags[pair_indices.1], &folder_path);

    let elapsed_time = now.elapsed().as_secs_f32();
    eprintln!("The function took: {elapsed_time} seconds to generate sequences with rowsums: {}, {}, {}, {}", rowsums[0], rowsums[1], rowsums[2], rowsums[3]);

    let side;
    match pair {
        1 => {side = EquationSide::LEFT},
        2 => {side = EquationSide::RIGHT},
        _ => {
            println!("ERROR: Invalid `pair` passed");
            return;
        }
    };

    let now = Instant::now();
    write_seq_pairs((&sequences_0, &sequences_1), (&tags[pair_indices.0], &tags[pair_indices.1]), str_to_seqtype(&folder), rs, p, &folder_path, side);
    let elapsed_time = now.elapsed().as_secs_f32();
    eprintln!("The function took: {elapsed_time} seconds to go through the two sets of pairs\n");
    
}

pub fn create_rowsum_dirs(folder : String, p : usize, rs : (isize, isize, isize, isize), pairing: Option<RowsumPairing>) {
    // This creates the rowsums.quad file as well as the rowsum_x_y_z_w directories, as well as the .pair files
    // For use when directories need to be known/iterated over, but have not been created yet
    // e.g., submitting SLURM jobs with dependencies

    let folder_path = "results/pairs/".to_string()+ &folder + &"/find_" + &p.to_string() + &"/rowsum_" + &(rs.0).to_string() + &"_" + &(rs.1).to_string() + &"_" + &(rs.2).to_string() + &"_" + &(rs.3).to_string();
    println!("{}",folder_path);
    fs::create_dir_all(&folder_path).expect("Error when creating the dir");

    let (_, indices) = sort(&rs); // we sort the rowsum in decreasing order, and we keep track of their original indices
    let tags : Vec<SequenceTag> = indices.iter().map(|i| index_to_tag(*i)).collect(); // we convert the indices to their respective tags

    let path1 : String;
    let path2 : String;

    match pairing {
        Some(RowsumPairing::XY) => {
            path1 = folder_path.clone() + &"/pair_" + &tag_to_string(&tags[0]) + &tag_to_string(&tags[1]) + ".pair";
            path2 = folder_path.clone() + &"/pair_" + &tag_to_string(&tags[2]) + &tag_to_string(&tags[3]) + ".pair";
        },
        Some(RowsumPairing::XZ) => {
            path1 = folder_path.clone() + &"/pair_" + &tag_to_string(&tags[0]) + &tag_to_string(&tags[2]) + ".pair";
            path2 = folder_path.clone() + &"/pair_" + &tag_to_string(&tags[1]) + &tag_to_string(&tags[3]) + ".pair";
        },
        Some(RowsumPairing::XW) => {
            path1 = folder_path.clone() + &"/pair_" + &tag_to_string(&tags[0]) + &tag_to_string(&tags[3]) + ".pair";
            path2 = folder_path.clone() + &"/pair_" + &tag_to_string(&tags[1]) + &tag_to_string(&tags[2]) + ".pair";
        },
        None => {panic!("Missing pairing arg")}
    };

    File::create(path1).expect("Invalid file ?");
    File::create(path2).expect("Invalid file ?");    
}

pub fn write_pairs(p : usize, seqtype : SequenceType, pairing: Option<RowsumPairing>) {
    // This is the starting point of the part of the algorithm that generates the possible sequences

    // all the possible rowsums of p
    let rowsums = generate_rowsums(p);
    for rs in &rowsums {
        eprintln!("{:?}", rs);
    }
    eprintln!("generated {} different rowsums", rowsums.len());

    let folder = seqtype.to_string();
    for rs in rowsums {
        write_pairs_rowsum(&folder, rs, p, pairing.clone());
    }
}

pub fn write_pairs_rowsum(folder : &str, rs : (isize, isize, isize, isize), p : usize, pairing: Option<RowsumPairing>) {
    // This function generates the sequences possible for specific rowsums and stores them
    
    let (rowsums, indices) = sort(&rs); // we sort the rowsum in decreasing order, and we keep track of their original indices
    let tags : Vec<SequenceTag> = indices.iter().map(|i| index_to_tag(*i)).collect(); // we convert the indices to their respective tags
    
    let folder_path = "results/pairs/".to_string()+ &folder + &"/find_" + &p.to_string() + &"/rowsum_" + &(rs.0).to_string() + &"_" + &(rs.1).to_string() + &"_" + &(rs.2).to_string() + &"_" + &(rs.3).to_string();
    println!("{}",folder_path);
    fs::create_dir_all(&folder_path).expect("Error when creating the dir");
    
    let now = Instant::now();
    // We generate all the sequences possible for each rowsums
    let sequences_0 = generate_sequences_with_rowsum(rowsums[0], p);
    let sequences_1 = generate_sequences_with_rowsum(rowsums[1], p);
    let sequences_2 = generate_sequences_with_rowsum(rowsums[2], p);
    let sequences_3 = generate_sequences_with_rowsum(rowsums[3], p);

    write_sequences(&sequences_0, &tags[0], &folder_path);
    write_sequences(&sequences_1, &tags[1], &folder_path);
    write_sequences(&sequences_2, &tags[2], &folder_path);
    write_sequences(&sequences_3, &tags[3], &folder_path);
    
    let elapsed_time = now.elapsed().as_secs_f32();
    eprintln!("The function took: {elapsed_time} seconds to generate sequences with rowsums: {}, {}, {}, {}", rowsums[0], rowsums[1], rowsums[2], rowsums[3]);

    let seqtype = str_to_seqtype(folder);

    let now = Instant::now();

    // Uses sequences to generate .pair files based on chosen pairing (default pairing is XW)
    match pairing {
        Some(RowsumPairing::XY) => {
            write_seq_pairs((&sequences_0, &sequences_1), (&tags[0], &tags[1]), seqtype.clone(), rs, p, &folder_path, EquationSide::LEFT);
            write_seq_pairs((&sequences_2, &sequences_3), (&tags[2], &tags[3]), seqtype.clone(), rs, p, &folder_path, EquationSide::RIGHT);
        },
        Some(RowsumPairing::XZ) => {
            write_seq_pairs((&sequences_0, &sequences_2), (&tags[0], &tags[2]), seqtype.clone(), rs, p, &folder_path, EquationSide::LEFT);
            write_seq_pairs((&sequences_1, &sequences_3), (&tags[1], &tags[3]), seqtype.clone(), rs, p, &folder_path, EquationSide::RIGHT);
        },
        Some(RowsumPairing::XW) | None => {
            write_seq_pairs((&sequences_0, &sequences_3), (&tags[0], &tags[3]), seqtype.clone(), rs, p, &folder_path, EquationSide::LEFT);
            write_seq_pairs((&sequences_1, &sequences_2), (&tags[1], &tags[2]), seqtype.clone(), rs, p, &folder_path, EquationSide::RIGHT);
        }
    };
    
    let elapsed_time = now.elapsed().as_secs_f32();
    eprintln!("The function took: {elapsed_time} seconds to go through the two sets of pairs\n");
}






pub fn join_pairs(p : usize, seqtype : SequenceType) -> Vec<QuadSeq>{
    // This is the starting point of the part of the algorithm that goes through the sorted files and finds valid QTS

    let mut result = vec![];

    let folder = seqtype.to_string();

    let find_i = fs::read_dir("./results/pairs/".to_string() + &folder + &"/find_".to_string() + &p.to_string()).unwrap();
    eprintln!("{}", "./results/pairs/".to_string() + &folder + &"/find_".to_string() + &p.to_string());

    for rowsum_x_y in find_i {
        let directory = rowsum_x_y.unwrap();

        if directory.metadata().unwrap().is_dir() {
            // We read the files in the directory to get back our 4 sets of sequences
            let sequences = get_sequences_from_dir(&directory);
    
            let (pathnames, order) = get_order_from_dir(&directory);
            eprintln!("Folder {} : sequences have order {order:?}", directory.file_name().into_string().unwrap());
    
            result.append(&mut join_pairs_files(&pathnames, seqtype.clone(), &order, &sequences));
        }
    }

    eprintln!("\ncount before equivalences {}", result.len());
    let reduced = reduce_to_equivalence(&result, seqtype);
    eprintln!("count after equivalences {}", reduced.len());

    reduced

}


pub fn get_sequences_from_dir(directory : &DirEntry) -> (Vec<Vec<i8>>,Vec<Vec<i8>>,Vec<Vec<i8>>,Vec<Vec<i8>>) {
    // This function reads the files from a directory and returns the sequences that are in the files ending in .seq

    let mut sequence_x = vec![];
    let mut sequence_y = vec![];
    let mut sequence_z = vec![];
    let mut sequence_w = vec![];

    for file in fs::read_dir(directory.path().display().to_string()).unwrap() {
        let f = file.unwrap();
        let pathname = f.path().display().to_string();
        if pathname.ends_with(".seq") {
            // We loop through files with extension .seq
            // eprintln!("Name: {}", pathname);

            let filename = pathname.split("/").last().expect("No last element ???");
            // eprintln!("Name: {}", filename);
            match filename {
                "seq_X.seq" => {sequence_x = file_to_sequences(&pathname)}
                "seq_Y.seq" => {sequence_y = file_to_sequences(&pathname)}
                "seq_Z.seq" => {sequence_z = file_to_sequences(&pathname)}
                "seq_W.seq" => {sequence_w = file_to_sequences(&pathname)}
                _ => {panic!("Unexpected file ending in .seq : {filename}")}
            }
        }
    }

    (sequence_x, sequence_y, sequence_z, sequence_w)
}


pub fn file_to_sequences(filename : &String) -> Vec<Vec<i8>> {
    // This function reads a file storing a set of sequences and returns them

    let mut seq = vec![];

    for line in read_lines(filename).expect("Error when reading file") {
        if let Ok(s) = line {
            seq.push(string_to_sequence(&s));
        }
    }

    seq
}

pub fn string_to_sequence(s : &String) -> Vec<i8>{
    // This file reads a sequence stored in string form and returns the sequence
    let mut res = vec![];

    for elm in s.chars() {
        match elm {
            '+' => {res.push(1);}
            '-' => {res.push(-1);}
            _ => {panic!("Unexpected entry during string to sequence conversion")}
        }
    }

    res
}



pub fn get_order_from_dir(directory : &DirEntry) -> ((String, String), (SequenceTag, SequenceTag, SequenceTag, SequenceTag)){
    // This function reads the generated files in a folder to determine what comparisons to make

    let mut filenames = vec![];
    let mut pathnames = vec![];

    for file in fs::read_dir(directory.path().display().to_string()).unwrap() {
        let f = file.unwrap();
        let pathname = f.path().display().to_string();
        if pathname.ends_with(".pair.sorted") {
            // We loop through files with extension .pair
            // eprintln!("Name: {}", pathname);
            pathnames.push(pathname.clone());

            let filename = pathname.split("/").last().expect("No last element ???");
            //eprintln!("Name: {}", filename);
            filenames.push(filename.to_string());
        }
    }

    assert!(filenames.len() == 2, "too many or not enough .pair.sorted files !");

    let (tag1, tag2) = get_tag_from_filename(&filenames[0]);
    let (tag3, tag4) = get_tag_from_filename(&filenames[1]);


    ((pathnames[0].to_string(), pathnames[1].to_string()),(tag1, tag2, tag3, tag4))
}

pub fn get_tag_from_filename(filename : &str) -> (SequenceTag, SequenceTag) {
    // This reads the name of a file and returns what are the corresponding tags

    let tag1 = match filename.chars().nth(5).expect("File name not long enough") {
        'X' => {SequenceTag::X}
        'Y' => {SequenceTag::Y}
        'Z' => {SequenceTag::Z}
        'W' => {SequenceTag::W}
        _ => {panic!("Unexpected character : {}",filename.chars().nth(5).unwrap())}
    };
    let tag2 = match filename.chars().nth(6).expect("File name not long enough") {
        'X' => {SequenceTag::X}
        'Y' => {SequenceTag::Y}
        'Z' => {SequenceTag::Z}
        'W' => {SequenceTag::W}
        _ => {panic!("Unexpected character : {}",filename.chars().nth(6).unwrap())}
    };

    (tag1, tag2)
}






pub fn join_pairs_files(filenames : &(String, String), seqtype : SequenceType, order : &(SequenceTag, SequenceTag, SequenceTag, SequenceTag), sequences : &(Vec<Vec<i8>>, Vec<Vec<i8>>, Vec<Vec<i8>>, Vec<Vec<i8>>)) -> Vec<QuadSeq> {
    // This function reads two sorted files of sequences and uses the order to determine what comparisons should be made, and returns the valid QTS

    let mut result = vec![];

    let (file12, file34) = filenames;

    let mut lines12 = read_lines(file12).expect("Invalid vile somehow ?");
    let mut lines34 = read_lines(file34).expect("Invalid vile somehow ?");

    let mut line12 = lines12.next();
    let mut line34 = lines34.next();

    while line12.is_some() && line34.is_some() {
        // We loop until there's no more lines to read

        let (mut seq12, mut indices12) = get_line_from(&line12);
        let (mut seq34, mut indices34) = get_line_from(&line34);
        
        let current_seq = seq12.clone();
        if seq12 == seq34 {

            // Store every sequence with the same values of auto/cross correlation
            let mut possible_matching_12 = vec![];
            while line12.is_some() && current_seq == seq12 {
                possible_matching_12.push(indices12);
                (seq12, indices12) = get_line_from(&line12);
                line12 = lines12.next();
            }

            // Store every sequence here as well
            let mut possible_matching_34 = vec![];
            while line34.is_some() && current_seq == seq34 {
                possible_matching_34.push(indices34);
                (seq34, indices34) = get_line_from(&line34);
                line34 = lines34.next();
            }

            // Loop through the possible matches
            for ((i1, i2),(i3, i4)) in iproduct!(possible_matching_12, possible_matching_34) {
                let indices = (i1, i2, i3, i4);
                // test if the sequence is of type seqtype, add them to the result files if it is
                let sequences = get_sequences(sequences, order, &indices);
                let mut quad_seq = QuadSeq::new(sequences.0.len());
                quad_seq.set_all_values(sequences);
                
                let condition: Box<dyn Fn(&QuadSeq) -> bool> = match seqtype {
                    SequenceType::QuaternionType => Box::new(|quad| quad.to_qs().is_perfect()),
                    SequenceType::WilliamsonType => Box::new(|quad| quad.verify_wts()),
                    SequenceType::Williamson => Box::new(|quad| quad.verify_ws()),
                    _ => Box::new(|_| false)
                };
                if condition(&quad_seq) {
                    result.push(quad_seq);
                }
            }
        }
        else if seq12 < seq34 {
            line12 = lines12.next();
        }
        else {
            line34 = lines34.next();
        }

    }

    result
}


pub fn get_line_from(seq : &Option<Result<String, Error>>) -> (String, (usize, usize)) {
    // This function parse a line from the sorted files

    assert!(seq.is_some());

    let s = seq.as_ref().unwrap().as_ref().expect("Error when reading file").clone();
    let mut s_parts = s.split("_:_");

    let line = s_parts.next().expect("Expected something before ':' !").to_string();
    let mut indices_parts = s_parts.next().expect("Expected something after ':' !").split("_");
    let indices = (indices_parts.next().expect("Expected something !").parse().expect("Expected a number !"), indices_parts.next().expect("Expected something !").parse().expect("Expected a number !"));

    (line, indices)
}


pub fn get_sequences<'a>(sequences : &'a (Vec<Vec<i8>>, Vec<Vec<i8>>, Vec<Vec<i8>>, Vec<Vec<i8>>), order : &'a (SequenceTag, SequenceTag, SequenceTag, SequenceTag), indices : &'a (usize, usize, usize, usize)) -> (&'a Vec<i8>, &'a Vec<i8>, &'a Vec<i8>, &'a Vec<i8>) {
    // This function returns the sequences corresponding to the indices in a specific order

    let seqx = get_sequence_aux(sequences, order, indices, SequenceTag::X);
    let seqy = get_sequence_aux(sequences, order, indices, SequenceTag::Y);
    let seqz = get_sequence_aux(sequences, order, indices, SequenceTag::Z);
    let seqw = get_sequence_aux(sequences, order, indices, SequenceTag::W);

    (seqx, seqy, seqz, seqw)
}


fn get_sequence_aux<'a>(sequences : &'a (Vec<Vec<i8>>, Vec<Vec<i8>>, Vec<Vec<i8>>, Vec<Vec<i8>>), order : &'a (SequenceTag, SequenceTag, SequenceTag, SequenceTag), indices : &(usize, usize, usize, usize), tag : SequenceTag)  -> &'a Vec<i8>{

    let (tag1, tag2, tag3 ,tag4) = order;

    let index = match tag {
        _ if tag == *tag1 => {indices.0}
        _ if tag == *tag2 => {indices.1}
        _ if tag == *tag3 => {indices.2}
        _ if tag == *tag4 => {indices.3}
        _ => {panic!("Problem with order !")}
    };

    match tag {
        SequenceTag::X => {&sequences.0[index]},
        SequenceTag::Y => {&sequences.1[index]},
        SequenceTag::Z => {&sequences.2[index]},
        SequenceTag::W => {&sequences.3[index]},
    }
}
