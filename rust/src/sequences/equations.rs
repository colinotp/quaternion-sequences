use super::{williamson::{SequenceTag, tag_to_string, periodic_autocorrelation, cross_correlation}, symmetries::SequenceType, rowsum::Quad};



pub enum OpType { // Type to know what operator to use when computing the values for the equations
    LeftPlus, LeftMinus, RightPlus, RightMinus
}

pub enum AddType {
    Plus, Minus
}


pub fn generate_equations(seq1 : &Vec<i8>, tag1 : &SequenceTag, seq2 : &Vec<i8>, tag2 : &SequenceTag, seqtype : &SequenceType, rowsum : &Quad) -> String{
    // generates all the equations given by two sequences and their positions, so that they return a specific type of sequence.

    match seqtype {
        SequenceType::QuaternionType => {
            equations_williamson_type(seq1, tag1, seq2, tag2, rowsum)
        }
        _ => {panic!()/* TODO */}
    }
}


fn equations_williamson_type(seq1 : &Vec<i8>, tag1 : &SequenceTag, seq2 : &Vec<i8>, tag2 : &SequenceTag, rowsum : &Quad) -> String {
    // generates the equations given by two sequences and their positions so that that they return Williamson-type sequences
    let mut result = "".to_string();
    
    match (tag1, tag2) { // ! In case of multiple sequences, make sure the same sequences are on the same side !!!
        (SequenceTag::W, SequenceTag::X) => {// x1 ... is D, ... xn is C
            result += &additional_comment(&SequenceTag::Z, &SequenceTag::Y);
            result += &equations_crosscorrelation(&seq1, OpType::RightMinus, &seq2, OpType::LeftMinus, AddType::Plus);
            result += &equations_crosscorrelation(&seq2, OpType::RightMinus, &seq1, OpType::RightMinus, AddType::Plus);
            result += &equations_nonlinear_crosscorrelation(&seq1, &seq2, OpType::LeftMinus, OpType::LeftMinus, AddType::Plus);
            result += &equations_rowsum(seq1.len(), rowsum.3, rowsum.2);
        }
        (SequenceTag::W, SequenceTag::Y) => {// x1 ... is D, ... xn is B
            result += &additional_comment(&SequenceTag::Z, &SequenceTag::X);
            result += &equations_crosscorrelation(&seq1, OpType::RightMinus, &seq2, OpType::RightMinus, AddType::Plus);
            result += &equations_crosscorrelation(&seq2, OpType::RightMinus, &seq1, OpType::LeftMinus, AddType::Plus);
            result += &equations_nonlinear_crosscorrelation(&seq1, &seq2, OpType::RightMinus, OpType::LeftMinus, AddType::Plus);
            result += &equations_rowsum(seq1.len(), rowsum.3, rowsum.1);
        }
        (SequenceTag::W, SequenceTag::Z) => {// x1 ... is B, ... xn is C
            result += &additional_comment(&SequenceTag::X, &SequenceTag::Y);
            result += &equations_crosscorrelation(&seq2, OpType::LeftMinus, &seq1, OpType::RightMinus, AddType::Plus);
            result += &equations_crosscorrelation(&seq1, OpType::LeftMinus, &seq2, OpType::LeftMinus, AddType::Plus);
            result += &equations_nonlinear_crosscorrelation(&seq1, &seq2, OpType::RightMinus, OpType::LeftMinus, AddType::Plus);
            result += &equations_rowsum(seq1.len(), rowsum.1, rowsum.2);
        }
        (SequenceTag::X, SequenceTag::Y) => {// x1 ... is D, ... xn is A
            result += &additional_comment(&SequenceTag::Z, &SequenceTag::W);
            result += &equations_crosscorrelation(&seq1, OpType::RightMinus, &seq2, OpType::LeftMinus, AddType::Plus);
            result += &equations_crosscorrelation(&seq2, OpType::RightMinus, &seq1, OpType::RightMinus, AddType::Plus);
            result += &equations_nonlinear_crosscorrelation(&seq1, &seq2, OpType::LeftMinus, OpType::LeftMinus, AddType::Plus);
            result += &equations_rowsum(seq1.len(), rowsum.3, rowsum.0);
        }
        (SequenceTag::X, SequenceTag::Z) => {// x1 ... is A, ... xn is C
            result += &additional_comment(&SequenceTag::W, &SequenceTag::Y);
            result += &equations_crosscorrelation(&seq2, OpType::LeftMinus, &seq1, OpType::LeftMinus, AddType::Plus);
            result += &equations_crosscorrelation(&seq1, OpType::RightMinus, &seq2, OpType::LeftMinus, AddType::Plus);
            result += &equations_nonlinear_crosscorrelation(&seq1, &seq2, OpType::RightMinus, OpType::RightMinus, AddType::Plus);
            result += &equations_rowsum(seq1.len(), rowsum.0, rowsum.2);
        }
        (SequenceTag::Y, SequenceTag::Z) => {// x1 ... is A, ... xn is B
            result += &additional_comment(&SequenceTag::W, &SequenceTag::X);
            result += &equations_crosscorrelation(&seq2, OpType::LeftMinus, &seq1, OpType::RightMinus, AddType::Plus);
            result += &equations_crosscorrelation(&seq1, OpType::LeftMinus, &seq2, OpType::LeftMinus, AddType::Plus);
            result += &equations_nonlinear_crosscorrelation(&seq1, &seq2, OpType::RightMinus, OpType::LeftMinus, AddType::Plus);
            result += &equations_rowsum(seq1.len(), rowsum.0, rowsum.1);
        }
        _ => {panic!("incorrect tags entered !")}
    }

    result += &equations_autocorrelation(&seq1, &seq2);

    result
}



 // TODO refactor these functions
fn op_left_plus(seq : &Vec<i8>, k : usize, t : usize, n : usize) -> i8{
    seq[((k+n) as isize - t as isize) as usize % n] + seq[(k+t) % n]
}

fn op_right_plus(seq : &Vec<i8>, k : usize, t : usize, n : usize) -> i8{
    seq[(k+t) % n] + seq[((k+n) as isize - t as isize) as usize % n]
}

fn op_left_minus(seq : &Vec<i8>, k : usize, t : usize, n : usize) -> i8{
    seq[((k+n) as isize - t as isize) as usize % n] - seq[(k+t) % n]
}
fn op_right_minus(seq : &Vec<i8>, k : usize, t : usize, n : usize) -> i8{
    seq[(k+t) % n] - seq[((k+n) as isize - t as isize) as usize % n]
}



pub fn equations_crosscorrelation(seq1 : &Vec<i8>, op_type1 : OpType, seq2 : &Vec<i8>, op_type2 : OpType, add_type : AddType) -> String {
    // generates the equations given by applying the crosscorrelation properties on each sequence with their specific operator: + or -

    let op1 = match op_type1 { // op1 is a function corresponding to either a + or a - operator
        OpType::LeftPlus => { op_left_plus}
        OpType::LeftMinus => { op_left_minus}
        OpType::RightPlus => { op_right_plus}
        OpType::RightMinus => { op_right_minus}
    };

    let op2 = match (op_type2, add_type) { // op2 is a function corresponding to either a + or a - operator
        (OpType::LeftPlus, AddType::Plus) => { |s,k,i,n| op_left_plus(s,k,i,n)}
        (OpType::LeftMinus, AddType::Plus) => { |s,k,i,n| op_left_minus(s,k,i,n)}
        (OpType::RightPlus, AddType::Plus) => { |s,k,i,n| op_right_plus(s,k,i,n)}
        (OpType::RightMinus, AddType::Plus) => { |s,k,i,n| op_right_minus(s,k,i,n)}
        (OpType::LeftPlus, AddType::Minus) => { |s,k,i,n|  -op_left_plus(s,k,i,n)}
        (OpType::LeftMinus, AddType::Minus) => { |s,k,i,n| -op_left_minus(s,k,i,n)}
        (OpType::RightPlus, AddType::Minus) => { |s,k,i,n| -op_right_plus(s,k,i,n)}
        (OpType::RightMinus, AddType::Minus) => { |s,k,i,n| -op_right_minus(s,k,i,n)}
    };

    let n = seq1.len();
    let mut values = vec![0 ; 2*n];
    let mut result = "".to_string();

    for t in 0..n { 
        // computes the n equations given by the crosscorrelation property
        
        for k in 0..n {
            values[k] = op1(seq1,k,t,n); 
        }
        
        for k in 0..n {
            values[k+n] = op2(seq2,k,t,n); 
        }
        
        let mut rightside_value = 0;
        for k in 0..2*n {
            // we convert the variables from {+-1} to {0,1} by applying x -> 2x-1
            // Then we divide all the values by 2
            rightside_value += (values[k]/2) as isize;
        }

        result += &generate_equation_from(&values, rightside_value);
    }

    result
}


pub fn generate_equation_from(coefficients : &Vec<i8>, rightside_value : isize) -> String{
    // generates the equation given by the specific values
    // The value on the right side of the equation is computed automatically unless specified otherwise


    let mut result = "".to_string();
    let mut variable_counter = 1; // keeps track of the number of variable
    let mut all_zero_coeff = true;


    for &elm in coefficients.iter() {
        if elm > 0 {
            if !all_zero_coeff {
                result += "+";
            }
            result += &(elm.to_string() + " x" + &variable_counter.to_string() + " ");
            all_zero_coeff = false;
        }
        if elm < 0 {
            result += &(elm.to_string() + " x" + &variable_counter.to_string() + " ");
            all_zero_coeff = false;
        }
        variable_counter += 1;
    }

    if all_zero_coeff {
        // the values were all 0, so nothing can be learned from it
        return "".to_string();
    }

    result += &("= ".to_owned() + &rightside_value.to_string() + &";"  + &"\n");


    result
}




fn equations_rowsum(n : usize, rs1 : isize, rs2 : isize) -> String {
    let mut result = "".to_string();
    let mut coeff1 = vec![0;2*n];
    let mut coeff2 = vec![0;2*n];
    for i in 0..n {
        coeff1[i] = 2;
        coeff2[2*n-1-i] = 2;
    }

    result += &(generate_equation_from(&coeff1, rs1 + n as isize));
    result += &(generate_equation_from(&coeff2, rs2 + n as isize));

    result
}







fn equations_autocorrelation(seq1 : &Vec<i8>, seq2 : &Vec<i8>) -> String {
    let mut result = "".to_string();
    let n = seq1.len();

    for offset in 1..n {
        let rightside_value = periodic_autocorrelation(seq1, offset) + periodic_autocorrelation(seq2, offset);   
        result += &generate_equation_autocorrelation(n, offset, -rightside_value);
    }

    result

}

fn generate_equation_autocorrelation(size : usize, offset : usize, mut rightside_value : isize) -> String {
    let mut result = "".to_string();
    let mut first_time = true;
    for index in 1..=size {

        if !first_time {
            result += "+";
        }

        let index_offset = (index + offset) % size;
        result += &("4 x".to_string() + &(index).to_string() + &" x" + &index_offset.to_string() + " ");
        result += &("-2 x".to_string() + &(index).to_string() + &" -2 x" + &index_offset.to_string() + " ");

        result += &("+4 x".to_string() + &(index + size).to_string() + &" x" + &(index_offset + size).to_string() + " ");
        result += &("-2 x".to_string() + &(index + size).to_string() + &" -2 x" + &(index_offset + size).to_string() + " ");

        rightside_value += -2;

        first_time = false;
    }
    result += &("= ".to_string() + &rightside_value.to_string() + &";\n");

    result
}








fn op_minus(seqa : &Vec<i8>, seqb : &Vec<i8>, offset: usize) -> isize{
    cross_correlation(seqa, seqb, offset) - cross_correlation(seqa, seqb, offset)
}

fn op_plus(seqa : &Vec<i8>, seqb : &Vec<i8>, offset: usize) -> isize{
    cross_correlation(seqa, seqb, offset) + cross_correlation(seqa, seqb, offset)
}



fn equations_nonlinear_crosscorrelation(seq1 : &Vec<i8>, seq2 : &Vec<i8>, op_type1 : OpType, op_type2 : OpType, add_type : AddType) -> String {
    let mut result = "".to_string();
    let n = seq1.len();


    let op = match op_type1 {
        OpType::LeftPlus | OpType::RightPlus => {|seq1, seq2, offset| op_plus(seq1, seq2, offset)}
        OpType::LeftMinus => {|seq1, seq2, offset| op_minus(seq1, seq2, offset)}
        OpType::RightMinus => {|seq1, seq2, offset| op_minus(seq2, seq1, offset)}
    };
    let op2 = match add_type {
        AddType::Minus => { |x : isize| -x}
        AddType::Plus => { |x : isize| x}
    };



    for offset in 1..n {
        let rightside_value = op2(-op(seq1,seq2, offset));   
        result += &generate_equation_non_linear_crosscorrelation(n, offset, &op_type2, -rightside_value);
    }

    result

}




fn generator_left_minus(index : usize, offset : usize, size : usize) -> String {
    let index_offset = (index + offset) % size;

    let mut result = "".to_string();
    result += &("4 x".to_string() + &(index).to_string() + &" x" + &(index_offset + size).to_string() + " ");
    result += &("-2 x".to_string() + &(index).to_string() + &" -2 x" + &(index_offset + size).to_string() + " ");

    result += &("-4 x".to_string() + &(index + size).to_string() + &" x" + &index_offset.to_string() + " ");
    result += &("+2 x".to_string() + &(index + size).to_string() + &" +2 x" + &index_offset.to_string() + " ");

    result
}

fn generator_right_minus(index : usize, offset : usize, size : usize) -> String {
    let index_offset = (index + offset) % size;

    let mut result = "".to_string();
    result += &("4 x".to_string() + &(index + size).to_string() + &" x" + &index_offset.to_string() + " ");
    result += &("-2 x".to_string() + &(index + size).to_string() + &" -2 x" + &index_offset.to_string() + " ");

    result += &("-4 x".to_string() + &(index).to_string() + &" x" + &(index_offset + size).to_string() + " ");
    result += &("+2 x".to_string() + &(index).to_string() + &" +2 x" + &(index_offset + size).to_string() + " ");

    result
}

fn generator_plus(index : usize, offset : usize, size : usize) -> String {
    let index_offset = (index + offset) % size;

    let mut result = "".to_string();
    result += &("4 x".to_string() + &(index).to_string() + &" x" + &(index_offset + size).to_string() + " ");
    result += &("-2 x".to_string() + &(index).to_string() + &" -2 x" + &(index_offset + size).to_string() + " ");

    result += &("+4 x".to_string() + &(index + size).to_string() + &" x" + &index_offset.to_string() + " ");
    result += &("-2 x".to_string() + &(index + size).to_string() + &" -2 x" + &index_offset.to_string() + " ");

    result
}





fn generate_equation_non_linear_crosscorrelation(size : usize, offset : usize, op_type : &OpType, rightside_value : isize) -> String {
    let mut result = "".to_string();
    let mut first_time = true;

    let line_generator = match *op_type {
        OpType::LeftMinus => {generator_left_minus}
        OpType::RightMinus => {generator_right_minus}
        OpType::LeftPlus => {generator_plus}
        OpType::RightPlus => {generator_plus}
    };


    for index in 1..=size {

        if !first_time {
            result += "+";
        }

        result += &line_generator(index, offset, size);

        first_time = false;
    }
    result += &("= ".to_string() + &rightside_value.to_string() + &";\n");

    result
}




fn additional_comment(tag1 : &SequenceTag, tag2 : &SequenceTag) -> String {
    "\n* x1 ... xn is ".to_string() + &tag_to_string(tag1) + &", xn+1 ... x2n is " + &tag_to_string(tag2) + &"\n"
}