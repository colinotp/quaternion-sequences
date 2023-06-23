

#[cfg(test)]
mod tests {
    use crate::sequences::equations::{generate_equation_from, equations_crosscorrelation, OpType, AddType};



    #[test]
    fn test_single_gen() {

        let seq = vec![1,1,1,-1,-1,1,-1,1,1,0,0,1,0,1];
        generate_equation_from(&seq, 5);
        // output:
        // 1 x1 +1 x2 +1 x3 -1 x4 -1 x5 +1 x6 -1 x7 +1 x8 +1 x9 +1 x12 +1 x14 = 5;

        let seq = vec![0,0,0,0,0,0];
        generate_equation_from(&seq, 5);
        // nothing should be printed

        let seq = vec![1,-1,-1,-1,-2,-3,4,2,5,9,10];
        generate_equation_from(&seq, 12);
        // output:
        // 1 x1 -1 x2 -1 x3 -1 x4 -2 x5 -3 x6 +4 x7 +2 x8 +5 x9 +9 x10 +10 x11 = 12;
        
        //panic!();
    }

    #[test]
    fn test_single_cross() {

        let seq1 = vec![1,1,-1,-1];
        let seq2 = vec![1,-1,1,1];
        equations_crosscorrelation(&seq1, OpType::LeftMinus, &seq2, OpType::RightMinus, AddType::Plus);
        // expects output:
        // -4 x1 +4 x2 +4 x3 -4 x4 -4 x5 +4 x7 = 0;
        // 4 x1 -4 x2 -4 x3 +4 x4 +4 x5 -4 x7 = 0;

        println!("=============");

        let seq1 = vec![1,1,-1];
        let seq2 = vec![1,-1,1];
        equations_crosscorrelation(&seq1, OpType::RightPlus, &seq2, OpType::LeftMinus, AddType::Minus);
        // expects result:
        // 4 x1 +4 x2 -4 x3 = 2;
        // 4 x3 -4 x4 +4 x6 = 2;
        // 4 x3 +4 x4 -4 x6 = 2;

        //panic!();
    }

}