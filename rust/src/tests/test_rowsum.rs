


#[cfg(test)]
mod tests {
    use crate::sequences::rowsum::*;

    #[test]
    fn test_rowsum(){

        let seq = vec![1;18];
        assert_eq!(rowsum(seq), 18);

        let seq = vec![1,-1,-1,-1,1,-1,1,1,-1,-1,1,1];
        assert_eq!(rowsum(seq), 0);

    }

    fn choose(n : usize, k : usize) -> usize{
        if k == 0 {return 1}
        return (n * choose(n - 1, k - 1)) / k
    }


    #[test]
    fn test_generate_sequence(){
        
        let tests = vec![(4,3), (6,4), (2,1), (4,0), (7,3)];
        
        for (n,k) in tests{
            let seq : Vec<i8> = vec![-1;n];
            let res = gen_seq_rec(&seq, k, 0);
            assert_eq!(res.len(), choose(n,k));
        }
    }


    #[test]
    fn test_rowsum_generate(){
        let tests = vec![(4,2), (6,4), (7,3), (4,0)];
        
        for (n,r) in tests{
            let res = generate_sequences_with_rowsum(r, n);
            assert_eq!(res.len(), choose(n,(n + r)/2));
        }

        
        let res = generate_sequences_with_rowsum(4, 3);
        assert_eq!(res.len(), 0);
    }


    #[test]
    fn test_four_squares(){

        let list = sum_of_four_squares(4*19);

        assert_eq!(list.len(), 5);
        assert_eq!(list[0], (0,2,6,6));
        assert_eq!(list[1], (1,1,5,7));
        assert_eq!(list[2], (1,5,5,5));
        assert_eq!(list[3], (2,2,2,8));
        assert_eq!(list[4], (3,3,3,7));

    }

    #[test]
    fn test_four_squares2(){

        let list = sum_of_four_squares(4*17);

        for elm in list {
            println!("{:?}", elm);
        }

        panic!();
    }


    #[test]
    fn test_gen_quads(){

        let quads = generate_other_quadruplets(&(1,1,5,7));
        assert_eq!(quads.len(), 3);
        assert_eq!(quads[0], (1,1,5,7));
        assert_eq!(quads[1], (1,5,1,7));
        assert_eq!(quads[2], (1,5,7,1));
        
        let quads = generate_other_quadruplets(&(1,3,3,7));
        assert_eq!(quads.len(), 3);
        assert_eq!(quads[0], (1,3,3,7));
        assert_eq!(quads[1], (1,3,7,3));
        assert_eq!(quads[2], (1,7,3,3));

        let quads = generate_other_quadruplets(&(3,3,3,8));
        assert_eq!(quads.len(), 1);
        assert_eq!(quads[0], (3,3,3,8));

        let quads = generate_other_quadruplets(&(3,3,5,5));
        assert_eq!(quads.len(), 3);
        assert_eq!(quads[0], (3,3,5,5));
        assert_eq!(quads[1], (3,5,3,5));
        assert_eq!(quads[2], (3,5,5,3));
    }



    #[test]
    fn test_possible_rowsums(){

        let rowsums = generate_rowsums(17);
        
        assert_eq!(rowsums.len(), 6);
        assert_eq!(rowsums[0], (1,3,3,7));
        assert_eq!(rowsums[1], (1,3,7,3));
        assert_eq!(rowsums[2], (1,7,3,3));
        assert_eq!(rowsums[3], (3,3,5,5));
        assert_eq!(rowsums[4], (3,5,3,5));
        assert_eq!(rowsums[5], (3,5,5,3));
    }



    #[test]
    fn test_aux(){
        
        let res = generate_sequences_with_rowsum(7, 11);

        for elm in res {
            println!("{:?}", elm);
        }
        
        panic!();
    }
}