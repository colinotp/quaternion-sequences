use cgmath::Quaternion;

use super::{sequence::{QS, quaternion_to_string}, williamson::QuadSeq, symmetries::SequenceType};



pub struct QHM {
    size : usize,
    matrix : Vec<Vec<Quaternion<f32>>>
}


impl QHM {

    pub fn new(size : usize) -> QHM {
        let matrix = vec![vec![Quaternion::new(1.,0.,0.,0.); size]; size];
        QHM {
            size,
            matrix
        }
    }


    pub fn get(&self, i : usize, j : usize) -> Quaternion<f32> {
        self.matrix[i][j]
    }

    // Returns reference to specified row (0 indexed)
    pub fn row(&self, row : usize) -> Vec<Quaternion<f32>> {
        self.matrix[row].clone()
    }

    pub fn col(&self, col : usize) -> Vec<Quaternion<f32>> {
        let mut col_vec = Vec::new();
        for i in 0..self.size() {
            col_vec.push(self.matrix[i][col]);
        }

        col_vec
    }


    


    pub fn from_pqs(pqs : QS) -> QHM {

        let size = pqs.size();
        let mut matrix = vec![];

        for row in 0..size {
            let mut row_vec = vec![];
            for col in 0..size {
                row_vec.push(pqs.values()[(col + size - row) % size].clone())
            }

            matrix.push(row_vec);
        }

        QHM {size, matrix}
    }

    pub fn size(&self) -> usize{
        self.size
    }




    pub fn dephase(&mut self) {
        // dephase the columns
        for col in 0..self.size {
            for inv_row in 0..self.size {
                let row = self.size - 1 - inv_row;
                self.matrix[row][col] = self.matrix[row][col]*self.matrix[0][col].conjugate();
            }
        }
        // dephase the rows
        for row in 0..self.size {
            for inv_col in 0..self.size {
                let col = self.size - 1 - inv_col;
                self.matrix[row][col] = self.matrix[row][0].conjugate()*self.matrix[row][col];
            }
        }
    }


    pub fn dephased(&self) -> QHM{
        let mut new_mat = vec![vec![Quaternion::<f32>::new(0.,0.,0.,0.);self.size]; self.size];

        // dephase the columns
        for col in 0..self.size {
            for inv_row in 0..self.size {
                let row = self.size - 1 - inv_row;
                new_mat[row][col] = self.matrix[row][col]*self.matrix[0][col].conjugate();
            }
        }
        // dephase the rows
        for row in 0..self.size {
            for inv_col in 0..self.size {
                let col = self.size - 1 - inv_col;
                new_mat[row][col] = new_mat[row][0].conjugate()*new_mat[row][col];
            }
        }

        QHM{size : self.size, matrix : new_mat}
    }


    pub fn contains_non_commuting_elements(&self) -> bool {
        let mut unique_elements : Vec<Quaternion<f32>> = Vec::new();
        
        // Iterate over elements in matrix
        for i in 1..self.size() {
            for j in 1..self.size() {
                let new = self.get(i,j);
                // Check if we have encountered this element before
                if !unique_elements.contains(&new) {
                    // If this is a new element, check if it commutes with all of the previous elements we have found
                    for elm in &unique_elements {
                        let left = new * elm;
                        let right = elm * new;

                        if left != right {
                            return true;
                        }
                    }
                    
                    // It has been shown to commute with all other elements so far, add it to unique elements list
                    unique_elements.push(new);
                }
            }
        }
        
        // All elements found commute
        false
    }

    // Verifies QHM property
    pub fn verify(&self) -> bool {
        let n = self.size();
        let f32_tolerance : f32 = f32::EPSILON.sqrt();

        // Take complex inner product of each row with each other row. Should get n when taking product with a row and itself, 0 otherwise.
        for row1 in 0..n {
            for row2 in 0..n {
                // Get two rows
                let row1_vec = self.row(row1);
                let row2_vec = self.row(row2);

                // Perform inner product
                let mut result = Quaternion::<f32>::new(0.0,0.0,0.0,0.0);

                for i in 0..n {
                    result = result + (row1_vec[i] * (row2_vec[i].conjugate()));
                }

                // Check if the real component is what it should be
                let real_match = if row1 == row2 {
                    (result.s - n as f32).abs() < f32_tolerance
                } else {
                    result.s.abs() < f32_tolerance
                };

                // Check that the imaginary components are all 0
                let imaginary_match = {
                    result.v.x.abs() < f32_tolerance &&
                    result.v.y.abs() < f32_tolerance &&
                    result.v.z.abs() < f32_tolerance
                };

                // If either are incorrect, then matrix is not QHM
                if !(real_match && imaginary_match) {
                    return false;
                }
            }
        }

        // If we have not yet found an entry that does not match with nI_n, then it is a QHM
        true
    }


    pub fn to_string(&self) -> String {

        let mut result = "".to_string();

        for row in 0..self.size {
            result += "| ";
            for col in 0..self.size {
                result += &(quaternion_to_string(&self.matrix[row][col]) + &" ");
            }
            result += "|\n";


        }

        result
    }

}


#[derive(PartialEq, Eq, Hash, Clone, Debug)]
pub struct HM {
    size : usize,
    matrix : Vec<Vec<i8>>
}

pub enum OpMat {
    NONE, MINUS, TRANSPOSE, MINUSTRANSPOSE
}


impl HM {
    
    pub fn new(size : usize) -> HM {
        let matrix = vec![vec![1; size]; size];
        HM {
            size,
            matrix
        }
    }


    pub fn get(&self, i : usize, j : usize) -> i8 {
        self.matrix[i][j]
    }

    pub fn set_value(&mut self, row : usize, col : usize, value : i8) {
        self.matrix[row][col] = value;
    }


    pub fn from_sequence(seq : &Vec<i8>) -> HM {

        let size = seq.len();
        let mut hm = HM::new(size);

        for (index, elm) in seq.iter().enumerate() {
            for i in 0..size {
                hm.set_value(i, (i+index) % size, *elm);
            }
        }

        hm
    }

    // Forms HM from QuadSeq 
    pub fn from_williamson(will : &QuadSeq, seqtype : SequenceType) -> HM {

        let size = will.size();
        let mut hm = HM::new(4*size);

        let matw = HM::from_sequence(&will.sequence(super::williamson::SequenceTag::W));
        let matx = HM::from_sequence(&will.sequence(super::williamson::SequenceTag::X));
        let maty = HM::from_sequence(&will.sequence(super::williamson::SequenceTag::Y));
        let matz = HM::from_sequence(&will.sequence(super::williamson::SequenceTag::Z));

        match seqtype {
            SequenceType::QuaternionType => {
                hm.copy_block_to(&matw, 0, 0, &OpMat::NONE);
                hm.copy_block_to(&matx, size, 0, &OpMat::NONE);
                hm.copy_block_to(&maty, 2*size, 0, &OpMat::NONE);
                hm.copy_block_to(&matz, 3*size, 0, &OpMat::NONE);
                hm.copy_block_to(&matx, 0, size, &OpMat::NONE);
                hm.copy_block_to(&matw, size, size, &OpMat::MINUS);
                hm.copy_block_to(&matz, 2*size, size, &OpMat::NONE);
                hm.copy_block_to(&maty, 3*size, size, &OpMat::MINUS);
                hm.copy_block_to(&maty, 0, 2*size, &OpMat::NONE);
                hm.copy_block_to(&matz, size, 2*size, &OpMat::MINUS);
                hm.copy_block_to(&matw, 2*size, 2*size, &OpMat::MINUS);
                hm.copy_block_to(&matx, 3*size, 2*size, &OpMat::NONE);
                hm.copy_block_to(&matz, 0, 3*size, &OpMat::NONE);
                hm.copy_block_to(&maty, size, 3*size, &OpMat::NONE);
                hm.copy_block_to(&matx, 2*size, 3*size, &OpMat::MINUS);
                hm.copy_block_to(&matw, 3*size, 3*size, &OpMat::MINUS);
            }
            _ => {panic!("Not Implemented yet !!!")}
        }

        hm
    }


    // Inserts a block matrix into self with upper-leftmost index at row_offset, col_offset
    // OpMat for applying operations to block when inserting (e.g, insert -X transpose)
    pub fn copy_block_to(&mut self, block : &HM, row_offset : usize, col_offset : usize, opmat : &OpMat) { 
        for row in 0..block.size {
            for col in 0..block.size {
                self.matrix[row + row_offset][col + col_offset] = block.get_with_op_mat(row, col, opmat)
            }
        }
    }

    pub fn get_qts(&self) -> QuadSeq {
        let len = self.size() / 4;
        let mut qts = QuadSeq::new(len);

        let seq_w : Vec<i8> = (0..len).map(|i| self.get(0, i)).collect();
        let seq_x : Vec<i8> = (0..len).map(|i| self.get(len, i)).collect();
        let seq_y : Vec<i8> = (0..len).map(|i| self.get(2*len, i)).collect();
        let seq_z : Vec<i8> = (0..len).map(|i| self.get(3*len, i)).collect();

        qts.set_all_values((&seq_w, &seq_x, &seq_y, &seq_z));

        qts
    }

    fn get_with_op_mat(&self, row : usize, col : usize, opmat : &OpMat) -> i8 {
        match opmat {
            OpMat::NONE => {self.matrix[row][col]}
            OpMat::MINUS => {-self.matrix[row][col]}
            OpMat::TRANSPOSE => {self.matrix[col][row]}
            OpMat::MINUSTRANSPOSE => {-self.matrix[col][row]}
        }
    }



    pub fn size(&self) -> usize{
        self.size
    }




    pub fn dephase(&mut self) {
        // dephase the rows
        for inv_row in 0..self.size {
            for col in 0..self.size {
                let row = self.size - 1 - inv_row;
                self.matrix[row][col] = self.matrix[0][col]*self.matrix[row][col];
            }
        }
        // dephase the columns
        for inv_col in 0..self.size {
            for row in 0..self.size {
                let col = self.size - 1 - inv_col;
                self.matrix[row][col] = self.matrix[row][0]*self.matrix[row][col];
            }
        }
    }


    pub fn dephased(&self) -> HM{
        let mut new_mat = vec![vec![1;self.size]; self.size];

        // dephase the rows
        for inv_row in 0..self.size {
            for col in 0..self.size {
                let row = self.size - 1 - inv_row;
                new_mat[row][col] = self.matrix[0][col]*self.matrix[row][col];
            }
        }
        // dephase the columns
        for inv_col in 0..self.size {
            for row in 0..self.size {
                let col = self.size - 1 - inv_col;
                new_mat[row][col] = self.matrix[row][0]*self.matrix[row][col];
            }
        }

        HM{size : self.size, matrix : new_mat}
    }






    pub fn to_string(&self) -> String {

        let mut result = "".to_string();

        for row in 0..self.size {
            result += "| ";
            for col in 0..self.size {
                if self.matrix[row][col] == 1 {
                    result += &"+ ";
                }
                else {
                    result += &"- "
                }
            }
            result += "|\n";


        }

        result
    }

    pub fn to_string_magma(&self) -> String {

        let mut result = "M := Matrix([".to_string();

        for row in 0..self.size-1 {
            result += "[";
            for col in 0..self.size-1 {
                if self.matrix[row][col] == 1 {
                    result += &"1,";
                }
                else {
                    result += &"-1,"
                }
            }
            if self.matrix[row][self.size-1] == 1 {
                result += &"1";
            }
            else {
                result += &"-1"
            }


            result += "],";
        }


        result += "[";
        for col in 0..self.size-1 {
            if self.matrix[self.size-1][col] == 1 {
                result += &"1,";
            }
            else {
                result += &"-1,"
            }
        }
        if self.matrix[self.size-1][self.size-1] == 1 {
            result += &"1";
        }
        else {
            result += &"-1"
        }

        result += "]";

        result += &"]);";

        result
    }

    // Returns reference to specified row (0 indexed)
    pub fn row(&self, row : usize) -> Vec<i8> {
        self.matrix[row].clone()
    }

    // Verifies Hadamard matrix property
    pub fn verify(&self) -> bool {
        let n = self.size();

        // Take inner product of each row with every other row
        for row1 in 0..n {
            for row2 in (row1+1)..n {
                // Get two rows
                let row1_vec = self.row(row1);
                let row2_vec = self.row(row2);

                // Perform inner product
                let mut result = 0;
                for i in 0..n {
                    result += row1_vec[i] * row2_vec[i];
                }

                if result != 0 {
                    return false;
                }
            }
        }

        return true;
    }

}
