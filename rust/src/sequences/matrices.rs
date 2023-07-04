use cgmath::Quaternion;

use super::sequence::{QS, quaternion_to_string};



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

    pub fn dephase(&mut self) {
        // dephase the rows
        for inv_row in 0..self.size {
            for col in 0..self.size {
                let row = self.size - 1 - inv_row;
                self.matrix[row][col] = self.matrix[0][col].conjugate()*self.matrix[row][col];
            }
        }
        println!("{}", self.to_string());
        // dephase the columns
        for inv_col in 0..self.size {
            for row in 0..self.size {
                let col = self.size - 1 - inv_col;
                self.matrix[row][col] = self.matrix[row][0].conjugate()*self.matrix[row][col];
            }
        }
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