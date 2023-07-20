use std::ops;

use cgmath::Quaternion;
use itertools::iproduct;
use num_complex::Complex;

use super::matrices::QHM;


pub fn make_operator(a : (f32,f32), b : (f32,f32), c : (f32,f32), d : (f32,f32)) -> Operator {
    let size = 2;
    let values = vec![vec![Complex {re : a.0 , im : a.1}, Complex {re : b.0 , im : b.1}],vec![Complex {re : c.0 , im : c.1},Complex {re : d.0 , im : d.1}]];
    Operator {size, values}
}



lazy_static! {
    
    pub static ref OP0 : Operator = make_operator((0.,0.), (0.,0.), (0.,0.), (0.,0.));
    pub static ref OP1 : Operator = make_operator((1.,0.), (0.,0.), (0.,0.), (1.,0.));
    pub static ref OPX : Operator = make_operator((0.,0.), (1.,0.), (1.,0.), (0.,0.));
    pub static ref OPY : Operator = make_operator((0.,0.), (0.,-1.), (0.,1.), (0.,0.));
    pub static ref OPZ : Operator = make_operator((1.,0.), (0.,0.), (0.,0.), (-1.,0.));

}




pub fn quaternion_to_operator(quat : &Quaternion<f32>) -> Operator {
    let (s,x,y,z) = (quat.s, quat.v.x, quat.v.y, quat.v.z);

    OP1.clone()*Complex::new(s, 0.) + &(OPX.clone()*Complex::new(0., x)) + &(OPY.clone()*Complex::new(0., -y)) + &(OPZ.clone()*Complex::new(0., z))
}


pub fn computational_basis(size : usize, j : usize, k : usize) -> Operator {
    let mut values = vec![vec![Complex::new(0.,0.); size]; size];
    values[j][k] = Complex::new(1.,0.);

    Operator {size, values}
}


#[derive(Clone, PartialEq, Debug)]
pub struct Operator {
    size : usize,
    values : Vec<Vec<Complex<f32>>>
}

impl Operator {
    pub fn size(&self) -> usize {
        self.size
    }

    pub fn empty(size : usize) -> Operator {
        Operator { size, values : vec![vec![Complex::new(0.,0.); size]; size]}
    }

    pub fn conjugate_transpose(&self) -> Operator {
        let mut values = vec![];
        
        for i in 0..self.size() {
            let mut line = vec![];
            for j in 0..self.size() {
                // We transpose AND conjugate
                line.push(self.values[j][i].conj())
            }
            values.push(line)
        }
        Operator {size : self.size(), values}
    }

    pub fn tensor(&self, other : &Operator) -> Operator {
        let size = self.size * other.size();
        let mut values = vec![vec![Complex::new(0.,0.); size]; size];

        for (i,j) in iproduct!(0..self.size(), 0..self.size()) {
            for (p,q) in iproduct!(0..other.size(), 0..other.size()) {
                values[i*other.size + p][j*other.size + q] = self.values[i][j] * other.values[p][q];
            }
        }

        Operator {size, values}
    }

    pub fn to_string(&self) -> String {

        let mut result = "".to_string();

        for row in 0..self.size {
            result += "| ";
            for col in 0..self.size {
                result += &(self.values[row][col].to_string() + &" ");
            }
            result += "|\n";

        }

        result

    }

}

impl ops::Add<&Operator> for Operator {
    type Output = Operator;
    
    fn add(self, rhs: &Operator) -> Self::Output {
        let mut values = vec![];
        for i in 0..self.size() {
            let mut line = vec![];
            for j in 0..self.size() {
                line.push(self.values[i][j] + rhs.values[i][j])
            }
            values.push(line);
        }
        Operator { size : self.size, values}
    }
}

impl ops::Mul<Complex<f32>> for Operator {
    type Output = Operator;
    
    fn mul(self, rhs: Complex<f32>) -> Self::Output {
        let mut values = vec![];
        for i in 0..self.size() {
            let mut line = vec![];
            for j in 0..self.size() {
                line.push(self.values[i][j] * rhs)
            }
            values.push(line);
        }
        Operator { size : self.size, values}
    }
}



#[derive(PartialEq, Debug)]
pub struct MUM {
    matrix_size : usize,
    sequence : Vec<Operator>
}


impl MUM {
    
    pub fn new(matrix_size : usize, sequence : &Vec<Operator>) -> MUM {
        MUM {matrix_size, sequence : sequence.clone()}
    }


    pub fn empty(length : usize, matrix_size : usize) -> MUM {
        let mat = Operator::empty(matrix_size);
        let sequence = vec![mat.clone();length];
        MUM {matrix_size, sequence}
    }

    pub fn matrix_size(&self) -> usize {
        self.matrix_size
    }

    pub fn length(&self)-> usize {
        self.sequence.len()
    }

    pub fn values(&self)-> Vec<Operator> {
        self.sequence.clone()
    }

    pub fn from_qhm(qhm : &QHM) -> MUM {
        let matrix_size = qhm.size()*2;
        let mut sequence = vec![];

        for b in 0..qhm.size() {
            let mut q_b = Operator::empty(matrix_size);

            for (j,k) in iproduct!(0..qhm.size(), 0..qhm.size()) {
                let quat = qhm.get(b,j) * qhm.get(b,k).conjugate();
                let basis_elm = computational_basis(qhm.size(), j, k);
                q_b = q_b + &quaternion_to_operator(&quat).tensor(&basis_elm);
            }
            // we normalize by 1/size
            sequence.push(q_b * Complex::new(1./matrix_size as f32, 0.));
        }

        MUM {matrix_size, sequence}
    }

    pub fn to_string(&self) -> String {
        let mut result = "".to_string();

        for mat in &self.sequence {
            result += &(mat.to_string() + &"\n\n");
        }

        result
    }

}