
#[derive(Clone)]
pub enum Symmetry{ // enum for the different types of Quaternion Sequences
    I, II, III, IV
}

#[derive(Clone)]
pub enum RowsumPairing{
    XY, XZ, XW
}

#[derive(Clone)]
pub enum SequenceType{ // enum for the different types of Quadruplets of sequences
    Williamson, WilliamsonType, ItoType, ExtraTypeI, ExtraTypeII, ExtraTypeIII
}