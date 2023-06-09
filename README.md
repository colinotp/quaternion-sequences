# Description

This repository will contain various approaches and algorithms to find specific quaternion sequences
The python folder contains python code, and the rust folder contains rust code.


## In the rust folder:
### The results folder contains all the data generated by the code.

The results folder has multiple folders corresponding to the different types of sequences the search has produced.
Here's what each of them stand for:
pqs     -> perfect quaternion sequence
ws      -> Williamson sequence
wts     -> Williamson-type sequence


The i.seq files contain a list of all the sequences of the specified type, of length i, starting with a 1.
The i.log contains information about the performances of the program when computing the i.seq file.


The Quaternion and Williamson sequences are stocked in the following format:
each line contains a sequence, and each element correspond to a character:
1       -> +
-1      -> -
i       -> i
j       -> j
k       -> k
q       -> q
qi     -> x
qj     -> y
qk     -> z

with q = (1+i+j+k)/2

The negatives of these quaternions (except 1 and -1) are stocked as the capital letter representing the quaternion.
For example, -qi -> X



### The src folder contains all of the code 

The main details what part of the code is running. You can change that by changing which function of the file it is running.
You can run the code with the command 'cargo run'

The find folder contains the code that finds and generates specific sequences

the sequences folder contains the code for the classes related to Quaternion or Williamson sequences.

The test folder contains tests of various parts of the code.
You can run the tests with the command 'cargo test'