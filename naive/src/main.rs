#[macro_use]
extern crate lazy_static;

use time::*;

pub mod sequence;
pub mod find;
pub mod symmetries;
use crate::find::find;
use crate::symmetries::*;

fn main() {
    let mut now;
    let mut elapsed_time;
    let mut count;

    for i in 1..18{
        now = Instant::now();
        count = find(i, Some(Symmetry::I));
        elapsed_time = now.elapsed().as_seconds_f32();

        println!("For n = {i}, the function took: {elapsed_time} seconds and found {count} sequences");
    }

}
