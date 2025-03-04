#![allow(unused_imports)]

use beralg::lattice::methods::timing::{enumeration_times_with_lll, plot_distance_diffs, plot_shortest_vector};

fn main() {
    let start = 7;
    // let end = 21;
    // plot_distance_diffs(start, end);
    let end = 26;
    plot_shortest_vector(start, end);
    let end = 19;
    enumeration_times_with_lll(start, end);

   // check_closest_vector_by_enumeration_limit(); 

}
