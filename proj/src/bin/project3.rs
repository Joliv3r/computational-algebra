#![allow(unused_imports)]
use beralg::lattice::methods::get_length_of_vector;
use beralg::lattice::methods::timing::check_closest_vector_by_enumeration_limit;
use beralg::lattice::methods::timing::cvp_statistics;
use beralg::lattice::methods::timing::enumeration_times_with_lll;
use beralg::lattice::methods::timing::increase_basis;
use beralg::lattice::methods::timing::plot_distance_diffs;
use beralg::lattice::methods::timing::plot_shortest_vector;
use beralg::lattice::Lattice;
use itertools::Itertools;
use ndarray::array;
use ndarray::Array2;
use ndarray_linalg::Solve;
use rand::thread_rng;

fn main() {
    // let top = 20;
    // cvp_statistics(top);

    // let top = 50;
    // babai_times(top);

    let start = 2;
    let end = 17;
    plot_distance_diffs(start, end);
    let end = 25;
    plot_shortest_vector(start, end);
    let end = 19;
    enumeration_times_with_lll(end);

   // check_closest_vector_by_enumeration_limit(); 

}
