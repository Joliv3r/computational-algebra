#![allow(unused_imports)]
use beralg::lattice::methods::get_length_of_vector;
use beralg::lattice::methods::timing::cvp_statistics;
use beralg::lattice::Lattice;
use itertools::Itertools;
use ndarray::array;
use ndarray::Array2;
use ndarray_linalg::Solve;

fn main() {
    // let basis = vec![
    //     array![329.,234.,643.,],
    //     array![5., 7., 2.,],
    //     array![164., 117., 321.,],
    // ];
    //
    // let lattice = Lattice::build_lattice_basis_from_vectors(&basis).unwrap();
    //
    // if let Ok(shortest_vector) = lattice.shortest_vector_by_enumeration() {
    //     println!("{}", shortest_vector);
    // }
    // let basis = vec![
    //     array![232., -12.2],
    //     array![-134., 11.1],
    // ];
    // let basis = vec![
    //     array![1., 1., 0.,],
    //     array![0., 1., 0.,],
    //     array![0., 0., 1.,]
    // ];
    //
    // let mut lattice = Lattice::build_lattice_basis_from_vectors(&basis).unwrap();
    // lattice.lll_reduction(0.75);
    // lattice.print_basis();
    // lattice.print_gram_schmidt_basis();
    // let a = array![-342.01, 2134.12];
    //
    // let (cvp, checked_points) = lattice.closest_vector_by_enumeration(&a).unwrap();
    // let babai = lattice.babai_nearest_plane(&a).unwrap();
    // let (cvp, checked_points) = lattice.shortest_vector_by_enumeration().unwrap();
    // println!("{:#?}", checked_points);
    // println!("{}", cvp);
    // println!("{}", babai);
    // let cvp_dist = get_length_of_vector(&(&a - &cvp));
    // let babai_dist = get_length_of_vector(&(&a - &babai));
    // println!("{}", cvp_dist);
    // println!("{}", babai_dist);
    // let combination = vec![6];
    // let current_vector = lattice.write_vector_with_gram_schmidt_vectors_from_reversed_basis_representation(&combination);
    // 
    // let point = (a[0], a[1]);
    // let size = 100;
    // let plot_size = 100;
    // lattice.print_lattice_around_point(point, plot_size, size, &checked_points);

    let top = 10;
    cvp_statistics(top);
}
