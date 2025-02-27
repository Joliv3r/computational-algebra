#![allow(unused_imports)]
use ndarray::{array, Array1, Array2, Slice};
use beralg::lattice::{self, methods::gram_schmidt, Lattice};
extern crate openblas_src;

fn main() {
    let matrix: Vec<Array1<f64>> = vec![
        array![902., 432., 83.],
        array![4., -5., 83.],
        array![5., 2., 2.],
    ];

    let lattice = Lattice::build_lattice_basis_from_vectors(&matrix).unwrap();
    let sbv = lattice.get_shortest_basis_vector().unwrap();

    println!("{:#?}", sbv);

    let svp = lattice.shortest_vector_by_enumeration().unwrap();
    println!("{}", svp);
}
