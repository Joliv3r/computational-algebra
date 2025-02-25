use ndarray::{array, Array1, Array2};
use ndarray_linalg::Solve;
use beralg::lattice;
extern crate openblas_src;

fn main() {
    let matrix: Array2<f64> = array![
        [3., 2., -1., 1.],
        [2., -2., 4., -4.],
        [-2., 1., -2., 2.]
    ];

    let zero = Array1::zeros(3);
    if let Ok(sol) = matrix.solve(&zero) {
        println!("{}", sol);
    } else {
        println!("Error");
    }
    //
    // let a: Array1<f64> = array![1., -1., 0.]; 
    // let p = Array1::from_vec(a.iter().map(|x| 5.*x).collect());
    // println!("{}", p);
    // let b: Array1<f64> = array![1., -0., 0.]; 
    //
    // let s = a-b;
    // println!("{}", s);
    //
    // if let Ok(sol) = matrix.solve(&b) {
    //     println!("{}", sol);
    // } else {
    //     println!("Error");
    // }

    // let a = vec![vec![1,2,3], vec![4,5,6], vec![7,8,9]];
    // let m = Array2::from(a);
    // println!("{}", m);

    // let mut s = vec![array![3.0, 1.0], array![2.0, 2.0]];
    // let gs = beralg::lattice::methods::gram_schmidt(&s);
    // let ip = gs.get(0).unwrap().dot(gs.get(1).unwrap());
    // println!("{}", ip);
    // println!("{:#?}", gs);
}
