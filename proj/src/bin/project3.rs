use ndarray::{array, Array1, Array2};
use ndarray_linalg::Solve;
use beralg::lattice::{self, methods::{babai_nearest_plane, gram_schmidt}, Lattice};
extern crate openblas_src;
use rand::{thread_rng, Rng};

fn main() {
    let matrix: Array2<f64> = array![
        [1., 0., 0.],
        [0., 1., 0.],
        [0., 0., 1.],
    ];

    let lattice = Lattice{ basis: matrix };

    let v = array![134.5, 3.5, 4.49];

    let p = match babai_nearest_plane(&lattice, &v) {
        Ok(p) => println!("{}", p),
        Err(e) => println!("Error: {}", e),
    };

    let mut rng = thread_rng();
    let mut matrix: Array2<f64> = Array2::zeros((7,7));
    matrix.map_inplace(|mut e| {*e = rng.gen_range(0..100) as f64/(rng.gen_range(1..100) as f64)});

    // let zero = Array1::zeros(3);
    // if let Ok(sol) = matrix.solve(&zero) {
    //     println!("{}", sol);
    // } else {
    //     println!("Error");
    // }
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
