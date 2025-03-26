use bit_matrix::BitMatrix;
use gauss_jordan_elimination::gauss_jordan_elimination_generic;
use ndarray::{array, linalg::Dot, Array1, Array2, ArrayBase};
use ndarray_linalg::LeastSquaresSvd;
use beralg::algebraic_structure::z2::Z2;

fn main() {
    let mut matrix: Vec<Vec<Z2>> = vec![
        vec![Z2(true),Z2(true),Z2(false)],
        vec![Z2(true),Z2(false),Z2(true)],
        vec![Z2(false),Z2(true),Z2(true)]
    ];

    // let b: Array1<Z2> = array![Z2(true), Z2(true), Z2(false)];
    
    gauss_jordan_elimination_generic(&mut matrix);

    // let mut matrix = BitMatrix::new(3, 3);

    // let points = &[
    //     (0, 0),
    //     (1, 1),
    // ];
    //
    // for &(i, j) in points {
    //     matrix.set(i, j, true);
    // }
    //

    println!("{:#?}", matrix);
}
