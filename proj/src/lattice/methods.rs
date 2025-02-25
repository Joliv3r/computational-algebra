use itertools::Itertools;
use ndarray::{array, Array1, Array2};

use super::Lattice;

// We expect this function to be called only with vectors of the same size.
pub fn make_matrix_from_column_vectors(vectors: &Vec<Array1<f64>>) -> Array2<f64> {
    let columns = vectors.len();
    let rows = vectors[0].shape()[0].clone();
    let flattened = vectors.into_iter().flatten().cloned().collect();

    Array2::from_shape_vec((columns, rows), flattened).unwrap().reversed_axes()
}

// TODO make sure this is a basis and not just a spanning of some vectors.
// For the moment we use only full rank lattices as a workaround.
pub fn make_into_basis_matrix(vectors: &Vec<Array1<f64>>) -> Result<Array2<f64>, String> {
    if !vectors.iter().map(|v| v.len()).all_equal() {
        return Err("All basis vectors are not of the same length.".to_string())
    }
    Ok(make_matrix_from_column_vectors(vectors))
}


// Uses the Gram-Schmidt algorithm, which gives an orthogonal set {v_1, ..., v_t} to a set {u_1, ..., u_t}
// where
//   u_i = v_i - sum_{j=1,...,i-1}( <u_j, v_i> (v_i) )
pub fn gram_schmidt(vectors: &Vec<Array1<f64>>) -> Vec<Array1<f64>> {
    let mut orthogonal = Vec::with_capacity(vectors.len());
    for i in 0..vectors.len() {
        let u = vectors.get(i).expect("Index should exists");
        let mut v = u.clone();
        for j in 0..i {
            let u_j: &Array1<f64> = orthogonal.get(j).expect("Index should exists");
            let scalar = u.dot(u_j)/u_j.dot(u_j);
            let diff = u_j*scalar;

            v = v - diff;
        }
        // println!("{}", u);
        orthogonal.push(v.clone());
    }
    orthogonal
}


fn collect_columns_in_vec(matrix: &Array2<f64>) -> Vec<Array1<f64>> {
    let mut vec: Vec<Array1<f64>> = Vec::with_capacity(matrix.rows().into_iter().len());
    
    matrix.map_axis(ndarray::Axis(1), |v| 
        vec.push(Array1::from_vec(v.to_vec()))
    );

    vec
}


fn gram_schmidt_columns(matrix: &Array2<f64>) -> Array2<f64> {
    let vec = collect_columns_in_vec(matrix);
    make_matrix_from_column_vectors(&gram_schmidt(&vec))
}


pub fn babai_nearest_plane(lattice: &Lattice, v: &Array1<f64>) -> Result<Array1<f64>, String> {
    let cols = lattice.get_length_of_basis_vectors();
    if cols != v.len() {
        return Err("Vector size not compatible with lattice".to_string())
    }

    let gram_schmidt_basis = gram_schmidt_columns(&lattice.basis);
    let mut w = v.clone();
    let dim = lattice.get_number_of_basis_vectors();
    let mut y = Array1::zeros(cols);

    for i in (0..dim).into_iter().rev() {
        let gs_i = Array1::from_vec(gram_schmidt_basis.column(i).to_vec());
        let b_i = lattice.get_basis_vector(i);
        let l_i = &w.dot(&gs_i)/gs_i.dot(&gs_i);
        let l_i_rounded = l_i.round();

        y = y + &b_i*l_i_rounded;

        w = w - gs_i*(l_i - l_i_rounded) - b_i*l_i_rounded;
    }

    Ok(y)
}


#[cfg(test)]
mod lattice_tests {
    use rand::{thread_rng, Rng};

    use super::*;

    #[test]
    fn test_make_matrix_from_column_vectors() {
        let a = vec![array![1.,2.,3.], array![4.,5.,6.], array![7.,8.,9.]];
        let m = make_matrix_from_column_vectors(&a);
        let result = array![
            [1.,4.,7.],
            [2.,5.,8.],
            [3.,6.,9.],
        ];

        assert_eq!(m, result);
    }


    #[test]
    fn test_gram_schmidt() {
        let tol = 1e-9;
        let mut rng = thread_rng();
        let mut matrix: Array2<f64> = Array2::zeros((7,7));
        matrix.map_inplace(|mut e| {*e = rng.gen_range(0..100) as f64/(rng.gen_range(1..100) as f64)});
        let b = gram_schmidt_columns(&matrix);
        
        for i in 0..7 {
            for j in i+1..7 {
                let ip = b.column(i).dot(&b.column(j));
                assert!(ip.abs() < tol, "ip = {}, between index {}, and {}", ip, i, j);
            }
        }
    }
}
