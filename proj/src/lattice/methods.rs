use itertools::Itertools;
use ndarray::{array, stack, Array, Array1, Array2, ArrayView2, Ix2, LinalgScalar};

// We expect this function to be called only with vectors of the same size.
pub fn make_matrix_from_vectors<T: Clone>(vectors: &Vec<Vec<T>>) -> Array2<T> {
    let columns = vectors.len();
    let rows = vectors[0].len();
    let flattened = vectors.into_iter().flatten().cloned().collect::<Vec<T>>();

    Array2::from_shape_vec((columns, rows), flattened).unwrap().reversed_axes()
}

// TODO make sure this is a basis and not just a spanning of some vectors.
// For the moment we use only full rank lattices as a workaround.
pub fn make_into_basis_matrix<T: Clone>(vectors: &Vec<Vec<T>>) -> Result<Array2<T>, String> {
    if !vectors.iter().map(|v| v.len()).all_equal() {
        return Err("All basis vectors are not of the same length.".to_string())
    }
    Ok(make_matrix_from_vectors(vectors))
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


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_make_matrix_from_vectors() {
        let a = vec![vec![1,2,3], vec![4,5,6], vec![7,8,9]];
        let m = make_matrix_from_vectors(&a);
        let result = array![
            [1,4,7],
            [2,5,8],
            [3,6,9],
        ];

        assert_eq!(m, result);
    }
}
