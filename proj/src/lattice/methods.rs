use std::ops::Mul;

use itertools::Itertools;
use ndarray::{Array1, Array2};
use ndarray_linalg::Solve;

use super::Lattice;
pub mod closest_vector;
pub mod shortest_vector;
pub mod basis_reduction;
pub mod timing;

// We expect this function to be called only with vectors of the same size.
pub fn make_matrix_from_column_vectors(vectors: &Vec<Array1<f64>>) -> Array2<f64> {
    let columns = vectors.len();
    let rows = vectors[0].shape()[0].clone();
    let flattened = vectors.into_iter().flatten().cloned().collect();

    Array2::from_shape_vec((columns, rows), flattened).unwrap().reversed_axes()
}

// TODO make sure this is a basis and not just a span of some vectors.
// For the moment we use only full rank lattices as a workaround.
pub fn make_into_basis_matrix(vectors: &Vec<Array1<f64>>) -> Result<Array2<f64>, String> {
    if !vectors.iter().map(|v| v.len()).all_equal() {
        return Err("All basis vectors are not of the same length.".to_string())
    }
    Ok(make_matrix_from_column_vectors(vectors))
}

// Workaround to ensure vectors are linearly independent.
// TODO: Make this more effective.
pub fn is_linearly_independent(vectors: &Vec<Array1<f64>>) -> bool {
    if let Ok(matrix) = make_into_basis_matrix(vectors) {
        let zero = Array1::zeros(matrix.shape()[0]);
        if let Err(_) = matrix.solve(&zero) {
            false
        } else {
            true
        }
    } else {
        false
    }
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
        orthogonal.push(v.clone());
    }
    orthogonal
}


fn collect_columns_in_vec(matrix: &Array2<f64>) -> Vec<Array1<f64>> {
    let mut vec: Vec<Array1<f64>> = Vec::with_capacity(matrix.rows().into_iter().len());
    
    matrix.map_axis(ndarray::Axis(0), |v| 
        vec.push(Array1::from_vec(v.to_vec()))
    );

    vec
}


pub fn gram_schmidt_columns(matrix: &Array2<f64>) -> Array2<f64> {
    let vec = collect_columns_in_vec(matrix);
    make_matrix_from_column_vectors(&gram_schmidt(&vec))
}


pub fn get_length_of_vector(vector: &Array1<f64>) -> f64 {
    vector.dot(vector)
}


impl Lattice {
    // Calculates mu_kj = <b_k, b*_j>/<b*_j, b*_j>
    fn get_mu(&self, k: usize, j: usize) -> f64 {
        let b_k = self.get_basis_vector(k).expect("Should exist.");
        let b_orth_j = self.get_gram_schmidt_basis_vector(j).expect("Should exist.");
        b_k.dot(&b_orth_j)/(b_orth_j.dot(&b_orth_j))
    }

    fn get_gram_schmidt_length(&self, index: usize) -> f64 {
        let b = self.get_gram_schmidt_basis_vector(index).expect("Should exist");
        b.dot(&b)
    }

    fn write_vector_with_gram_schmidt_vectors(&self, vector: &Array1<f64>) -> Array1<f64> {
        let gram_schmidt_matrix = make_into_basis_matrix(&self.gram_schmidt_basis).expect("Gram-Schmidt basis should be a square matrix.");
        gram_schmidt_matrix.solve(vector).expect("Gram-Schmidt matrix is of full rank.")
    }

    fn write_vector_with_basis_vectors(&self, vector: &Array1<f64>) -> Array1<f64> {
        let basis_matrix = make_into_basis_matrix(&self.basis).expect("Basis matrix should be a square matrix.");
        basis_matrix.solve(&vector).expect("Basis matrix is of full rank.")
    }
}


#[cfg(test)]
mod lattice_tests {
    use rand::{thread_rng, Rng};
    use ndarray::array;

    use super::*;
    use super::timing::{generate_random_basis, generate_random_vector};

    #[test]
    fn test_matrix_from_column_vectors() {
        let vectors = vec![array![1.,2.,3.], array![4.,5.,6.], array![7.,8.,9.]];
        let matrix = array![
            [1.,4.,7.],
            [2.,5.,8.],
            [3.,6.,9.],
        ];
        let made_matrix = make_matrix_from_column_vectors(&vectors);
        let made_vectors = collect_columns_in_vec(&matrix);
        assert_eq!(made_matrix, matrix);
        assert_eq!(made_vectors, vectors);
    }

    #[test]
    #[allow(unused_must_use)]
    fn test_lattice_building() {
        let vectors = vec![array![1.,2.,3.], array![4.,5.,6.], array![7.,8.,8.]];

        let mut lattice = Lattice::build_lattice_basis_from_vectors(&vectors).unwrap();
        let shortest_vector = array![1.,2.,3.];

        assert_eq!(lattice.get_shortest_basis_vector().unwrap(), shortest_vector.view());

        let combination = vec![4, 6, -2];

        assert_eq!(lattice.get_lattice_point(&combination).unwrap(), array![14., 22., 32.]);

        let new_vectors = vec![
            array![7.,8.,8.],
            array![4.,5.,6.]
        ];
        
        lattice.swap_basis_vectors(1, 2);
        assert_eq!(lattice.get_basis_columns(1, 3).unwrap(), new_vectors);

        lattice.update_basis_vector(2, &shortest_vector);
        assert_eq!(lattice.get_basis_vector(2).unwrap(), shortest_vector);
    }

    #[test]
    fn test_simple_lattice_operations() {
        let mut rng = thread_rng();
        let dimension = rng.gen_range(2..14);
        let basis = generate_random_basis(dimension);
        let lattice = Lattice::build_lattice_basis_from_vectors(&basis).unwrap();

        let index = rng.gen_range(0..dimension);
        let column = lattice.get_basis_vector(index).unwrap();
        let correct_column = basis.get(index).unwrap();

        assert_eq!(column, correct_column);
        assert_eq!(lattice.get_length_of_basis_vectors(), dimension);
        assert_eq!(lattice.columns(), dimension);
        assert!(lattice.index_exists(dimension-1));
        assert!(!lattice.index_exists(dimension));
    }

    #[test]
    fn test_gram_schmidt() {
        let dimension = 100;
        let tol = 1e-9;
        let basis = generate_random_basis(dimension);
        let b = gram_schmidt(&basis);
        
        for i in 0..7 {
            let b_i = b.get(i).unwrap();
            for j in i+1..7 {
                let b_j = b.get(j).unwrap();
                let ip = b_i.dot(b_j);
                assert!(ip.abs() < tol, "ip = {}, between index {}, and {}", ip, i, j);
            }
        }
    }

    #[test]
    #[allow(unused_must_use)]
    fn problem_solvers_runs_without_crashing() {
        let dimension = 5;
        let basis = generate_random_basis(dimension);
        let vector = generate_random_vector(dimension, 1.);
        let mut lattice = Lattice::build_lattice_basis_from_vectors(&basis).unwrap();

        let delta = 0.75;
        lattice.lll_reduction(delta);

        lattice.shortest_vector_by_enumeration();
        lattice.babai_nearest_plane(&vector);
        lattice.closest_vector_by_enumeration(&vector);
    }

    #[test]
    fn test_linearly_independent() {
        let mut vectors = vec![
            array![1., 4., 7.],
            array![2., 5., 8.],
            array![3., 6., 9.]
        ];

        assert!(!is_linearly_independent(&vectors));

        vectors[2] = array![3., 6., 8.];
        assert!(is_linearly_independent(&vectors));

        vectors[1] = array![1., 5.];
        assert!(!is_linearly_independent(&vectors));
    }

    #[test]
    fn test_cvp_by_enumeration() {
        let rng = thread_rng();
        let dimension = 15;
        let loops = 20;
        let basis = generate_random_basis(dimension);
        let vector = generate_random_vector(dimension, 1.);

        let mut lattice = Lattice::build_lattice_basis_from_vectors(&basis).unwrap();
        lattice.lll_reduction(0.75);

        for _ in 0..loops {
            let cvp_enumeration = lattice.closest_vector_by_enumeration(&vector).unwrap();
            let cvp_babai = lattice.babai_nearest_plane(&vector).unwrap();

            let len_enum = get_length_of_vector(&(&cvp_enumeration - &vector));
            let len_babai = get_length_of_vector(&(&cvp_babai - &vector));

            if len_enum > len_babai {
                panic!("Enumeration did not find the closest vector.");
            }
        }
    }
}
