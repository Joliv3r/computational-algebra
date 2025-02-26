use std::thread::current;

use itertools::Itertools;
use ndarray::{array, Array1, Array2, ArrayView1, Slice};

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
    
    matrix.map_axis(ndarray::Axis(0), |v| 
        vec.push(Array1::from_vec(v.to_vec()))
    );

    vec
}


pub fn gram_schmidt_columns(matrix: &Array2<f64>) -> Array2<f64> {
    let vec = collect_columns_in_vec(matrix);
    make_matrix_from_column_vectors(&gram_schmidt(&vec))
}


impl Lattice {
    pub fn babai_nearest_plane(&self, v: &Array1<f64>) -> Result<Array1<f64>, String> {
        let rows = self.get_length_of_basis_vectors();
        if rows != v.len() {
            return Err("Vector size not compatible with lattice".to_string())
        }

        let mut w = v.clone();
        let dim = self.columns();
        let mut y = Array1::zeros(rows);

        for i in (0..dim).into_iter().rev() {
            let gs_i = Array1::from_vec(self.get_gram_schmidt_basis_vector(i).to_vec());
            let b_i = Array1::from_vec(self.get_basis_vector(i).to_vec());
            let l_i = &w.dot(&gs_i)/gs_i.dot(&gs_i);
            let l_i_rounded = l_i.round();

            y = y + &b_i*l_i_rounded;

            w = w - gs_i*(l_i - l_i_rounded) - b_i*l_i_rounded;
        }

        Ok(y)
    }


    // Say v = sum_{i=1}^{n} x_j b_j, is the shortest vector. We use that
    //   -(M_1 + M_2) =< x_i => M_1 - M_2
    // where
    //   M_1 = sqrt{ ( A - sum_{j = i+1}^{n} x_j^2 B_j )/B_i }
    //   M_2 = sum_{j = i+1}^{n} µ_{j,i} x_j
    // with A > ||v||^2 and B_j = ||b_j||^2 and µ_{j,i} = <b_i, b*_j>/||b*_j||^2.
    pub fn shortest_vector_by_enumeration(&self) -> Array1<f64> {
        let mut shortest_vector: Array1<f64> = self.get_shortest_basis_vector().to_vec().into();
        let mut shortest_length = shortest_vector.dot(&shortest_vector);

        let b_orth_n = self.get_gram_schmidt_basis_vector(self.columns()-1);
        let bound_n = (shortest_length/b_orth_n.dot(&b_orth_n)).sqrt().floor() as i32;

        println!("Bound of n is: {}", bound_n);


        for x_n in 0..bound_n {
            let combination = vec![x_n];
            let (candidate_vector, candidate_length) = self.shortest_vector_enumeration_steps(1, &combination, &shortest_vector.to_vec().into(), shortest_length);
            if candidate_length < shortest_length && candidate_length != 0. {
                shortest_vector = candidate_vector;
                shortest_length = candidate_length;
            }
        }

        shortest_vector
    }

    fn shortest_vector_enumeration_steps(&self, depth: usize, mut combination: &Vec<i32>, current_shortes_vector: &Array1<f64>, current_shortest_length: f64) -> (Array1<f64>, f64) {
        println!("Checking combination: {:#?}", combination);
        if depth == self.columns() {
            let mut current_vector = Array1::zeros(self.get_length_of_basis_vectors());
            for (x_i, b_i) in combination.iter().zip(self.basis.columns().into_iter().rev()) {
                current_vector = current_vector + b_i.to_owned() * (*x_i as f64);
            }
            let current_length = current_vector.dot(&current_vector);
            return (current_vector, current_length)
        }

        let mut shortest_vector = current_shortes_vector.clone();
        let mut shortest_length = current_shortest_length.to_owned();
        let (lower_bound, upper_bound) = self.get_enumeration_bounds(combination, depth, shortest_length);
        println!("We have the bounds: {}, {}", lower_bound, upper_bound);
        for i in lower_bound..=upper_bound {
            let mut new_combination = combination.clone();
            new_combination.push(i);
            let (candidate_vector, candidate_length) = self.shortest_vector_enumeration_steps(depth+1, &new_combination, &shortest_vector, shortest_length);
            if candidate_length < shortest_length {
                shortest_vector = candidate_vector;
                shortest_length = candidate_length;
            }
        }

        (shortest_vector, shortest_length)
    }

    #[allow(non_snake_case)]
    fn get_enumeration_bounds(&self, combination: &Vec<i32>, basis_number: usize, A: f64) -> (i32, i32) {
        let mut sum: f64 = 0.;
        let mut M_2: f64 = 0.;
        let b_i = self.get_basis_vector(self.columns()-basis_number);
        let slice: Slice = Slice::from(self.columns()-basis_number-1..self.columns());
        for (x_j, (b_j, b_orth_j)) in combination.iter().zip(self.get_basis_columns(slice).columns().into_iter().zip(self.get_gram_schmidt_basis_columns(slice).columns().into_iter()).rev()) {
            println!("Summing with: {}, {}, {}", x_j, b_j, b_orth_j);
            let B_j = b_orth_j.dot(&b_orth_j);
            let mu_ij = b_i.dot(&b_orth_j)/b_orth_j.dot(&b_orth_j);
            sum += (x_j.pow(2) as f64)*B_j;
            M_2 += mu_ij*(*x_j as f64);
        }
        println!("We have: A={}, sum={}", A, sum);
        let M_1 = ((A-sum)/b_i.dot(&b_i)).sqrt();
        println!("We got: M1 = {}, M2 = {}", M_1, M_2);

        (-(M_1+M_2).ceil() as i32, (M_1-M_2).floor() as i32)
    }
}


#[cfg(test)]
mod lattice_tests {
    use rand::{thread_rng, Rng};

    use super::*;

    fn generate_random_matrix(dimension: usize) -> Array2<f64> {
        let mut rng = thread_rng();
        let mut matrix: Array2<f64> = Array2::zeros((dimension,dimension));
        matrix.map_inplace(|mut e| {*e = rng.gen_range(0..100) as f64/(rng.gen_range(1..100) as f64)});
        matrix
    }

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
    fn test_lattice_building() {
        let vectors = vec![array![1.,2.,3.], array![4.,5.,6.], array![7.,8.,9.]];
        let basis = array![
            [1.,4.,7.],
            [2.,5.,8.],
            [3.,6.,9.]
        ];

        let m1 = Lattice::build_from_basis(&basis);
        let m2 = Lattice::build_lattice_basis_from_vectors(&vectors).unwrap();
        let shortest_vector = array![1.,2.,3.];

        assert_eq!(m1.basis, m2.basis);
        assert_eq!(m1.get_shortest_basis_vector(), shortest_vector.view());
    }

    #[test]
    fn test_simple_lattice_operations() {
        let mut rng = thread_rng();
        let dimension = rng.gen_range(2..14);
        let basis = generate_random_matrix(dimension);
        let lattice = Lattice::build_from_basis(&basis);

        let index = rng.gen_range(0..dimension);
        let column = lattice.get_basis_vector(index);
        let correct_column = basis.column(index);

        assert_eq!(column, correct_column);
        assert_eq!(lattice.get_length_of_basis_vectors(), dimension);
        assert_eq!(lattice.columns(), dimension);

    }

    #[test]
    fn test_gram_schmidt() {
        let dimension = 7;
        let tol = 1e-9;
        let matrix = generate_random_matrix(dimension);
        let b = gram_schmidt_columns(&matrix);
        
        for i in 0..7 {
            for j in i+1..7 {
                let ip = b.column(i).dot(&b.column(j));
                assert!(ip.abs() < tol, "ip = {}, between index {}, and {}", ip, i, j);
            }
        }
    }
}
