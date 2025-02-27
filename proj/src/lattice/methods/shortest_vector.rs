use crate::lattice::Lattice;
use ndarray::Array1;

impl Lattice {
    // Say v = sum_{i=1}^{n} x_j b_j, is the shortest vector. We use that
    //   -(M_1 + M_2) =< x_i => M_1 - M_2
    // where
    //   M_1 = sqrt{ ( A - sum_{j = i+1}^{n} x_j^2 B_j )/B_i }
    //   M_2 = sum_{j = i+1}^{n} µ_{j,i} x_j
    // with A > ||v||^2 and B_j = ||b_j||^2 and µ_{j,i} = <b_i, b*_j>/||b*_j||^2.
    pub fn shortest_vector_by_enumeration(&self) -> Result<Array1<f64>, String> {
        if self.columns() == 0 {
            return Err("Lattice is empty.".to_string())
        }
        let mut shortest_vector: Array1<f64> = self.get_shortest_basis_vector().expect("Should exist.").to_vec().into();
        let mut shortest_length = shortest_vector.dot(&shortest_vector);

        let b_orth_n = self.get_gram_schmidt_basis_vector(self.columns()-1).expect("Should exist.");
        let bound_n = (shortest_length/b_orth_n.dot(&b_orth_n)).sqrt().floor() as i32;



        for x_n in 0..bound_n {
            let combination = vec![x_n];
            let (candidate_vector, candidate_length) = self.shortest_vector_enumeration_steps(1, &combination, &shortest_vector.to_vec().into(), shortest_length);
            if candidate_length < shortest_length && candidate_length != 0. {
                shortest_vector = candidate_vector;
                shortest_length = candidate_length;
            }
        }

        Ok(shortest_vector)
    }

    fn shortest_vector_enumeration_steps(&self, depth: usize, combination: &Vec<i32>, current_shortes_vector: &Array1<f64>, current_shortest_length: f64) -> (Array1<f64>, f64) {
        if depth == self.columns() {
            let mut current_vector = Array1::zeros(self.get_length_of_basis_vectors());
            for (x_i, b_i) in combination.iter().zip(self.get_basis_columns(0, self.columns()).expect("Should exist").into_iter().rev()) {
                current_vector = current_vector + b_i * (*x_i as f64);
            }
            let current_length = current_vector.dot(&current_vector);
            return (current_vector, current_length)
        }

        let mut shortest_vector = current_shortes_vector.clone();
        let mut shortest_length = current_shortest_length.to_owned();
        let (lower_bound, upper_bound) = self.get_enumeration_bounds(combination, depth, shortest_length);
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
        let start = self.columns()-basis_number+1;
        let end = self.columns();
        let i = self.columns()-basis_number;
        for (x_j, j) in combination.iter().zip((start..end).into_iter().rev()) {
            let B_j = self.get_gram_schmidt_length(j);
            let mu_ij = self.get_mu(i, j);
            sum += (x_j.pow(2) as f64)*B_j;
            M_2 += mu_ij*(*x_j as f64);
        }
        let M_1 = ((A-sum)/self.get_gram_schmidt_length(i)).sqrt();

        (-(M_1+M_2).ceil() as i32, (M_1-M_2).floor() as i32)
    }
}
