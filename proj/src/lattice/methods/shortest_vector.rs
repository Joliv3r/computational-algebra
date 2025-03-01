use crate::lattice::Lattice;
use ndarray::Array1;
use itertools::Itertools;

impl Lattice {
    // Say v = sum_{i=1}^{n} x_j b_j, is the shortest vector. We use that
    //   -(M_1 + M_2) =< x_i => M_1 - M_2
    // where
    //   M_1 = sqrt{ ( A - sum_{j = i+1}^{n} x_j^2 B_j )/B_i }
    //   M_2 = sum_{j = i+1}^{n} µ_{j,i} x_j
    // with A > ||v||^2 and B_j = ||b_j||^2 and µ_{j,i} = <b_i, b*_j>/||b*_j||^2.
    pub fn shortest_vector_by_enumeration(&self) -> Result<(Array1<f64>, Vec<(f64, f64)>), String> {
        if self.columns() == 0 {
            return Err("Lattice is empty.".to_string())
        }
        let mut shortest_vector: Array1<f64> = self.get_basis_vector(self.columns()-1).expect("Should exist.").to_vec().into();
        let mut shortest_length = shortest_vector.dot(&shortest_vector);
        let mut checked_points: Vec<(f64, f64)> = vec![(shortest_vector[0], shortest_vector[1])];

        let mut combination = vec![0; self.columns()];
        let bound_n = (shortest_length/self.get_gram_schmidt_length(self.columns()-1)).sqrt().floor() as i64;

        println!("Bound on n is {}", bound_n);

        for x_n in 0..bound_n {
            combination[self.columns()-1] = x_n;
            let (candidate_vector, candidate_length) = self.shortest_vector_enumeration_steps(1, &mut combination, &shortest_vector.to_vec().into(), shortest_length, &mut checked_points);
            if candidate_length < shortest_length && candidate_length != 0. {
                shortest_vector = candidate_vector;
                shortest_length = candidate_length;
            }
        }

        Ok((shortest_vector, checked_points))
    }

    fn shortest_vector_enumeration_steps(&self, depth: usize, combination: &mut Vec<i64>, current_shortes_vector: &Array1<f64>, current_shortest_length: f64, checked_points: &mut Vec<(f64, f64)>) -> (Array1<f64>, f64) {
        if depth == self.columns() {
            let current_vector = self.get_lattice_point(combination).expect("Should exist.");
            let current_length = current_vector.dot(&current_vector);
            println!("Depth reached.");
            checked_points.push((current_vector[0], current_vector[1]));
            return (current_vector, current_length)
        }

        let mut shortest_vector = current_shortes_vector.clone();
        let mut shortest_length = current_shortest_length;
        let (lower_bound, upper_bound) = self.get_svp_enumeration_bounds(combination, depth, shortest_length);
        println!("Checking in {}..{}", lower_bound, upper_bound);
        for i in lower_bound..=upper_bound {
            combination[self.columns()-1-depth] = i;
            let (candidate_vector, candidate_length) = self.shortest_vector_enumeration_steps(depth+1, combination, &shortest_vector, shortest_length, checked_points);
            if candidate_length < shortest_length && candidate_length != 0. {
                shortest_vector = candidate_vector;
                shortest_length = candidate_length;
                println!("New shortest distance found: {}, with {}", shortest_length, shortest_vector);
            }
        }

        (shortest_vector, shortest_length)
    }

    #[allow(non_snake_case)]
    fn get_svp_enumeration_bounds(&self, combination: &Vec<i64>, basis_number: usize, A: f64) -> (i64, i64) {
        let mut sum: f64 = 0.;
        let mut M_2: f64 = 0.;
        let start = self.columns()-basis_number+1;
        let end = self.columns();
        let i = self.columns()-basis_number;
        for j in (start..end).into_iter().rev() {
            let B_j = self.get_gram_schmidt_length(j);
            let mu_ji = self.get_mu(j, i);
            sum += (combination[j].pow(2) as f64)*B_j;
            M_2 += mu_ji*(combination[j] as f64);
        }
        let M_1 = ((A-sum)/self.get_gram_schmidt_length(i)).sqrt();

        (-(M_1+M_2).ceil() as i64, (M_1-M_2).floor() as i64)
    }
}
