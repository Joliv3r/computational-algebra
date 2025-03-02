use ndarray::Array1;
use rand::{thread_rng, Rng};
use crate::lattice::{methods::{get_length_of_vector, shortest_vector}, Lattice};

use super::is_linearly_independent;

pub fn generate_random_basis(dimension: usize) -> Vec<Array1<f64>> {
    let mut basis: Vec<Array1<f64>> = Vec::with_capacity(dimension);
    for _ in 0..dimension {
        basis.push(generate_random_vector(dimension, 1.));
    }

    while !is_linearly_independent(&basis) {
        for i in 0..dimension {
            basis[i] = generate_random_vector(dimension, 1.);
        }
    }
    basis
}

pub fn generate_random_vector(dimension: usize, scaling: f64) -> Array1<f64> {
    let mut rng = thread_rng();
    let mut vector: Array1<f64> = Array1::zeros(dimension);
    vector.map_inplace(|e| {*e = scaling*(rng.gen_range(-1000..1000) as f64)});
    vector
}


pub fn cvp_statistics(top: usize) {
    for dimension in 2..top {
        println!("For dimension {}, we found the following:", dimension);
        
        let scale = 20.;
        let basis = generate_random_basis(dimension);
        let vector = generate_random_vector(dimension, scale);
        let mut lattice = Lattice::build_lattice_basis_from_vectors(&basis).expect("Basis is square.");
        
        let babai_pre_lll = lattice.babai_nearest_plane(&vector).expect("Should be well-defined.");
        let babai_pre_lll_distance = get_length_of_vector(&(&vector - &babai_pre_lll));

        let delta = 0.75;
        lattice.lll_reduction(delta);
        let babai_post_lll = lattice.babai_nearest_plane(&vector).expect("Should be well-defined.");
        let babai_post_lll_distance = get_length_of_vector(&(&vector - &babai_post_lll));

        let closest_vector = lattice.closest_vector_by_enumeration(&vector).expect("Should be well-defined.");
        let closest_vector_distance = get_length_of_vector(&(&vector - &closest_vector));

        let shortest_vector = lattice.shortest_vector_by_enumeration().expect("Should be well-defined.");
        let shortest_vector_distance = get_length_of_vector(&shortest_vector);

        let shortest_basis_vector = lattice.get_basis_vector(0).expect("Should exist.");
        let shortest_basis_vector_distance = get_length_of_vector(&shortest_basis_vector);

        println!("distance for babai_pre_lll:          {}", babai_pre_lll_distance);
        println!("distance for babai_post_lll:         {}", babai_post_lll_distance);
        println!("distance for closest_vector:         {}", closest_vector_distance);
        println!("distance for shortest_vector:        {}", shortest_vector_distance);
        println!("distance for shortest_basis_vector:  {}", shortest_basis_vector_distance);
        println!();
    }
}
