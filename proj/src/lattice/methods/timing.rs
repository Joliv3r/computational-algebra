use ndarray::Array1;
use rand::{thread_rng, Rng};
use crate::lattice::Lattice;

use super::is_linearly_independent;

pub fn generate_random_basis(dimension: usize) -> Vec<Array1<f64>> {
    let mut basis: Vec<Array1<f64>> = Vec::with_capacity(dimension);
    for _ in 0..dimension {
        basis.push(generate_random_vector(dimension));
    }

    while !is_linearly_independent(&basis) {
        for i in 0..dimension {
            basis[i] = generate_random_vector(dimension);
        }
    }
    basis
}

pub fn generate_random_vector(dimension: usize) -> Array1<f64> {
    let mut rng = thread_rng();
    let mut vector: Array1<f64> = Array1::zeros(dimension);
    vector.map_inplace(|e| {*e = rng.gen_range(0..1000) as f64/(rng.gen_range(1..100) as f64)});
    vector
}


pub fn cvp_statistics(top: usize) {
    for dimension in 2..top {
        println!("For dimension {}, we found the following:", dimension);
        
        let basis = generate_random_basis(dimension);
        let vector = generate_random_vector(dimension);
        let mut lattice = Lattice::build_lattice_basis_from_vectors(&basis).expect("Basis is square.");
        
        let babai_pre_lll = lattice.babai_nearest_plane(&vector).expect("Should be well-defined.");

        let delta = 0.75;
        lattice.lll_reduction(delta);
        let babai_post_lll = lattice.babai_nearest_plane(&vector).expect("Should be well-defined.");

        let (closest_vector, _) = lattice.closest_vector_by_enumeration(&vector).expect("Should be well-defined.");

        let (shortest_vector, _) = lattice.shortest_vector_by_enumeration().expect("Should be well-defined.");

        println!("babai_pre_lll:    {}", babai_pre_lll);
        println!("babai_post_lll:   {}", babai_post_lll);
        println!("closest_vector:   {}", closest_vector);
        println!("shortest_vector:  {}", shortest_vector);
    }
}
