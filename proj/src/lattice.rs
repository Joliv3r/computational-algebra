extern crate openblas_src;
use ndarray::{Array2, Array1};

pub mod methods;

#[derive(Debug)]
// Due to complications, this is implemented for only f64, mostly because of time constraints.
// TODO: Make Lattice<T> for some generic type.
pub struct Lattice {
    pub basis: Array2<f64>
}

impl Lattice {
    pub fn build_lattice_basis_from_vectors(vectors: &Vec<Array1<f64>>) -> Option<Lattice> {
        if let Ok(basis) = methods::make_into_basis_matrix(vectors) {
            Some(Lattice {
                basis,
            })
        } else {
            None
        }
    }

    pub fn get_length_of_basis_vectors(&self) -> usize {
        self.basis.shape()[0]
    }

    pub fn get_number_of_basis_vectors(&self) -> usize {
        self.basis.shape()[1]
    }

    pub fn get_basis_vector(&self, index: usize) -> Array1<f64> {
        self.basis.column(index).to_vec().into()
    }
}
