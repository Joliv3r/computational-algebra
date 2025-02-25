extern crate openblas_src;
use methods::gram_schmidt;
use ndarray::{Array2, Array1};

pub mod methods;

#[derive(Debug)]
// Due to complications, this is implemented for only f64, mostly because of time constraints.
// TODO: Make Lattice<T> for some generic type.
pub struct Lattice {
    pub basis: Array2<f64>,
    pub gram_schmidt_basis: Option<Array2<f64>>,
}

impl Lattice {
    pub fn build_lattice_basis_from_vectors(vectors: &Vec<Array1<f64>>) -> Option<Lattice> {
        if let Ok(basis) = methods::make_into_basis_matrix(vectors) {
            Some(Lattice {
                basis,
                gram_schmidt_basis: None
            })
        } else {
            None
        }
    }

    pub fn build_from_basis(basis: Array2<f64>) -> Lattice {
        Lattice { 
            basis,
            gram_schmidt_basis: None
        }
    }

    fn calculate_gram_schmidt_basis(&mut self) {
        if self.gram_schmidt_basis.is_none() {
            self.gram_schmidt_basis = Some(methods::gram_schmidt_columns(&self.basis));
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

    pub fn get_shortest_basis_vector_length(&self) -> f64 {
        let first_vec = self.basis.column(0);
        let mut length = first_vec.dot(&first_vec);
        for column in self.basis.columns() {
            let vector_length = column.dot(&column);
            if vector_length < length {
                length = vector_length;
            }
        }
        length
    }
}
