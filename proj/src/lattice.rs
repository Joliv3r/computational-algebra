extern crate openblas_src;

use itertools::Itertools;
use ndarray::{Array1, Axis, Slice};

pub mod methods;

#[derive(Debug)]
// Due to complications, this is implemented for only f64, mostly because of time constraints.
// TODO: Make Lattice<T> for some generic type.
pub struct Lattice {
    basis: Vec<Array1<f64>>,
    gram_schmidt_basis: Vec<Array1<f64>>,
}

impl Lattice {
    pub fn build_lattice_basis_from_vectors(basis: &Vec<Array1<f64>>) -> Option<Lattice> {
        if basis.iter().map(|v| v.len()).all_equal() {
            Some(Lattice {
                basis: basis.clone(),
                gram_schmidt_basis: methods::gram_schmidt(basis)
            })
        } else {
            None
        }
    }

    pub fn get_basis_columns(&self, start: usize, end: usize) -> Option<Vec<Array1<f64>>> {
        // self.basis.slice_axis(Axis(1), slice)
        if self.columns() < end {
            return None
        }
        Some(self.basis[start..end].to_vec())
    }

    // Should only be called if Gram-Schmidt basis is already calculated
    pub fn get_gram_schmidt_basis_columns(&self, start: usize, end: usize) -> Option<Vec<Array1<f64>>> {
        // self.gram_schmidt_basis.slice_axis(Axis(1), slice)
        if self.columns() < end {
            return None
        }
        Some(self.gram_schmidt_basis[start..end].to_vec())
    }

    pub fn get_length_of_basis_vectors(&self) -> usize {
        self.basis[0].len()
    }

    pub fn columns(&self) -> usize {
        self.basis.len()
    }

    pub fn get_basis_vector(&self, index: usize) -> Option<Array1<f64>> {
        self.basis.get(index).cloned()
    }

    pub fn get_gram_schmidt_basis_vector(&self, index: usize) -> Option<Array1<f64>> {
        self.gram_schmidt_basis.get(index).cloned()
    }

    pub fn get_shortest_basis_vector(&self) -> Option<Array1<f64>> {
        if self.columns() == 0 {
            return None
        }
        let mut shortest_basis_vector = self.get_basis_vector(0).expect("Should exist.");
        let length = shortest_basis_vector.dot(&shortest_basis_vector);
        for column in self.get_basis_columns(1, self.columns()).expect("Should exist.") {
            let vector_length = column.dot(&column);
            if vector_length < length {
                shortest_basis_vector = column;
            }
        }
        Some(shortest_basis_vector)
    }
}
