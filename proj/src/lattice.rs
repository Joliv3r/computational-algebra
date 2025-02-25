extern crate openblas_src;
use ndarray::{Array2, ArrayBase, Dimension};
use std::sync::Arc;

pub mod methods;

#[derive(Debug)]
pub struct LatticeBasis<T> {
    matrix: Array2<T>,
}

pub struct Lattice<T> {
    dimension: usize,
    basis: Arc<LatticeBasis<T>>,
}

impl<T: Clone> LatticeBasis<T> {
    pub fn build_lattice_basis_from_vectors(vectors: &Vec<Vec<T>>) -> Option<LatticeBasis<T>> {
        if let Ok(matrix) = methods::make_into_basis_matrix(vectors) {
            Some(LatticeBasis {
                matrix,
            })
        } else {
            None
        }
    }
}
