extern crate openblas_src;
use ndarray::{Array1, Array2, ArrayView1, ArrayView2, Axis, Slice};

pub mod methods;

#[derive(Debug)]
// Due to complications, this is implemented for only f64, mostly because of time constraints.
// TODO: Make Lattice<T> for some generic type.
// NOTE: When first implementing this struct my idea was to represent it by a matrix, so you could
//       do matrix operations. This however was never done, and the representation of a matrix is
//       not really useful. However too much code have been written for me to do the big rewrite
//       yet.
// TODO: Write basis to Vec<Array1<_>> instead of Array2<_>.
pub struct Lattice {
    basis: Array2<f64>,
    gram_schmidt_basis: Array2<f64>,
}

impl Lattice {
    pub fn build_lattice_basis_from_vectors(vectors: &Vec<Array1<f64>>) -> Option<Lattice> {
        if let Ok(basis) = methods::make_into_basis_matrix(vectors) {
            Some(Lattice {
                basis,
                gram_schmidt_basis: methods::make_matrix_from_column_vectors(&methods::gram_schmidt(vectors))
            })
        } else {
            None
        }
    }

    pub fn build_from_basis(basis: &Array2<f64>) -> Lattice {
        Lattice { 
            basis: basis.clone(),
            gram_schmidt_basis: methods::gram_schmidt_columns(basis).to_owned()
        }
    }

    pub fn get_basis_columns(&self, slice: Slice) -> ArrayView2<f64> {
        self.basis.slice_axis(Axis(1), slice)
    }

    // Should only be called if Gram-Schmidt basis is already calculated
    pub fn get_gram_schmidt_basis_columns(&self, slice: Slice) -> ArrayView2<f64> {
        self.gram_schmidt_basis.slice_axis(Axis(1), slice)
    }

    pub fn get_length_of_basis_vectors(&self) -> usize {
        self.basis.shape()[0]
    }

    pub fn columns(&self) -> usize {
        self.basis.shape()[1]
    }

    pub fn get_basis_vector(&self, index: usize) -> ArrayView1<f64> {
        self.basis.column(index)
    }

    pub fn get_gram_schmidt_basis_vector(&self, index: usize) -> ArrayView1<f64> {
        self.gram_schmidt_basis.column(index)
    }

    pub fn get_shortest_basis_vector(&self) -> ArrayView1<f64> {
        let first_vec = self.basis.column(0);
        let mut shortest_basis_vector = first_vec;
        let length = first_vec.dot(&first_vec);
        for column in self.basis.columns() {
            let vector_length = column.dot(&column);
            if vector_length < length {
                shortest_basis_vector = column;
            }
        }
        shortest_basis_vector
    }
}
