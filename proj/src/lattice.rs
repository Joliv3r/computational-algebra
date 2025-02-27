extern crate openblas_src;

use itertools::Itertools;
use methods::gram_schmidt;
use ndarray::Array1;

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

    pub fn print_basis(&self) {
        for i in &self.basis {
            println!("{}", i);
        }
    }

    pub fn print_gram_schmidt_basis(&self) {
        for i in &self.gram_schmidt_basis {
            println!("{}", i);
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
        if !self.index_exists(end-1) {
            return None
        }
        Some(self.gram_schmidt_basis[start..end].to_vec())
    }

    fn index_exists(&self, index: usize) -> bool {
        if index >= self.columns() {
            return false
        }
        true
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

    pub fn update_basis_vector(&mut self, index: usize, new_vector: &Array1<f64>) -> Result<(), String> {
        if !self.index_exists(index) {
            return Err("Index out of range.".to_string())
        }

        self.basis[index] = new_vector.clone();
        Ok(())
    }

    // TODO: Take indexes and update only these to make fewer unnecessary compuations.
    fn update_gram_schmidt_basis(&mut self) {
        self.gram_schmidt_basis = gram_schmidt(&self.basis);
    }

    pub fn swap_basis_vectors(&mut self, i: usize, j: usize) -> Result<(), String> {
        if !self.index_exists(i) || !self.index_exists(j) {
            return Err("Index out of range".to_string())
        }
        self.basis.swap(i, j);
        self.update_gram_schmidt_basis();
        Ok(())
    }
}
