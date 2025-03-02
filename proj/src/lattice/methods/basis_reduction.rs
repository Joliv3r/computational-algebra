use ndarray_linalg::Scalar;

use crate::lattice::Lattice;
use std::cmp::max;

impl Lattice {
    // This procedure is directly copied from Galbraith
    //   https://www.math.auckland.ac.nz/~sgal018/crypto-book/ch17.pdf
    pub fn lll_reduction(&mut self, delta: f64) -> Result<(), String> {
        if delta >= 1. || delta <= 0.25 {
            return Err("Given delta is out of the allowed range (1/4, 1)".to_string())
        }
        let mut k = 1;
        while k < self.columns() {
            for j in (0..k).rev() {
                let b_k = self.get_basis_vector(k).expect("Should exist.");

                let mu_kj = self.get_mu(k, j);
                if mu_kj.abs() > 0.5 {
                    let b_j = self.get_basis_vector(j).expect("Should exist.");
                    let new_bk = &b_k - &b_j*(mu_kj.round());
                    if let Err(e) = self.update_basis_vector(k, &new_bk) {
                        return Err(e)
                    }
                }
            }
            if self.get_gram_schmidt_length(k) > (delta - self.get_mu(k, k-1).powi(2))*self.get_gram_schmidt_length(k-1) {
                k += 1;
            } else {
                if let Err(e) = self.swap_basis_vectors(k, k-1) {
                    return Err(e)
                }
                k = max(1, k-1);
            }
        }
        Ok(())
    }
}
