use crate::lattice::Lattice;
use ndarray::Array1;

impl Lattice {
    pub fn babai_nearest_plane(&self, vector: &Array1<f64>) -> Result<Array1<f64>, String> {
        let rows = self.get_length_of_basis_vectors();
        if rows != vector.len() {
            return Err("Vector size not compatible with lattice".to_string())
        }

        let dim = self.columns();
        if dim == 0 {
            return Err("Lattice is empty.".to_string())
        }
        let mut w = vector.clone();
        let mut y = Array1::zeros(rows);

        for i in (0..dim).into_iter().rev() {
            let gs_i = self.get_gram_schmidt_basis_vector(i).expect("Vector should exist.");
            let b_i = self.get_basis_vector(i).expect("Vector should exist.");
            let l_i = &w.dot(&gs_i)/gs_i.dot(&gs_i);
            let l_i_rounded = l_i.round();

            y = y + &b_i*l_i_rounded;

            w = w - gs_i*(l_i - l_i_rounded) - b_i*l_i_rounded;
        }

        Ok(y)
    }
}
