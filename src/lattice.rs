extern crate openblas_src;

use itertools::Itertools;
use methods::gram_schmidt;
use ndarray::Array1;
use plotters::prelude::*;

pub mod methods;

#[derive(Debug)]
// Due to complications, this is implemented for only f64, mostly because of time constraints.
// TODO: Make Lattice<T> for some generic type.
// We require basis to be of full rank.
pub struct Lattice {
    basis: Vec<Array1<f64>>,
    gram_schmidt_basis: Vec<Array1<f64>>,
}

impl Lattice {
    pub fn build_lattice_basis_from_vectors(basis: &Vec<Array1<f64>>) -> Option<Lattice> {
        if !methods::is_linearly_independent(&basis) {
            return None
        }
                
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

    pub fn get_lattice_point(&self, combination: &Vec<i64>) -> Option<Array1<f64>> {
        if combination.len() != self.columns() {
            return None
        }
        let mut vector = Array1::zeros(self.columns());
        for (x_i, b_i) in combination.iter().zip(&self.basis) {
            vector = vector + b_i*(*x_i as f64);
        }
        Some(vector)
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

    fn update_basis_vector(&mut self, index: usize, new_vector: &Array1<f64>) -> Result<(), String> {
        if !self.index_exists(index) {
            return Err("Index out of range.".to_string())
        }

        self.basis[index] = new_vector.clone();
        self.update_gram_schmidt_basis();
        Ok(())
    }

    // TODO: Take indexes and update only these to make fewer unnecessary compuations.
    fn update_gram_schmidt_basis(&mut self) {
        self.gram_schmidt_basis = gram_schmidt(&self.basis);
    }

    fn swap_basis_vectors(&mut self, i: usize, j: usize) -> Result<(), String> {
        if !self.index_exists(i) || !self.index_exists(j) {
            return Err("Index out of range".to_string())
        }
        self.basis.swap(i, j);
        self.update_gram_schmidt_basis();
        Ok(())
    }

    fn get_lattice_points_2d(&self, n: i64) -> Vec<(f64, f64)>{
        let mut points = Vec::with_capacity(n.pow(2) as usize);
        for i in -n..n {
            for j in -n..n {
                let point = self.get_lattice_point(&vec![i, j]).expect("Function is only run when lattice is 2d").to_vec();
                points.push((point[0], point[1]));
            }
        }
        points
    }

    pub fn print_lattice_around_point(&self, point: (f64, f64), plot_size: usize, linear_combinations: i64, extra_points: &Vec<(f64, f64)>) -> Result<(), Box<dyn std::error::Error>> {
        // if self.columns() != 2 {
        //     return Err("Not 2d lattice")
        // }

        let endpoints = plot_size as f64;
        let points = self.get_lattice_points_2d(linear_combinations);

        let root = SVGBackend::new("images/lattice.svg", (1000, 1000)).into_drawing_area();
        root.fill(&WHITE)?;

        let mut chart = ChartBuilder::on(&root)
            .caption("Lattice", ("computer-modern", 30).into_font())
            .margin(40)
            .x_label_area_size(30)
            .y_label_area_size(50)
            .build_cartesian_2d(-endpoints+point.0..endpoints+point.0, -endpoints+point.1..endpoints+point.1)?;

        chart.configure_mesh()
            .x_desc("x-axis")
            .x_label_style(("computer-modern", 12).into_font())
            .y_desc("y-axis")
            .y_label_style(("computer-modern", 12).into_font())
            .draw()?;

        chart
            .draw_series(PointSeries::of_element(
                // (-50..=50).map(|x| x as f32 / 50.0).map(|x| (x, x * x)),
                points,
                3,
                &RED,
                &|c, s, st| {
                    return EmptyElement::at(c)    // We want to construct a composed element on-the-fly
                    + Circle::new((0,0),s,st.filled()) // At this point, the new pixel coordinate is established
                    // + Text::new(format!("{:?}", c), (10, 0), ("sans-serif", 10).into_font());
                },
            ))?;
            // .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

        chart
            .draw_series(PointSeries::of_element(
                // (-50..=50).map(|x| x as f32 / 50.0).map(|x| (x, x * x)),
                extra_points.clone(),
                3,
                &BLUE,
                &|c, s, st| {
                    return EmptyElement::at(c)    // We want to construct a composed element on-the-fly
                    + Circle::new((0,0),s,st.filled()) // At this point, the new pixel coordinate is established
                    // + Text::new(format!("{:?}", c), (10, 0), ("sans-serif", 10).into_font());
                },
            ))?;
            // .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

        chart
            .draw_series(PointSeries::of_element(
                // (-50..=50).map(|x| x as f32 / 50.0).map(|x| (x, x * x)),
                vec![point],
                3,
                &BLACK,
                &|c, s, st| {
                    return EmptyElement::at(c)    // We want to construct a composed element on-the-fly
                    + Circle::new((0,0),s,st.filled()) // At this point, the new pixel coordinate is established
                    // + Text::new(format!("{:?}", c), (10, 0), ("sans-serif", 10).into_font());
                },
            ))?;
            // .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

        root.present()?;

        Ok(())
    }
}
