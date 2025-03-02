use std::time::{Instant, SystemTime};

use plotters::{coord::combinators::IntoLogRange, style::full_palette::{BLUE, GREEN, ORANGE}};
use itertools::{max, min};
use plotters::prelude::*;
use ndarray::{array, Array1, Axis};
use ndarray_linalg::generate;
use rand::{thread_rng, Rng};
use rug::az::CastFrom;
use crate::lattice::{methods::{get_length_of_vector, shortest_vector}, Lattice};

use super::is_linearly_independent;

pub fn generate_random_basis(dimension: usize) -> Vec<Array1<f64>> {
    let mut basis: Vec<Array1<f64>> = Vec::with_capacity(dimension);
    for _ in 0..dimension {
        basis.push(generate_random_vector(dimension, 1.));
    }

    while !is_linearly_independent(&basis) {
        for i in 0..dimension {
            basis[i] = generate_random_vector(dimension, 1.);
        }
    }
    basis
}

pub fn generate_random_vector(dimension: usize, scaling: f64) -> Array1<f64> {
    let mut rng = thread_rng();
    let mut vector: Array1<f64> = Array1::zeros(dimension);
    vector.map_inplace(|e| {*e = scaling*(rng.gen_range(-1000..1000) as f64)});
    vector
}


pub fn cvp_statistics(top: usize) {
    for dimension in 2..top {
        println!("For dimension {}, we found the following:", dimension);
        
        let scale = 20.;
        let basis = generate_random_basis(dimension);
        let vector = generate_random_vector(dimension, scale);
        let mut lattice = Lattice::build_lattice_basis_from_vectors(&basis).expect("Basis is square.");
        
        let babai_pre_lll = lattice.babai_nearest_plane(&vector).expect("Should be well-defined.");
        let babai_pre_lll_distance = get_length_of_vector(&(&vector - &babai_pre_lll));

        let delta = 0.75;
        lattice.lll_reduction(delta);
        let babai_post_lll = lattice.babai_nearest_plane(&vector).expect("Should be well-defined.");
        let babai_post_lll_distance = get_length_of_vector(&(&vector - &babai_post_lll));

        let closest_vector = lattice.closest_vector_by_enumeration(&vector).expect("Should be well-defined.");
        let closest_vector_distance = get_length_of_vector(&(&vector - &closest_vector));

        let shortest_vector = lattice.shortest_vector_by_enumeration().expect("Should be well-defined.");
        let shortest_vector_distance = get_length_of_vector(&shortest_vector);

        let shortest_basis_vector = lattice.get_basis_vector(0).expect("Should exist.");
        let shortest_basis_vector_distance = get_length_of_vector(&shortest_basis_vector);

        println!("distance for babai_pre_lll:          {}", babai_pre_lll_distance);
        println!("distance for babai_post_lll:         {}", babai_post_lll_distance);
        println!("distance for closest_vector:         {}", closest_vector_distance);
        println!("distance for shortest_vector:        {}", shortest_vector_distance);
        println!("distance for shortest_basis_vector:  {}", shortest_basis_vector_distance);
        println!();
    }
}


pub fn increase_basis(basis: &mut Vec<Array1<f64>>) {
    let mut rng = thread_rng();
    for vector in basis.iter_mut() {
        let append = rng.gen_range(-1000..1000) as f64;
        vector.append(Axis(0), array![append].view());
    }
    let dimension = basis[0].len();
    basis.push(generate_random_vector(dimension, 1.));
    while !is_linearly_independent(&basis) {
        basis[dimension-1] = generate_random_vector(dimension, 1.);
    }
}


pub fn babai_times(top: usize) {
    let mut times_babai_pre_lll = Vec::with_capacity(top-2);
    let mut times_babai_post_lll = Vec::with_capacity(top-2);
    let mut basis = generate_random_basis(2);
    let mut max_y = 0;

    for dimension in 2..top {
        let vector = generate_random_vector(dimension, 20.);
        let mut lattice = Lattice::build_lattice_basis_from_vectors(&basis).expect("Should be a square matrix.");

        let mut loop_num = 10;
        let now = Instant::now();
        for _ in 0..loop_num {
            lattice.babai_nearest_plane(&vector);
        }
        let elapsed = (now.elapsed()/loop_num as u32).as_micros() as u64;
        if elapsed > max_y {
            max_y = elapsed;
        }
        times_babai_pre_lll.push((dimension as u64, elapsed));

        lattice.lll_reduction(0.75);
        let now = Instant::now();
        for _ in 0..loop_num {
            lattice.babai_nearest_plane(&vector);
        }
        let elapsed = (now.elapsed()/loop_num as u32).as_micros() as u64;
        if elapsed > max_y {
            max_y = elapsed;
        }
        times_babai_post_lll.push((dimension as u64, elapsed));

        increase_basis(&mut basis);
    }

    plot_time(vec![(times_babai_pre_lll, "Before LLL-reduction"), (times_babai_post_lll, "After LLL-reduction")], "Babai's Nearest Plane Method", max_y, top as u64, "babai-times");
}


pub fn enumeration_times(top: usize) {
    let mut times_svp = Vec::with_capacity(top-2);
    let mut times_cvp = Vec::with_capacity(top-2);
    let mut basis = generate_random_basis(2);
    let mut max_y = 0;
    for dimension in 2..top {
        let vector = generate_random_vector(dimension, 20.);
        let index = dimension - 2;
        let loops = vec![15, 15, 15, 15, 10, 10, 10, 10, 10, 5, 5, 5, 5, 3, 3, 3, 3, 2, 2];
        let mut lattice = Lattice::build_lattice_basis_from_vectors(&basis).expect("Should be a square matrix.");
        lattice.lll_reduction(0.75);

        let mut loop_num = 1;
        if loops.len() > index {
            loop_num = loops[index];
        }

        let now = Instant::now();
        for _ in 0..loop_num {
            lattice.shortest_vector_by_enumeration();
        }
        let elapsed = (now.elapsed()/loop_num as u32).as_micros() as u64;
        if elapsed > max_y {
            max_y = elapsed;
        }
        times_svp.push((dimension as u64, elapsed));

        let now = Instant::now();
        for _ in 0..loop_num {
            lattice.closest_vector_by_enumeration(&vector);
        }
        let elapsed = (now.elapsed()/loop_num as u32).as_micros() as u64;
        if elapsed > max_y {
            max_y = elapsed;
        }
        times_cvp.push((dimension as u64, elapsed));

        increase_basis(&mut basis);
    }

    let times = vec![(times_svp, "SVP Enumeration"), (times_cvp, "CVP Enumeration")];

    plot_time(times, "Enumeration Methods", max_y, top as u64, "enumeration-methods");
}


fn plot_time(times: Vec<(Vec<(u64, u64)>, &str)>, caption: &str, max_y: u64, max_x: u64, filename: &str) -> Result<(), Box<dyn std::error::Error>> {
    let file: String = "images/".to_string() + filename + ".svg";

    let root = SVGBackend::new(&file, (1200, 800)).into_drawing_area();
    root.fill(&WHITE)?;

    let palette = [
        BLACK,
        ORANGE,
        BLUE,
        GREEN,
    ];

    // let y = (1..max_y).log_scale();
    let y = 1..max_y;
    let x = 1..max_x;

    let mut chart = ChartBuilder::on(&root)
        .caption(caption.to_string(), ("computer-modern", 50).into_font())
        .margin(40)
        .x_label_area_size(50)
        .y_label_area_size(70)
        .build_cartesian_2d(x, y)?;
        // .build_cartesian_2d(1..top, (min_y..max_y).log_scale())?;

    chart.configure_mesh()
        .x_desc("Dimension of lattice")
        .x_label_style(("computer-modern", 24).into_font())
        .y_desc("Microseconds")
        // .y_label_formatter(&|y| format!("10^{}", y.ilog10()))
        .y_label_style(("computer-modern", 24).into_font())
        .draw()?;

    for (i, (t, label)) in times.iter().enumerate() {
        let color = palette.get(i).expect("Should exist.");
        chart
            .draw_series(LineSeries::new(
                // (-50..=50).map(|x| x as f32 / 50.0).map(|x| (x, x * x)),
                t.clone(),
                &color,
            ))?;
        chart
            .draw_series(PointSeries::of_element(
                // (-50..=50).map(|x| x as f32 / 50.0).map(|x| (x, x * x)),
                t.clone(),
                3,
                &color,
                &|c, s, st| {
                    return EmptyElement::at(c)    // We want to construct a composed element on-the-fly
                    + Circle::new((0,0),s,st.filled()) // At this point, the new pixel coordinate is established
                    // + Text::new(format!("{:?}", c), (10, 0), ("sans-serif", 10).into_font());
                },
            ))?
            .label(label.to_string())
            // .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));
            .legend(move |(x, y)| Circle::new((x, y), 3, color.filled()));
    }

    chart
        .configure_series_labels()
        .position(SeriesLabelPosition::UpperLeft)
        .label_font(("computer-modern", 24).into_font())
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .legend_area_size(12)
        .draw()?;

    root.present()?;

    Ok(())
}


pub fn check_closest_vector_by_enumeration_limit() {
    let now = Instant::now();
    let mut dimension = 30;
    let mut basis = generate_random_basis(dimension);     

    loop {
        let vector = generate_random_vector(dimension, 20.);
        let mut lattice = Lattice::build_lattice_basis_from_vectors(&basis).expect("Should be a square matrix.");
        lattice.lll_reduction(0.75);
        
        let cvp = match lattice.closest_vector_by_enumeration(&vector) {
            Err(err) => println!("Error in dimension {}, Err: {}", dimension, err),
            Ok(_) => println!("Finished cvp for dimension {}, after {} seconds", dimension, now.elapsed().as_secs()),
        };

        dimension += 1;        
        increase_basis(&mut basis);
    }
}
