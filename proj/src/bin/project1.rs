use std::{sync::Arc, time::Duration};
use std::time::Instant;
use beralg::{algebraic_structure::{finite_field::FiniteField, Element}, integers::integer_computations::naive_pow};
use rug::integer::IsPrime;
use rug::ops::PowAssign;
use rug::Integer;
use plotters::prelude::*;
use plotters::coord::combinators::IntoLogRange;
use beralg::random::{randint_bits, randint_digits};


fn check_timing_against_rug(a: &Integer, b: &Integer, p: &Integer) -> (Duration, Duration) {
    let f = Arc::new(FiniteField::new(p.clone()).unwrap());
    let a_elem = Element::new(f.clone(), a.clone());
    let b_exp = b.clone();

    let now = Instant::now();
    a_elem.pow(&b_exp);
    let elapsed_elem = now.elapsed();

    let now = Instant::now();
    a.clone().pow_mod(&b, &p).unwrap();
    let elapsed = now.elapsed();

    (elapsed_elem, elapsed)
}


fn check_timing_naive(a: &Integer, b: &Integer, p: &Integer, n: usize) -> Duration {
    let now = Instant::now();
    for _ in 0..n {
        naive_pow(a, b, p);
    }
    now.elapsed()/50
}


fn check_timing_square(a: &Integer, b: &Integer, p: &Integer, n: usize) -> Duration {
    let f = Arc::new(FiniteField::new(p.clone()).unwrap());
    let a_elem = Element::new(f.clone(), a.clone());
    
    let now = Instant::now();
    for _ in 0..n {
        a_elem.pow(b);
    }
    now.elapsed()/50
}


fn check_timing_against_naive(a: &Integer, b: &Integer, p: &Integer, n: usize) -> (Duration, Duration) {
    (check_timing_naive(a, b, p, n), check_timing_square(a, b, p, n))
}


fn plot_timing_naive_square(n: usize, m: usize) -> Result<(), Box<dyn std::error::Error>> {
    let mut q = Integer::from(17);
    let mut p = Integer::from(2);
    p.pow_assign(127);
    p = p-1;
    assert!(p.is_probably_prime(20) != IsPrime::No);
    // let mut rng = RandState::new();
    let mut naive_vec: Vec<(u128, u128)> = Vec::new();
    let mut square_vec: Vec<(u128, u128)> = Vec::new();
    let mut max_time_naive = 0;
    // let mut last_bit_number = p.significant_bits();
    for _ in 0..n {
        // while p.significant_bits() == last_bit_number {
        //     p.next_prime_mut();
        // }
        // last_bit_number = p.significant_bits();
        for _ in 0..5 {
            q.next_prime_mut();
        }

        let a = randint_digits(q.to_string().len());
        // let a = p.clone() - Integer::ONE.clone();
        // let b = randint_digits(p.significant_digits::<usize>());
        let b = q.clone();
        let (elapsed_naive, elapsed_square) = check_timing_against_naive(&a, &b, &p, m);
        naive_vec.push((q.to_u128().unwrap(), elapsed_naive.as_nanos()));
        square_vec.push((q.to_u128().unwrap(), elapsed_square.as_nanos()));
        if elapsed_naive.as_nanos() > max_time_naive {
            max_time_naive = elapsed_naive.as_nanos();
        }
    }

    // for i in &square_vec {
    //     println!("{}", i.1);
    // }

    let root = SVGBackend::new("../latex/proj1/images/naive-square.svg", (600, 400)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("Runtime of Naive vs Square-Multiply", ("computer-modern", 30).into_font())
        .margin(40)
        .x_label_area_size(30)
        .y_label_area_size(50)
        .build_cartesian_2d(0..q.to_u128().unwrap_or(u128::MAX), 0..max_time_naive)?;

    chart.configure_mesh()
        .x_desc("Prime number")
        .x_label_style(("computer-modern", 12).into_font())
        .y_desc("Nanoseconds")
        .y_label_style(("computer-modern", 12).into_font())
        .draw()?;

    chart
        .draw_series(PointSeries::of_element(
            // (-50..=50).map(|x| x as f32 / 50.0).map(|x| (x, x * x)),
            naive_vec,
            3,
            &RED,
            &|c, s, st| {
                return EmptyElement::at(c)    // We want to construct a composed element on-the-fly
                + Circle::new((0,0),s,st.filled()) // At this point, the new pixel coordinate is established
                // + Text::new(format!("{:?}", c), (10, 0), ("sans-serif", 10).into_font());
            },
        ))?
        .label("Naive")
        .legend(|(x, y)| Circle::new((x, y), 3, RED.filled()));
        // .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    chart
        .draw_series(PointSeries::of_element(
            // (-50..=50).map(|x| x as f32 / 50.0).map(|x| (x, x * x)),
            square_vec,
            3,
            &BLUE,
            &|c, s, st| {
                return EmptyElement::at(c)    // We want to construct a composed element on-the-fly
                + Circle::new((0,0),s,st.filled()) // At this point, the new pixel coordinate is established
                // + Text::new(format!("{:?}", c), (10, 0), ("sans-serif", 10).into_font());
            },
        ))?
        .label("Square-Multiply")
        // .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));
        .legend(|(x, y)| Circle::new((x, y), 3, BLUE.filled()));

    chart
        .configure_series_labels()
        .label_font(("computer-modern", 12).into_font())
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .legend_area_size(12)
        .draw()?;

    root.present()?;


    Ok(())
}


fn plot_timing_naive(n: usize, m: usize) -> Result<(), Box<dyn std::error::Error>> {
    let mut p = Integer::ONE.clone();
    let mut q = Integer::from(2);
    q.pow_assign(127);
    q = q-1;
    assert!(q.is_probably_prime(20) != IsPrime::No);
    // let mut rng = RandState::new();
    let mut naive_vec: Vec<(u128, u128)> = Vec::new();
    let mut max_time_naive = 0;
    // let mut last_bit_number = p.significant_bits();
    for _ in 1..=n {
        // while p.significant_bits() == last_bit_number {
        //     p.next_prime_mut();
        // }
        // last_bit_number = p.significant_bits();
        for _ in 0..400 {
            p.next_prime_mut();
        }
        // p = randint_digits(i as usize).next_prime();

        let a = randint_digits(p.to_string().len());
        // let a = p.clone() - Integer::ONE.clone();
        // let b = randint_digits(p.to_string().len());
        let b = p.clone();
        println!("{}", b);
        let elapsed_naive = check_timing_naive(&a, &b, &q, m);
        naive_vec.push((p.to_u128().unwrap(), elapsed_naive.as_micros()));
        if elapsed_naive.as_micros() > max_time_naive {
            max_time_naive = elapsed_naive.as_micros();
        }
        // println!("{}, {}, {}", &p, elapsed_square.as_nanos(), elapsed_naive.as_nanos());
    }

    let root = SVGBackend::new("../latex/proj1/images/naive.svg", (600, 400)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("Runtime of Naive Exponentiation", ("computer-modern", 30).into_font())
        .margin(40)
        .x_label_area_size(30)
        .y_label_area_size(50)
        .build_cartesian_2d(0..p.to_u128().unwrap_or(u128::MAX), 0..max_time_naive)?;

    chart.configure_mesh()
        .x_desc("Prime number")
        .x_label_style(("computer-modern", 12).into_font())
        .y_desc("Microseconds")
        .y_label_style(("computer-modern", 12).into_font())
        .draw()?;

    chart
        .draw_series(PointSeries::of_element(
            // (-50..=50).map(|x| x as f32 / 50.0).map(|x| (x, x * x)),
            naive_vec,
            3,
            &RED,
            &|c, s, st| {
                return EmptyElement::at(c)    // We want to construct a composed element on-the-fly
                + Circle::new((0,0),s,st.filled()) // At this point, the new pixel coordinate is established
                // + Text::new(format!("{:?}", c), (10, 0), ("sans-serif", 10).into_font());
            },
        ))?;
        // .label("Naive approach")
        // .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));
        // .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    // chart
    //     .configure_series_labels()
    //     .background_style(&WHITE.mix(0.8))
    //     .border_style(&BLACK)
    //     .draw()?;

    root.present()?;
    
    Ok(())
}


fn plot_timing_square(n: usize, m: usize) -> Result<(), Box<dyn std::error::Error>> {

    let mut q = Integer::from(2);
    q.pow_assign(127);
    q = q-1;
    assert!(q.is_probably_prime(20) != IsPrime::No);
    let mut p = Integer::ONE.clone();
    // let mut rng = RandState::new();
    let mut square_vec: Vec<(u64, u64)> = Vec::new();
    let mut max_time_square = 0;
    // let mut last_bit_number = p.significant_bits();
    for i in 1..=n {
        // while p.significant_bits() == last_bit_number {
        //     p.next_prime_mut();
        // }
        // last_bit_number = p.significant_bits();
        // for _ in 0..400 {
        //     p.next_prime_mut();
        // }
        p = randint_bits(i as usize).next_prime();

        let a = randint_digits(p.to_string().len());
        // let b = randint_digits(p.significant_digits::<usize>());
        // let a = p.clone() - Integer::ONE.clone();
        let b = p.clone();
        let elapsed_square = check_timing_square(&a, &b, &q, m);
        square_vec.push((p.to_u64().unwrap(), elapsed_square.as_nanos() as u64));
        if elapsed_square.as_nanos() > max_time_square {
            max_time_square = elapsed_square.as_nanos();
        }
        // println!("{}, {}, {}", &p, elapsed_square.as_nanos(), elapsed_naive.as_nanos());
    }


    let root = SVGBackend::new("../latex/proj1/images/square.svg", (600, 400)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("Runtime of Square-Multiply", ("computer-modern", 30).into_font())
        .margin(30)
        .x_label_area_size(30)
        .y_label_area_size(50)
        .build_cartesian_2d((0..p.to_u64().unwrap_or(u64::MAX)).log_scale(), 0..max_time_square as u64)?;

    chart.configure_mesh()
        .x_desc("Prime number")
        .x_label_formatter(&|x| format!("10^{}", x.ilog10()))
        .x_label_style(("computer-modern", 12).into_font())
        .y_desc("Nanoseconds")
        .y_label_style(("computer-modern", 12).into_font())
        .draw()?;

    chart
        .draw_series(PointSeries::of_element(
            // (-50..=50).map(|x| x as f32 / 50.0).map(|x| (x, x * x)),
            square_vec,
            3,
            &BLUE,
            &|c, s, st| {
                return EmptyElement::at(c)    // We want to construct a composed element on-the-fly
                + Circle::new((0,0),s,st.filled()) // At this point, the new pixel coordinate is established
                // + Text::new(format!("{:?}", c), (10, 0), ("sans-serif", 10).into_font());
            },
        ))?;
        // .label("Square approach")
        // .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));

    // chart
    //     .configure_series_labels()
    //     .background_style(&WHITE.mix(0.8))
    //     .border_style(&BLACK)
    //     .draw()?;

    root.present()?;

    Ok(())
}


fn main() {
    let loops = 10;
    let naive_square_points: usize = 20;
    let naive_points: usize = 50;
    let square_points: usize = 50;
    plot_timing_naive_square(naive_square_points, loops).expect("Should not fail");
    plot_timing_naive(naive_points, loops).expect("Should not fail");
    plot_timing_square(square_points, loops).expect("Should not fail");
}
