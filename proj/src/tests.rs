use std::{sync::Arc, time::Duration};
use std::time::Instant;
use crate::{algebraic_structure::{finite_field::FiniteField, Element}, integer_computations::naive_pow};
use rug::Integer;
use plotters::prelude::*;
use plotters::coord::combinators::IntoLogRange;
use rand::{thread_rng, Rng};


pub fn randint_bits(bits: usize) -> Integer {
    let mut n = Integer::from(1);
    let mut rng = thread_rng();
    for _ in 1..bits {
        n <<= 1;
        n += rng.gen_range(0..=1);
    }
    n
}


pub fn randint_digits(digits: usize) -> Integer {
    let mut rng = thread_rng();
    let mut n: Integer = Integer::from(rng.gen_range(1..=9));

    for _ in 0..digits {
        n *= 10;
        n += rng.gen_range(0..=9);
    }
    n
}


pub fn check_timing_against_rug(a: &Integer, b: &Integer, p: &Integer) -> (Duration, Duration) {
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


pub fn check_timing_against_naive(a: &Integer, b: &Integer, p: &Integer, n: usize) -> (Duration, Duration) {
    (check_timing_naive(a, b, p, n), check_timing_square(a, b, p, n))
}


pub fn plot_timing_naive_square(n: usize, m: usize) -> Result<(), Box<dyn std::error::Error>> {
    let mut p = Integer::from(17);
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
        for _ in 0..400 {
            p.next_prime_mut();
        }

        // let a = randint_digits(p.significant_digits::<usize>());
        let a = p.clone() - Integer::ONE.clone();
        // let b = randint_digits(p.significant_digits::<usize>());
        // let a = p.clone();
        let b = p.clone();
        let (elapsed_naive, elapsed_square) = check_timing_against_naive(&a, &b, &p, m);
        naive_vec.push((p.to_u128().unwrap(), elapsed_naive.as_nanos()));
        square_vec.push((p.to_u128().unwrap(), elapsed_square.as_nanos()));
        if elapsed_naive.as_nanos() > max_time_naive {
            max_time_naive = elapsed_naive.as_nanos();
        }
    }

    let root = SVGBackend::new("plots/naive-square.svg", (600, 400)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("Runtime in nanoseconds", ("sans-serif", 25).into_font())
        .margin(40)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(0..p.to_u128().unwrap_or(u128::MAX), 0..max_time_naive)?;

    chart.configure_mesh().draw()?;

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
        .label("Naive approach")
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
        .label("Square approach")
        // .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));
        .legend(|(x, y)| Circle::new((x, y), 3, BLUE.filled()));

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    root.present()?;

    Ok(())
}


pub fn plot_timing_naive(n: usize, m: usize) -> Result<(), Box<dyn std::error::Error>> {
    let mut p = Integer::ONE.clone();
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

        // let a = randint_digits(p.significant_digits::<usize>());
        // let b = randint_digits(p.significant_digits::<usize>());
        let a = p.clone() - Integer::ONE.clone();
        let b = p.clone();
        let elapsed_naive = check_timing_naive(&a, &b, &p, m);
        naive_vec.push((p.to_u128().unwrap(), elapsed_naive.as_nanos()));
        if elapsed_naive.as_nanos() > max_time_naive {
            max_time_naive = elapsed_naive.as_nanos();
        }
        // println!("{}, {}, {}", &p, elapsed_square.as_nanos(), elapsed_naive.as_nanos());
    }

    let root = SVGBackend::new("plots/naive.svg", (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("Runtime in nanoseconds", ("sans-serif", 50).into_font())
        .margin(40)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(0..p.to_u128().unwrap_or(u128::MAX), 0..max_time_naive)?;

    chart.configure_mesh().draw()?;

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
        .label("Naive approach")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));
        // .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    root.present()?;
    
    Ok(())
}


pub fn plot_timing_square(n: usize, m: usize) -> Result<(), Box<dyn std::error::Error>> {

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

        // let a = randint_digits(p.significant_digits::<usize>());
        // let b = randint_digits(p.significant_digits::<usize>());
        let a = p.clone() - Integer::ONE.clone();
        let b = p.clone();
        let elapsed_square = check_timing_square(&a, &b, &p, m);
        square_vec.push((p.to_u64().unwrap(), elapsed_square.as_nanos() as u64));
        if elapsed_square.as_nanos() > max_time_square {
            max_time_square = elapsed_square.as_nanos();
        }
        // println!("{}, {}, {}", &p, elapsed_square.as_nanos(), elapsed_naive.as_nanos());
    }


    let root = SVGBackend::new("plots/square.svg", (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("Runtime in nanoseconds", ("sans-serif", 50).into_font())
        .margin(40)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d((0..p.to_u64().unwrap_or(u64::MAX)).log_scale(), 0..max_time_square as u64)?;

    chart.configure_mesh().draw()?;

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
        .label("Square approach")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    root.present()?;

    Ok(())
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebraic_structure::Element;
    use crate::algebraic_structure::finite_field::MultiplicativeGroup;
    use rug::{Complete, rand::RandState};
    use std::sync::Arc;

    #[test]
    fn test_algebraic_structure_arithmetic() {
        let n: u32 = 200;
        let mut prime: Integer = Integer::from(2);
        let mut rng = RandState::new();

        for _ in 2..n {
            let p = prime.clone();
            let f: Arc<FiniteField> = Arc::new(FiniteField::new(p.clone()).unwrap());
            let a_rand: Integer = Integer::from(rng.bits(32));
            let b_rand: Integer = Integer::from(rng.bits(32));
            let a = Element::new(f.clone(), a_rand.clone());
            let b = Element::new(f.clone(), b_rand.clone());

            assert_eq!(a.get_rep(), &(&a_rand % &p).complete());
            assert_eq!(b.get_rep(), &(&b_rand % &p).complete());

            let added = (a.add_ref(&b)).get_rep().clone();
            let subtracted = (a.sub_ref(&b)).get_rep().clone();
            let multiplied = (a.mul_ref(&b)).get_rep().clone();

            let added_check: Integer = ((&a_rand % &p).complete() + (&b_rand % &p).complete()) % &p;
            let sub_check: Integer = ((&a_rand % &p).complete() - (&b_rand % &p).complete()).modulo(&p);
            let mul_check: Integer = ((a_rand % &p) * (b_rand % &p)) % &p;

            assert_eq!(added, added_check, "Failed for: {} + {}, and got {} instead of {}", a.get_rep(), b.get_rep(), added, added_check);
            assert_eq!(subtracted, sub_check, "Failed for: {} + {}, and got {} instead of {}", a.get_rep(), b.get_rep(), subtracted, sub_check);
            assert_eq!(multiplied, mul_check, "Failed for: {} * {}, and got {} instead of {}", a.get_rep(), b.get_rep(), multiplied, mul_check);

            let g: Arc<MultiplicativeGroup> = Arc::new(MultiplicativeGroup::from_finite_field(&f));

            // To not create zero
            let modification: Integer = Integer::from(rng.bits(32));
            let a_rand: Integer = (Integer::from(rng.bits(32)).modulo(&(&p-&Integer::ONE.clone()).complete())) + 1 + (&p * &modification);
            let b_rand: Integer = (Integer::from(rng.bits(32)).modulo(&(&p-&Integer::ONE.clone()).complete())) + 1 + (&p * modification);

            let a = Element::new(g.clone(), a_rand.clone());
            let b = Element::new(g.clone(), b_rand.clone());

            assert_eq!(a.get_rep(), &(&a_rand % &p).complete(), "Failed representation of {} in Z_{}*,  got {} instead of {}", &a_rand, &p, a.get_rep(), (&a_rand % &p).complete());
            assert_eq!(b.get_rep(), &(&b_rand % &p).complete(), "Failed representation of {} in Z_{}*,  got {} instead of {}", &b_rand, &p, b.get_rep(), (&b_rand % &p).complete());

            let multiplied = a.mul_ref(&b);
            let mul_check: Integer = ((&a_rand % &p).complete() * (&b_rand % &p).complete()) % &p;

            assert_eq!(multiplied.get_rep(), &mul_check);

            let divided = a.div_ref(&b);

            assert_eq!(divided.mul_ref(&b).get_rep(), &(&a_rand % &p).complete(), "Failed division {}/{} in Z_{}*", &a_rand, &b_rand, &p);


            prime.next_prime_mut();
        }
    }

    #[test]
    #[should_panic(expected = "Zero Division")]
    fn test_zero_division() {
        let mut rng = RandState::new();
        let p: Integer = Integer::from(rng.bits(32)).next_prime();
        let a_rand: Integer = Integer::from(rng.bits(32));
        let f: Arc<FiniteField> = Arc::new(FiniteField::new(p).unwrap());
        let a_elem: Element<FiniteField> = Element::new(f.clone(), a_rand);
        let b_elem: Element<FiniteField> = Element::new(f, Integer::ZERO);

        a_elem.div_ref(&b_elem);
    }

    #[test]
    fn test_extended_euclidean() {
        let n: u32 = 200;
        let mut prime: Integer = Integer::from(2);
        let mut rng = RandState::new();

        for _ in 2..n {
            let p = prime.clone();
            let f = Arc::new(FiniteField::new(p.clone()).unwrap());
            let a_rand = (Integer::from(rng.bits(32)).modulo(&(&p-&Integer::ONE.clone()).complete())) + 1 ;
            let a = Element::new(f.clone(), a_rand);

            let a_inv = a.mul_inv();
            assert_eq!(Integer::ONE, a_inv.mul_ref(&a).get_rep(), "We have a: {}, a_inv: {}, p: {}", a.get_rep(), a_inv.get_rep(), p);

            let g = Arc::new(MultiplicativeGroup::from_finite_field(&f));
            let a_rand = (Integer::from(rng.bits(32)).modulo(&(&p-&Integer::ONE.clone()).complete())) + 1 ;
            let a = Element::new(g.clone(), a_rand);
            let a_inv = a.mul_inv();
            assert_eq!(Integer::ONE, a_inv.mul_ref(&a).get_rep(), "We have a: {}, a_inv: {}, p: {}", a.get_rep(), a_inv.get_rep(), p);


            prime.next_prime_mut();

        }
    }


    #[test]
    fn test_exponentiation() {
        let n: u32 = 200;
        let mut rng = RandState::new();
        let mut prime: Integer = Integer::from(2);


        for _ in 2..n {
            let p = prime.clone();
            let f = Arc::new(FiniteField::new(p.clone()).unwrap());

            let a_rand = Integer::from(rng.bits(32));
            let a = Element::new(f.clone(), a_rand.clone());

            let x_rand = Integer::from(rng.bits(32));

            // let a_exp_check = a_rand.exp_residue(&x_rand, &p);
            let a_exp_check = a_rand.pow_mod_ref(&x_rand, &p).unwrap().complete();
            let a_exp = a.pow(&x_rand);

            assert_eq!(a_exp.get_rep(), &a_exp_check);

            prime.next_prime_mut();
        }
    }
}
