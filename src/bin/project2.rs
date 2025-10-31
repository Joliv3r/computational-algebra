#![allow(dead_code)]
use beralg::integers::prime::{find_prime_in_interval_with_sieving, find_prime_with_bit_length, find_prime_with_bit_length_using_sieving, find_prime_with_bit_length_using_trial_division, is_likely_prime_with_trial_division, rabin_miller_is_prime};
use beralg::random::{randint_bits, randint_bits_odd};
use rug::Integer;
use std::time::Instant;
use std::io::{BufRead, BufReader, Write};
use std::fs;
use std::str::FromStr;


fn count_candidates_in_find_prime_with_bit_length(bits: usize, t: usize, bound: usize) -> usize {
    let mut p: Integer = randint_bits_odd(bits);
    let mut count = 1;
    while !is_likely_prime_with_trial_division(&p, t, bound) {
        count += 1;
        p = randint_bits_odd(bits);
    }
    count
}


fn count_candidates_in_find_prime_with_bit_length_avg(bits: usize, n: usize, t: usize, bound: usize) -> f64 {
    let mut count = 0;
    for _ in 0..n {
        count += count_candidates_in_find_prime_with_bit_length(bits, t, bound);
    }
    (count as f64)/(n as f64)
}


fn width_in_random_interval_search(bits: usize, probability: f64) -> usize {
    let mut d = 1;
    while (1f64 - 1f64/((bits as f64)*2f64.ln())).powf(d as f64) > 1f64 - probability {
        d += 1;
    }
    d
}


fn expected_candidates(bits: usize) -> f64 {
    (bits as f64)*2f64.ln()
}


fn expected_candidates_with_filtration(bits: usize, bound: usize) -> f64 {
    let small_primes = fs::File::open("small-primes").expect("small-primes file should have been generated");
    let reader = BufReader::new(small_primes);
    let mut prod: f64 = 1.;
    
    for prime in reader.lines() {
        let p: Integer = Integer::from_str(&prime.unwrap()).expect("All entries of small-primes should be integers.");
        if p < bound {
            prod *= (p.to_f64()-1.)/p.to_f64();
        }
    }

    (bits as f64)*2f64.ln()*prod
}


fn time_finding_primes() {
    let bits = 300;
    let loops = 80;
    let t = 20;

    let mut number_of_primes_for_trial_division = 0;
    let increase = 50;
    let number_of_increases = 40;
    let mut timing_vector = Vec::with_capacity(number_of_increases);


    for _ in 0..number_of_increases {
        number_of_primes_for_trial_division += increase;
        // generate_small_primes(number_of_primes_for_trial_division);
        let now = Instant::now();
        for _ in 0..loops {
            find_prime_with_bit_length(bits, t);
        }
        let elapsed = now.elapsed()/loops;
        timing_vector.push((number_of_primes_for_trial_division, elapsed.as_micros()))
    }

    for (num_of_precomputations, elapsed_as_millis) in timing_vector.iter() {
        println!("{}, {}", num_of_precomputations, elapsed_as_millis);
    }
}


fn generate_list_of_composite_until_prime(bits: usize, append: bool) {
    println!("Opening file ./trial-division");
    let mut file = fs::OpenOptions::new()
        .write(true)
        .truncate(!append)
        .append(append)
        .open("trial-division")
        .unwrap();

    let t = 30;

    loop {
        let p = randint_bits(bits);

        if let Err(e) = writeln!(file, "{}", &p) {
            eprintln!("Couldn't write to file: {}", e)
        }
        
        if rabin_miller_is_prime(&p, t) {
            break;
        }
    }
}


fn generate_random_numbers(bits: usize, number: usize) {
    let mut file = fs::OpenOptions::new()
        .write(true)
        .truncate(true)
        .open("random-numbers")
        .unwrap();

    for _ in 0..number {
        let a = randint_bits(bits);
        if let Err(e) = writeln!(file, "{}", &a) {
            eprintln!("Couldn't write to file: {}", e)
        }
    }
}


fn find_good_number_of_precomputed_primes(bits: usize) {
    let increment = 50;
    let stop = 1000;
    let loops = 40;
    let t = 30;
    let security = 5;

    generate_list_of_composite_until_prime(bits, false);
    for _ in 0..security {
        generate_list_of_composite_until_prime(bits, true);
    }

    let mut bound = 0;

    while bound <= stop {
        let now = Instant::now();
        for _ in 0..loops {
            let list_until_prime = fs::File::open("trial-division").expect("trial-division file should have been generated");
            let reader = BufReader::new(list_until_prime);
            for prime in reader.lines() {
                let p = Integer::from_str(&prime.expect("Line exists")).expect("Every line should be an integer");
                is_likely_prime_with_trial_division(&p, t, bound);
            }
        }
        let elapsed = now.elapsed()/(loops*security) as u32;
        println!("Bound: {}, Elapsed: {}", bound, elapsed.as_micros());
        bound += increment;
    }
}


fn find_good_number_of_primes_for_sieving(bits: usize) {
    let increment = 1000;
    let start = 2000;
    let stop = 17000;
    let loops = 40;
    let t = 30;
    let number_of_generated = 40;

    let mut bound = start;
    let d = width_in_random_interval_search(bits, 0.95);

    println!("Found width to be {}", d);

    generate_random_numbers(bits, number_of_generated);

    println!("Done generating lists");

    while bound <= stop {
        // println!("Checking loops now");
        let now = Instant::now();
        for _ in 0..loops {
            // find_prime_with_bit_length_using_sieving(bits, t, bound);
            let list_until_prime = fs::File::open("random-numbers").expect("trial-division file should have been generated");
            let reader = BufReader::new(list_until_prime);
            for number in reader.lines() {
                let a = Integer::from_str(&number.expect("Line should exist")).expect("Lines should be integers");
                // println!("Checking for a={}", a);
                find_prime_in_interval_with_sieving(&a, d, t, bound);
            }
        }
        let elapsed = now.elapsed()/(loops*number_of_generated) as u32;
        println!("Bound: {}, Elapsed: {}", bound, elapsed.as_micros());
        bound += increment;
    }
}


fn time_finding_primes_trial_division_vs_sieving() {
    let bits = 500;
    let loops = 1000;
    let bound_td = 150;
    let bound_sieving = 12000;
    let t = 30;

    let now = Instant::now();
    for _ in 0..loops {
        find_prime_with_bit_length_using_trial_division(bits, t, bound_td);
    }
    let elapsed_td = now.elapsed().as_micros()/loops;
    println!("trial-division: {}", elapsed_td);

    let now = Instant::now();
    for _ in 0..loops {
        find_prime_with_bit_length_using_sieving(bits, t, bound_sieving);
    }
    let elapsed_sieving = now.elapsed().as_micros()/loops;

    println!("sieving: {}", elapsed_sieving);
}


fn main() {
    time_finding_primes_trial_division_vs_sieving();
}
