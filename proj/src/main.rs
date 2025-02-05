// This is in addition to computational algebra a personal learning experience with rust.

use std::{env, time::Instant};
use std::io::{BufReader, BufRead, Write};
use std::fs;

use prime::{count_candidates_in_find_prime_with_bit_length_avg, find_prime_with_bit_length, find_prime_with_bit_length_using_interval, generate_primes, rabin_miller_is_prime};
use rug::Float;
use rug::{ops::{AddFrom, PowAssign}, rand::RandState, Complete, Integer};
use std::str::FromStr;

pub mod algebraic_structure;
pub mod tests;
pub mod integer_computations;
pub mod prime;
pub mod random;


fn options_handler(args: Vec<String>) -> bool {
    let mut used = false;
    for arg in &args[1..] {
        if arg == "generate-test-lists" {
            prime::generate_primes();
            prime::generate_non_primes();
            used = true;
        } else if arg == "generate-small-primes" {
            let n = 500;
            prime::generate_small_primes(n);
            used = true;
        } else if arg == "proj1" {
            let loops = 10;
            let naive_square_points: usize = 20;
            let naive_points: usize = 50;
            let square_points: usize = 50;
            tests::plot_timing_naive_square(naive_square_points, loops).expect("Should not fail");
            tests::plot_timing_naive(naive_points, loops).expect("Should not fail");
            tests::plot_timing_square(square_points, loops).expect("Should not fail");
        } else if arg == "proj2" {
            continue;
        } else {
            println!("Argument {} is not valid.", arg);
        }
    }
    used
}


fn main() {
    if !options_handler(env::args().collect()) {
        // let loops = 1000;
        // let mut count = 0;
        //
        // let bits = 500;
        // let d = 240;
        // let t = 50;
        // 
        // for _ in 0..loops {
        //     if let Some(_) = find_prime_with_bit_length_using_interval(bits, d, t) {
        //         count += 1;
        //     }
        // }
        // println!("We found {} primes out of {} tries.", count, loops);

        // let bits = 500;
        // let loops = 1000;
        // let t = 50;
        // let avg = count_candidates_in_find_prime_with_bit_length_avg(bits, loops, t);
        // println!("The average amount of candidates tried for n={}, was {}", bits, avg);


        let small_primes = fs::File::open("small-primes").expect("small-primes file should have been generated");
        let reader = BufReader::new(small_primes);
        let mut prod: f64 = 1.;
        
        for prime in reader.lines() {
            let p: Integer = Integer::from_str(&prime.unwrap()).expect("All entries of small-primes should be integers.");
            if p < 2000 {
                prod *= (p.to_f64()-1.)/p.to_f64();
            }
        }
        println!("{}", prod);
    }
}
