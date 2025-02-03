// This is in addition to computational algebra a personal learning experience with rust.

use std::env;

use rand::{thread_rng, Rng};

pub mod algebraic_structure;
pub mod tests;
pub mod integer_computations;
pub mod prime;
pub mod random;


fn options_handler(args: Vec<String>) {
    for arg in &args[1..] {
        if arg == "generate-test-lists" {
            prime::generate_primes();
            prime::generate_non_primes();
        } else if arg == "generate-small-primes" {
            prime::generate_small_primes();
        } else {
            println!("Argument {} is not valid.", arg);
        }

    }
}

fn main() {
    options_handler(env::args().collect());
}
