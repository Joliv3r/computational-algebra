use beralg::factor::{file_handler::choose_primes, random_squares::{factorization_by_random_squares, find_factors_by_random_squares, find_multiple_relations, find_squares_by_relations, find_two_real_factors_by_random_squares, merge_tuples}};
use rug::Integer;

fn main() {
    let n = Integer::from(2003u64*1064u64*3539u64*8539u64);
    let number_of_relations = 200;

    let factors = factorization_by_random_squares(&n, number_of_relations, 0);
    println!("{} factors into: ", &n);
    for factor in factors {
        print!("{}, ", factor);
    }
    println!();
}
