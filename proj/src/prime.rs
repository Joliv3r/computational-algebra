use rand::{thread_rng, Rng};
use rug::integer::IsPrime;
use rug::{Complete, Integer};
use std::ops::Range;
use std::fs;
use std::io::{BufReader, BufRead, Write};
use std::str::FromStr;
use crate::random::{randint_bits, randint_digits};
use crate::integer_computations::pow_rug;


pub fn generate_primes() {
    println!("Opening file ./primes");
    // Rewrites file.
    let mut file = fs::OpenOptions::new()
        .write(true)
        .truncate(true)
        .open("primes")
        .unwrap();
    
    
    let mut p: Integer = Integer::ONE.clone();
    let number_of_primes = 1000;
    let space_of_primes = 100000;

    println!("Starting prime generating...");
    for _ in 0..number_of_primes {
        for _ in 0..space_of_primes {
            p.next_prime_mut();
        }
        if let Err(e) = writeln!(file, "{}", &p) {
            eprintln!("Couldn't write to file: {}", e)
        }
    }
}


pub fn generate_non_primes() {
    println!("Opening file ./non-primes");
    // Rewrites file.
    let mut file = fs::OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(true)
        .open("non-primes")
        .unwrap();

    let primes = fs::File::open("primes").expect("File should exist");
    let reader = BufReader::new(primes);

    println!("Generating list of non-primes...");
    for line in reader.lines() {
        let mut rng = thread_rng();
        let p: Integer = Integer::from_str(&line.unwrap()).expect("All entries of file should be numbers"); 
        let mut n: Integer = (&p + rng.gen_range::<usize, Range<usize>>(1..1000)).complete();
        while n.is_probably_prime(30) != IsPrime::No {
            n = (&p + rng.gen_range::<usize, Range<usize>>(1..1000)).complete();
        }

        if let Err(e) = writeln!(file, "{}", &n) {
            eprintln!("Couldn't write to file: {}", e);
        }
    }
}


pub fn generate_small_primes() {
    println!("Opening file ./small-primes");
    let mut file = fs::OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(true)
        .open("small-primes")
        .unwrap();

    let n = 200;
    let mut p: Integer = Integer::ONE.clone();
    
    println!("Generating list of small primes...");
    for _ in 0..n {
        if let Err(e) = writeln!(file, "{}", &p) {
            eprintln!("Couldn't write to file: {}", e);
        }
        p.next_prime_mut();
    }
}


pub fn is_likely_prime(candidate: &Integer, n: usize) -> bool {
    let small_primes = fs::File::open("small-primes").expect("small-primes file should have been generated");
    let reader = BufReader::new(small_primes);
    for prime in reader.lines() {
        let p: Integer = Integer::from_str(&prime.unwrap()).expect("All entries of small-primes should be integers.");
        if candidate > &p && p%candidate == 0 {
            return false
        }         
    }

    fermat_is_prime(candidate, n)
}


pub fn fermat_is_prime(p: &Integer, n: usize) -> bool {
    let mut rng = thread_rng();
    for _ in 0..n {
        let length = p.to_string().len();
        let digits = rng.gen_range(1..=length);
        let a: Integer = randint_digits(digits);
        if pow_rug(&a, &(p-Integer::ONE).complete(), &p) != 1 {
            return false
        }
    }
    true
}


pub fn find_prime_with_bit_length(bits: usize, n: usize) -> Integer {
    let mut p: Integer = randint_bits(bits);
    while !is_likely_prime(&p, n) {
        p = randint_bits(bits);
    }
    p
}


#[cfg(test)]
mod test {
    use std::{io::{BufRead, BufReader}, str::FromStr};

    use super::*;

    #[test]
    fn test_primality_test_for_primes() {
        let file = fs::File::open("primes").unwrap();
        let reader = BufReader::new(file);

        for line in reader.lines() {
            let p: Integer = Integer::from_str(&line.unwrap()).unwrap();
            assert_eq!(true, is_likely_prime(&p, 30), "Identified {} as non-prime", &p);
        }
    }

    #[test]
    fn test_primality_test_for_non_primes() {
        let file = fs::File::open("non-primes").unwrap();
        let reader = BufReader::new(file);

        for line in reader.lines() {
            let n: Integer = Integer::from_str(&line.unwrap()).unwrap();
            assert_eq!(false, is_likely_prime(&n, 30), "Identified {} as prime", &n);
        }
    }

    #[test]
    fn test_find_prime_with_bit_length() {
        let reps = 50;
        let mut rng = thread_rng();
        for _ in 0..reps {
            let bits = rng.gen_range(1..200);
            let p = find_prime_with_bit_length(bits, 10);
            assert!(p.is_probably_prime(50) != IsPrime::No);
            assert!(p.significant_bits() == bits as u32);
        }
    }
}
