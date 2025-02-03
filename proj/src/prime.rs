use rand::{thread_rng, Rng};
use rug::integer::IsPrime;
use rug::rand::RandState;
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
    let mut p: Integer = Integer::from(2);
    
    println!("Generating list of small primes...");
    for _ in 0..n {
        if let Err(e) = writeln!(file, "{}", &p) {
            eprintln!("Couldn't write to file: {}", e);
        }
        p.next_prime_mut();
    }
}


pub fn is_likely_prime_with_trial_division(candidate: &Integer, n: usize) -> bool {
    let small_primes = fs::File::open("small-primes").expect("small-primes file should have been generated");
    let reader = BufReader::new(small_primes);
    for prime in reader.lines() {
        let p: Integer = Integer::from_str(&prime.unwrap()).expect("All entries of small-primes should be integers.");
        if candidate > &p && candidate%p == 0 {
            return false
        }         
    }

    fermat_is_prime(candidate, n)
}


pub fn fermat_is_prime(n: &Integer, reps: usize) -> bool {
    let mut rng = RandState::new();
    for _ in 0..reps {
        let a = Integer::ONE + (n-Integer::ONE.clone()).random_below_ref(&mut rng).complete();
        if pow_rug(&a, &(n-Integer::ONE).complete(), &n) != 1 {
            return false
        }
    }
    true
}


// For n an odd prime with n-1 = 2^s * r with r odd and a in [1, n-1] such that n and a are relatively prime we have: 
//      a^r = 1 (mod n)    or    a^(2^j * r) = -1 (mod n), for j in [0, s-1]
pub fn rabin_miller_is_prime(n: &Integer, reps: usize) -> bool {
    if !n.get_bit(0) {
        return false;
    }

    let mut rng = RandState::new();
    let mut r: Integer = n.clone() - Integer::ONE;
    let mut s = 0;
    while !r.get_bit(0) {
        r = r >> 1;
        s += 1;
    }

    for _ in 0..reps {
        let a = Integer::from(2) + (n-Integer::from(4)).random_below(&mut rng);
        let mut y = pow_rug(&a, &r, n);

        if &y != Integer::ONE && y != (n-Integer::ONE).complete() {
            let mut j = 1;

            while j <= s-1 && y != (n-Integer::ONE).complete() {
                y = pow_rug(&y, &Integer::from(2), n);
                if y == 1 {
                    return false;
                }
                j += 1;
            }

            if y != (n-Integer::ONE).complete() {
                return false;
            }
        }
    }

    true
}


pub fn find_prime_with_bit_length(bits: usize, n: usize) -> Integer {
    let mut p: Integer = randint_bits(bits);
    while !is_likely_prime_with_trial_division(&p, n) {
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
        let t = 30;

        for line in reader.lines() {
            let p: Integer = Integer::from_str(&line.unwrap()).unwrap();
            assert_eq!(true, is_likely_prime_with_trial_division(&p, t), "Identified {} as non-prime", &p);
            assert_eq!(true, fermat_is_prime(&p, t), "Identified {} as non-prime using Fermat-test", &p);
            assert_eq!(true, rabin_miller_is_prime(&p, t), "Identified {} as non-prime using Rabin-Miller", &p);
        }
    }

    #[test]
    fn test_primality_test_for_non_primes() {
        let file = fs::File::open("non-primes").unwrap();
        let reader = BufReader::new(file);
        let t = 30;

        for line in reader.lines() {
            let n: Integer = Integer::from_str(&line.unwrap()).unwrap();
            assert_eq!(false, is_likely_prime_with_trial_division(&n, t), "Identified {} as prime", &n);
            assert_eq!(false, fermat_is_prime(&n, t), "Identified {} as prime using Fermat-test", &n);
            assert_eq!(false, rabin_miller_is_prime(&n, t), "Identified {} as prime using Rabin-Miller", &n);
        }
    }

    #[test]
    fn test_find_prime_with_bit_length() {
        let reps = 50;
        let t = 30;
        let mut rng = thread_rng();
        for _ in 0..reps {
            let bits = rng.gen_range(1..200);
            let p = find_prime_with_bit_length(bits, t);
            assert!(p.is_probably_prime(t as u32) != IsPrime::No, "Found {} as prime", &p);
            assert!(p.significant_bits() == bits as u32);
        }
    }
}
