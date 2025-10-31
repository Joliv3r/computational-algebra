use rand::{thread_rng, Rng};
use rug::rand::RandState;
use rug::{Complete, Integer};
use std::fs;
use std::io::{BufReader, BufRead};
use std::str::FromStr;
use crate::random::{randint_bits_odd, randint_bits};
use crate::integers::integer_computations::pow_rug;


pub fn is_likely_prime_with_trial_division(candidate: &Integer, n: usize, bound: usize) -> bool {
    if bound == 0 {
        return rabin_miller_is_prime(candidate, n);
    }
    let small_primes = fs::File::open("small-primes").expect("small-primes file should have been generated");
    let reader = BufReader::new(small_primes);
    for prime in reader.lines() {
        let p: Integer = Integer::from_str(&prime.unwrap()).expect("All entries of small-primes should be integers.");
        if p > bound {
            break;
        }
        if candidate > &p && candidate%p == 0 {
            return false
        }         
    }

    rabin_miller_is_prime(candidate, n)
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


// For n an odd prime with n-1 = 2^s * r with r odd and a in [1, n-1] we have: 
//      a^r = 1 (mod n)    or    a^(2^j * r) = -1 (mod n), for j in [0, s-1]
pub fn rabin_miller_is_prime(n: &Integer, reps: usize) -> bool {
    if *n == 2 {
        return true;
    } else if *n == 3 {
        return true;
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


pub fn find_prime_with_bit_length(bits: usize, t: usize) -> Integer {
    let mut p: Integer = randint_bits_odd(bits);
    while !rabin_miller_is_prime(&p, t) {
        p = randint_bits_odd(bits);
    }
    p
}


pub fn find_prime_with_bit_length_using_trial_division(bits: usize, t: usize, bound: usize) -> Integer {
    let mut p: Integer = randint_bits_odd(bits);
    while !is_likely_prime_with_trial_division(&p, t, bound) {
        p = randint_bits_odd(bits);
    }
    p
}


pub fn find_prime_with_bit_length_using_interval(bits: usize, d: usize, t: usize, bound: usize) -> Option<Integer> {
    let mut n = randint_bits(bits);
    if is_likely_prime_with_trial_division(&n, t, bound) {
        return Some(n)
    }
    for _ in 0..d {
        n += 1;
        if is_likely_prime_with_trial_division(&n, t, bound) {
            return Some(n)
        }
    }

    None
}


pub fn find_prime_in_interval_with_sieving(a: &Integer, d: usize, t: usize, bound: usize) -> Option<Integer> {
    let small_primes = fs::File::open("small-primes").expect("small-primes file should have been generated");
    let reader = BufReader::new(small_primes);

    let mut vec: Vec<bool> = vec![true; d];
    let mut capacity = d;

    for prime in reader.lines() {
        let p = prime.expect("Line should exist").parse::<usize>().expect("Lines should be primes");
        if p > bound {
            break;
        }

        let off_set: usize = (p - (a%p).complete()).to_usize().expect("Should be small enough");
        let mut count = 0;
        loop {
            let index: usize = count*p + &off_set;
            if index >= d {
                break;
            }

            if vec[index] {
                vec[index] = false;
                capacity -= 1;
            }
            count += 1;
        }
    }

    let mut sieving_vec: Vec<usize> = Vec::with_capacity(capacity);   

    for (i, j) in vec.iter().enumerate() {
        if *j {
            sieving_vec.push(i);
        }
    }

    if sieving_vec.len() == 0 {
        return None
    }

    let mut rng = thread_rng();
    let mut index = rng.gen_range(0..capacity);
    let mut p: Integer = (a + sieving_vec[index]).into();

    while !rabin_miller_is_prime(&p, t) {
        sieving_vec.remove(index);
        if sieving_vec.len() == 0 {
            return None
        }
        capacity -= 1;
        p = (a + sieving_vec[rng.gen_range(0..capacity)]).into();
        index = rng.gen_range(0..capacity);
    }
    Some(p)
}


fn approx_width_in_random_interval_search(bits: usize, probability: f64) -> usize {
    ( - probability/( (1f64 - 1f64/((bits as f64)*2f64.ln())) ).ln() ) as usize
}


pub fn find_prime_with_bit_length_using_sieving(bits: usize, t: usize, bound: usize) -> Integer {
    if bound == 0 {
        find_prime_with_bit_length(bits, t);
    }
    let probability = 0.95;
    let d = approx_width_in_random_interval_search(bits, probability);
    
    loop {
        let a = randint_bits(bits);
        if let Some(p) = find_prime_in_interval_with_sieving(&a, d, t, bound) {
            return p
        }
    }
}



