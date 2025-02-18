#[cfg(test)]
mod test {
    use std::{io::{BufRead, BufReader}, str::FromStr, fs};
    use rug::{Integer, integer::IsPrime};
    use beralg::integers::prime::*;
    use rand::{thread_rng, Rng};

    #[test]
    fn test_primality_test_for_primes() {
        let file = fs::File::open("primes").unwrap();
        let reader = BufReader::new(file);
        let t = 30;
        let bound = 2000;

        for line in reader.lines() {
            let p: Integer = Integer::from_str(&line.unwrap()).unwrap();
            assert_eq!(true, is_likely_prime_with_trial_division(&p, t, bound), "Identified {} as non-prime", &p);
            assert_eq!(true, fermat_is_prime(&p, t), "Identified {} as non-prime using Fermat-test", &p);
            assert_eq!(true, rabin_miller_is_prime(&p, t), "Identified {} as non-prime using Rabin-Miller", &p);
        }
    }

    #[test]
    fn test_primality_test_for_non_primes() {
        let file = fs::File::open("non-primes").unwrap();
        let reader = BufReader::new(file);
        let t = 30;
        let bound = 2000;

        for line in reader.lines() {
            let n: Integer = Integer::from_str(&line.unwrap()).unwrap();
            assert_eq!(false, is_likely_prime_with_trial_division(&n, t, bound), "Identified {} as prime", &n);
            assert_eq!(false, fermat_is_prime(&n, t), "Identified {} as prime using Fermat-test", &n);
            assert_eq!(false, rabin_miller_is_prime(&n, t), "Identified {} as prime using Rabin-Miller", &n);
        }
    }

    #[test]
    fn test_find_prime_with_bit_length() {
        let reps = 25;
        let t = 30;
        let bound_td = 200;
        let bound_sieving = 30;
        let mut rng = thread_rng();
        for _ in 0..reps {
            let bits = rng.gen_range(5..200);
            let p = find_prime_with_bit_length(bits, t);
            assert!(p.is_probably_prime(t as u32) != IsPrime::No, "Found {} as prime", &p);
            assert_eq!(p.significant_bits(), bits as u32);
        }

        for _ in 0..reps {
            let bits = rng.gen_range(5..200);
            let trial_division = find_prime_with_bit_length_using_trial_division(bits, t, bound_td);
            assert!(trial_division.is_probably_prime(t as u32) != IsPrime::No, "Found {} as prime with trial_division", &trial_division);
            assert_eq!(trial_division.significant_bits(), bits as u32);
        }

        for _ in 0..reps {
            let bits = rng.gen_range(5..200);
            let sieving = find_prime_with_bit_length_using_sieving(bits, t, bound_sieving);
            assert!(sieving.is_probably_prime(t as u32) != IsPrime::No, "Found {} as prime with sieving", &sieving);
            assert_eq!(sieving.significant_bits(), bits as u32);
        }
    }
}
