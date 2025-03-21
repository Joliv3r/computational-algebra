use rug::{rand::RandState, Complete, Integer};
use std::{collections::HashMap, fs, io::{BufRead, BufReader}};


fn find_relation(n: &Integer) -> (Integer, Vec<(u64, u64)>) {
    let mut rng = RandState::new();
    let mut t = n.random_below_ref(&mut rng).complete();

    let t_factors = loop {
        if let Some(temp_factors) = trial_division(&t) {
            break temp_factors
        }
        t = n.random_below_ref(&mut rng).complete();
    };

    (t, t_factors)
}


fn find_relations(n: &Integer, m: u64) -> HashMap<Integer, Vec<(u64, u64)>> {
    let mut hashmap: HashMap<Integer, Vec<(u64, u64)>> = HashMap::with_capacity(m as usize);
    for _ in 0..m {
        let (t, factors) = find_relation(n);
        hashmap.insert(t, factors);
    }
    hashmap
}


fn check_if_relations_gives_squares(relations: &HashMap<Integer, Vec<(u64, u64)>>) -> bool {
    for 
    todo!()
}


fn trial_division(t: &Integer) -> Option<Vec<(u64, u64)>> {
    let chosen_primes_path = "chosen-primes";
    if !fs::exists(chosen_primes_path).expect("Can't check existence of file.") {
        panic!("There must exist a generated set of small-primes");
    }

    let mut t_clone = t.clone();
    let mut factors: Vec<(u64, u64)> = Vec::new();

    let chosen_primes = fs::File::open(chosen_primes_path).expect("small-primes file should have been generated");
    let reader = BufReader::new(chosen_primes);
    for prime in reader.lines() {
        let p = prime.expect("Line should exist").parse::<u64>().expect("Should be a positive integer");
        while (&t_clone%p).complete() == 0 {
            if factors.len() == 0 {
                factors.push((p, 1));
            } else if factors.last().expect("Vector is not empty").0 == p {
                factors.last_mut().expect("Vector is not empty").1 += 1;
            } else {
                factors.push((p, 1));
            }
            t_clone = t_clone/p;
        }
    }
    
    if t_clone == 1 {
        Some(factors)
    } else {
        None
    }
}


#[cfg(test)]
mod tests {
    use std::process::{Command, Output};
    use crate::factor::file_handler::largest_chosen_prime;

    use rand::{thread_rng, Rng};

    use super::*;

    fn largest_factor(t: &Integer) -> u64 {
        let factors: Output = Command::new("factor")
            .arg(t.to_string())
            .output()
            .expect("should always work");
        let binding = String::from_utf8(factors.stdout).expect("Should exist");
        let binding = binding.strip_suffix("\n").expect("The string always ends with endline.");
        let binding: Vec<&str> = binding.split(" ").collect();
        let largest_factor: u64 = binding
            .into_iter()
            .filter_map(|s| s.parse::<u64>().ok())
            .into_iter()
            .max()
            .expect("All entries are comparable.");

        largest_factor
    }

    fn check_factors(t: &Integer, factors: &Vec<(u64, u64)>) {
        let mut product = Integer::ONE.clone();

        for (p, i) in factors {
            product = product * p.pow(*i as u32);
        }

        assert_eq!(&product, t);
    }

    #[test]
    fn test_trial_division() {
        let loops = 200;
        let mut rng = thread_rng();

        let largest_prime = largest_chosen_prime();

        for _ in 0..loops {
            let t = Integer::from(rng.gen_range(1..u64::MAX));

            if let Some(factors) = trial_division(&t) {
                check_factors(&t, &factors);
            } else {
                let largest_factor = largest_factor(&t);
                if largest_factor <= largest_prime {
                    panic!("Should have been able to factor {}", t);
                }
            }
        }
    }
}
