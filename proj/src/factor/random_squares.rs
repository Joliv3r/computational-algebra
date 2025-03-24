use gauss_jordan_elimination::gauss_jordan_elimination_generic;
use itertools::Itertools;
use ndarray::{array, Array2};
use rug::{rand::RandState, Complete, Integer};
use std::{collections::HashMap, fs::{self, File}, hash::Hash, io::{BufRead, BufReader}, ops::{AddAssign, Index}};
use crate::algebraic_structure::z2::Z2;


// We do not care for the factorization, and therefore only return the factors that are repeated an
// odd number of times.
fn find_one_relation(n: &Integer) -> (Integer, Vec<u64>) {
    let mut rng = RandState::new();
    let mut t = n.random_below_ref(&mut rng).complete().pow_mod(&Integer::from(2), n).expect("The square should exist.");

    let t_factors = loop {
        if let Some(temp_factors) = trial_division(&t) {
            break temp_factors.iter().filter(|(_, exponent)| exponent % 2 == 1).map(|(factor, _)| *factor).collect_vec()
        }
        t = n.random_below_ref(&mut rng).complete();
    };

    (t, t_factors)
}


fn find_multiple_relations(n: &Integer, m: u64) -> HashMap<Integer, Vec<u64>> {
    let mut hashmap: HashMap<Integer, Vec<u64>> = HashMap::with_capacity(m as usize);
    for _ in 0..m {
        // NOTE: might give the same integer and factors twice, however this is unlikely when n is
        // large.
        let (t, factors) = find_one_relation(n);
        hashmap.insert(t, factors);
    }
    hashmap
}


fn merge_tuples<S, T>(vector: &Vec<(S, T)>) -> Vec<(S, T)>
where 
    S: Eq + Hash + Copy,
    T: AddAssign + Default + Copy
{
    let mut map: HashMap<S, T> = HashMap::with_capacity(vector.len());

    for (u, v) in vector.iter() {
        *map.entry(*u).or_insert_with(Default::default) += *v;
    }

    map.into_iter().collect()
}


fn find_first_non_pivot(matrix: &Vec<Vec<Z2>>) -> Option<usize> {
    for (i, row) in matrix.iter().enumerate() {
        if !row[i].0 {
            return Some(i);
        }
    }
    if matrix.len() < matrix[0].len() {
        Some(matrix.len())
    } else {
        None
    }
}


fn find_a_nonzero_solution(matrix: &Vec<Vec<Z2>>) -> Option<Vec<Z2>> {
    if let Some(non_pivot) = find_first_non_pivot(matrix) {
        let mut vector = vec![Z2::zero(); matrix[0].len()];
        for i in 0..non_pivot {
            if matrix[i][non_pivot].0 {
                vector[i] = Z2::one();
            }
        }
        vector[non_pivot] = Z2::one();
        Some(vector)
    } else {
        None
    }
}

fn find_squares_by_relations(relations: &HashMap<Integer, Vec<u64>>) -> Option<Vec<Integer>> {
    let integers = relations.keys().collect_vec();
    let primes: Vec<u64> = relations.values().into_iter().flatten().unique().map(|n| *n as u64).collect_vec();
    
    let mut choice_matrix: Vec<Vec<Z2>> = vec![vec![Z2::zero(); integers.len()]; primes.len()];

    for (i, integer) in integers.iter().enumerate() {
        for prime in relations.get(integer).expect("Key should exist.") {
            let index = primes.iter().position(|p| p == prime).expect("Prime should be in the vector.");
            choice_matrix[index][i] = Z2::one();
        }
    }

    gauss_jordan_elimination_generic(&mut choice_matrix);
    
    let relation_vector = Vec::new();

    if let Some(choice_vector) = find_a_nonzero_solution(&choice_matrix) {
        for (i, choice) in choice_vector.iter().enumerate() {
            if choice.0 {
                relation_vector.push(integers[i]);           
            }
        }
        Some(relation_vector)
    } else {
        None
    }
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
    use std::{collections::HashSet, process::{Command, Output}};
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

    // #[test]
    // fn test_find_relation() {
    //     let loops = 1;
    //     let n = 5;
    //     let mut rng = thread_rng();
    //
    //     for _ in 0..loops {
    //         let t = Integer::from(rng.gen_range(1..u64::MAX));
    //         let relations = find_multiple_relations(&t, n);
    //
    //         for (number, factors) in relations.into_iter() {
    //             check_factors(&number, &factors);
    //             if number > t {
    //                 panic!("{} should be smaller than {}", number, t)
    //             }
    //         }
    //     }
    // }

    // #[test]
    // fn test_check_if_relation_is_square() {
    //     let mut rng = thread_rng();
    //     let n = 5;
    //     let mut relations: HashMap<Integer, Vec<(u64, u64)>> = HashMap::with_capacity(n);
    //     let mut integers: Vec<Integer> = Vec::with_capacity(n);
    //
    //     for _ in 0..n {
    //         let mut s = Integer::from(rng.gen_range(1..u32::MAX)).square();
    //         let factors = loop {
    //             if let Some(factors) = trial_division(&s) {
    //                 break factors
    //             }
    //             s = Integer::from(rng.gen_range(1..u32::MAX)).square();
    //         };
    //         integers.push(s.clone());
    //         relations.insert(s.clone(), factors);
    //     }
    //
    //     assert!(check_if_relations_gives_squares(&integers, &relations))
    // }

    #[test]
    fn test_merge_tuples() {
        let vec = vec![(1, 2), (1, 5), (5, 2), (7, 2), (5, 8), (1, 3), (7, 7), (5, 0)];
        let answer = vec![(1, 10), (5, 10), (7, 9)];
        let merged_vec = merge_tuples(&vec);

        let answer_unordered: HashSet<_> = answer.iter().collect();
        let merged_unordered: HashSet<_> = merged_vec.iter().collect();

        assert_eq!(merged_unordered, answer_unordered);
    }
}
