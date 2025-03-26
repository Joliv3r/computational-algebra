use gauss_jordan_elimination::gauss_elimination_generic;
use itertools::Itertools;
use rug::{integer::IsPrime, rand::RandState, Complete, Integer};
use std::{collections::HashMap, fs, hash::Hash, io::{BufRead, BufReader}, ops::AddAssign, time::{SystemTime, UNIX_EPOCH}};
use crate::algebraic_structure::z2::Z2;
use num::traits::{Zero, One};


// Returns factors of the square.
fn find_one_relation(n: &Integer, rng: &mut RandState) -> (Integer, Vec<(u64, u64)>) {
    // println!("Trying to find a relation...");
    let mut t = n.random_below_ref(rng).complete();

    let t_factors = loop {
        if let Some(temp_factors) = trial_division(&t.clone().pow_mod(&Integer::from(2), n).expect("Square exists.")) {
            // break temp_factors.iter().filter(|(_, exponent)| exponent % 2 == 1).map(|(factor, _)| *factor).collect_vec()
            break temp_factors;
        }
        // println!("For some reason trial division does not give factors.");
        t = n.random_below_ref(rng).complete();
    };

    (t, t_factors)
}


pub fn find_multiple_relations(n: &Integer, m: usize) -> HashMap<Integer, Vec<(u64, u64)>> {
    let mut rng = RandState::new();
    rng.seed(&Integer::from(SystemTime::now().duration_since(UNIX_EPOCH).expect("We are not time travelling").as_secs()));
    let mut hashmap: HashMap<Integer, Vec<(u64, u64)>> = HashMap::with_capacity(m);
    // println!("Starting to find relations...");
    for _ in 0..m {
        // NOTE: might give the same integer and factors twice, however this is unlikely when n is
        // large.
        let (t, factors) = find_one_relation(n, &mut rng);
        hashmap.insert(t, factors);
    }
    hashmap
}


pub fn merge_tuples<S, T>(vector: &Vec<(S, T)>) -> Vec<(S, T)>
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


fn find_all_pivots(matrix: &Vec<Vec<Z2>>) -> Vec<bool> {
    let cols = matrix[0].len();
    let rows = matrix.len();
    let mut pivots: Vec<bool> = vec![false; cols];
    // let mut count = 0;
    for i in 0..rows {
        for j in 0..cols {
            if matrix[i][j].0 {
                pivots[j] = true;
                break
            }
        }
    }
    pivots
}


fn find_a_nonzero_solution(matrix: &Vec<Vec<Z2>>, size_of_solution: usize) -> Option<Vec<Z2>> {
    let mut solution = vec![Z2::zero(); size_of_solution];
    let pivots = find_all_pivots(matrix);

    let first_non_pivot = pivots.iter().position(|x| !x).unwrap_or(matrix.len());

    if first_non_pivot >= size_of_solution {
        return None;
    }

    solution[first_non_pivot] = Z2::one();
    let mut ones: Vec<usize> = vec![first_non_pivot];

    for i in (0..matrix.len()).rev() {
        let mut sum = false;
        for j in &ones {
            if matrix[i][*j].0 {
                sum ^= true;
            }
        }
        if sum {
            solution[i] = Z2::one();
            ones.push(i);
        }
    }

    Some(solution)

    // for i in 0..*pivots.last().expect("Pivot vector is non-empty.") {
    //     if i >= matrix.len() {
    //         break;
    //     }
    //     let mut count = false;
    //     for pivot in pivots {
    //         if matrix[i][*pivot].0 {
    //             count ^= true;
    //         }
    //     }
    //     if count {
    //         vector[i] = Z2::one();
    //     }
    // }
    // for pivot in pivots {
    //     vector[*pivot] = Z2::one();
    // }
    // Some(vector)
}

pub fn find_squares_by_relations(relations: &HashMap<Integer, Vec<(u64, u64)>>) -> Option<Vec<&Integer>> {
    let integers = relations.keys().collect_vec();
    let primes: Vec<u64> = relations.values().into_iter().flatten().filter(|(_, exp)| exp % 2 == 1 ).map(|(n, _)| *n as u64).unique().collect_vec();
    let mut cols = integers.len();
    let rows = primes.len();

    if primes.len() > integers.len() {
        cols = primes.len();
    }
    
    let mut choice_matrix: Vec<Vec<Z2>> = vec![vec![Z2::zero(); cols]; rows];
    for (i, integer) in integers.iter().enumerate() {
        for prime in relations.get(integer).expect("Key should exist.") {
            if prime.1 % 2 == 1 {
                if let Some(index) = primes.iter().position(|&p| p == prime.0) {
                    choice_matrix[index][i] = Z2::one();
                }
            }
        }
    }

    gauss_elimination_generic(&mut choice_matrix, gauss_jordan_elimination::GaussEliminationOption::PrepareReduce);

    if let Some(choice_vector) = find_a_nonzero_solution(&choice_matrix, integers.len()) {
        let mut relation_vector = Vec::new();
        for (i, choice) in choice_vector.iter().enumerate() {
            if i >= integers.len() {
                break;
            }
            if choice.0 {
                relation_vector.push(integers[i]);           
            }
        }
        Some(relation_vector)
    } else {
        None
    }
}


pub fn find_factors_by_random_squares(n: &Integer, number_of_relations: usize) -> Integer {
    let mut relations = find_multiple_relations(&n, number_of_relations);
    let square_vec = loop {
        // println!("Trying to find squares from relations.");
        if let Some(square) = find_squares_by_relations(&relations) {
            break square;
        }
        relations = find_multiple_relations(&n, number_of_relations);
    };

    for integer in relations.keys() {
        let mut product = Integer::ONE.clone();
        for factor in relations.get(integer).expect("Exists.") {
            product = product * factor.0.pow(factor.1 as u32);
        }
        assert_eq!(integer.clone().pow_mod(&Integer::from(2), &n).expect("Square exists."), product, "Failed checking {}, which then doesn't factor to {:#?}", integer, relations.get(integer).expect("Exists"))
    }

    let mut square1 = Integer::ONE.clone();
    let mut square2_vec: Vec<(u64, u64)> = Vec::new();
    for factor in square_vec {
        square1 = square1.clone() * factor;
        square1 = square1.clone() % n;

        for prime_exp in relations.get(factor).expect("Key exists.") {
            square2_vec.push(*prime_exp);
        }
    }

    square2_vec = merge_tuples(&square2_vec);

    let mut square2 = Integer::ONE.clone();
    for (prime, exp) in &square2_vec {
        square2 = square2 * Integer::from(*prime).pow_mod(&Integer::from(*exp), &n).expect("Square exists.");
        square2 = square2 % n;
    }

    square2_vec = square2_vec.iter().map(|(prime, exp)| (*prime, exp/2)).collect();

    let mut square2 = Integer::ONE.clone();

    for (prime, exp) in square2_vec {
        square2 = square2 * prime.pow(exp as u32);
        square2 = square2 % n;
    }

    // println!("Found: {}^2 = {}^2  (mod {})", &square1, &square2, &n);
    n.clone().gcd(&(&square2 - &square1).complete())
}


pub fn find_two_real_factors_by_random_squares(n: &Integer, number_of_relations: usize) ->  (Integer, Integer) {
    let mut factor = find_factors_by_random_squares(n, number_of_relations);
    while &factor == n || &factor == Integer::ONE {
        // println!("Unsuccessfully found a factor.");
        factor = find_factors_by_random_squares(n, number_of_relations);
    }
    
    let factor2 = (n / &factor).complete();
    (factor, factor2)
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


pub fn factorization_by_random_squares(n: &Integer, number_of_relations: usize, depth: usize) -> Vec<Integer> {
    println!("Entered depth {} and we are now factorizing {}", depth, n);
    if n.is_probably_prime(30) != IsPrime::No {
        return vec![n.clone()];
    }
    if let Some(trial_division_factors) = trial_division(n) {
        let mut factors = Vec::new();
        for (factor, exp) in trial_division_factors {
            factors.append(&mut vec![Integer::from(factor); exp as usize]);
        }
        return factors;
    }

    let (factor1, factor2) = find_two_real_factors_by_random_squares(n, number_of_relations);

    let mut factors1 = factorization_by_random_squares(&factor1, number_of_relations, depth+1);
    let mut factors2 = factorization_by_random_squares(&factor2, number_of_relations, depth+1);

    factors1.append(&mut factors2);
    factors1
}


#[cfg(test)]
mod tests {
    use std::{collections::HashSet, fs::File, hash::RandomState, os::unix::thread, process::{Command, Output}};
    use crate::factor::file_handler::largest_chosen_prime;

    use rand::{thread_rng, Rng};
    use rug::integer::IsPrime;

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


    #[test]
    fn test_find_square_from_relations() {
        let mut relations: HashMap<Integer, Vec<(u64, u64)>> = HashMap::new();
        relations.insert(Integer::from(34), vec![(2, 1), (17, 1)]);
        relations.insert(Integer::from(3247), vec![(17, 1), (191, 1)]);
        relations.insert(Integer::from(142), vec![(2, 1), (71, 1)]);
        relations.insert(Integer::from(1458), vec![(2, 1), (3, 6)]);
        relations.insert(Integer::from(1207), vec![(17, 1), (71, 1)]);

        let square = find_squares_by_relations(&relations).unwrap();

        let num1 = Integer::from(34);
        let num2 = Integer::from(142);
        let num3 = Integer::from(1207);
        let answer: Vec<&Integer> = vec![&num1, &num2, &num3];

        let answer_unordered: HashSet<_> = answer.iter().collect();
        let square_unordered: HashSet<_> = square.iter().collect();

        assert_eq!(square_unordered, answer_unordered);

        relations.remove_entry(&num1);

        let square = find_squares_by_relations(&relations);

        assert!(square.is_none());
    }

    
    #[test]
    fn test_find_factor_by_random_square() {
        let mut rng = thread_rng();
        let loops = 4;
        let max = 1e6 as u64;
        let number_of_relations = 40;

        for _ in 0..loops {
            let n = Integer::from(rng.gen_range(2..max));
            let factor = find_factors_by_random_squares(&n, number_of_relations);
            assert!((&n % &factor).complete() == 0, "Found a non-factor {} of {}", factor, n);
        }

    }


    #[test]
    fn test_factoring_into_two_by_random_squares() {
        let mut rng = RandState::new();
        let loops = 5;
        let max = Integer::from(1e9 as u64);
        let number_of_relations = 150;

        for _ in 0..loops {
            let n = loop {
                let m = max.clone().random_below(&mut rng);
                if m.is_probably_prime(30) == IsPrime::No {
                    break m;
                }
            };
            let (factor1, factor2) = find_two_real_factors_by_random_squares(&n, number_of_relations);
            assert!(&factor1 != Integer::ONE, "Found 1 as a factor");
            assert!(&factor1 != &n, "Found n as a factor");
            assert_eq!(&(&factor1*&factor2).complete(), &n, "Product of factors {}, {} was not {}", &factor1, &factor2, &n);
        }

    }


    fn get_random_prime(max: usize) -> u64 {
        let mut rng = thread_rng();
        let filepath = "small-primes";
        let file = File::open(filepath).expect("File should exist.");
        let reader = BufReader::new(file);
        let lines: Vec<u64> = reader.lines().collect_vec().iter().map(|x| x.as_ref().expect("Every line should exist.").parse::<u64>().expect("Every line should be an integer.")).collect_vec();
        // let number_of_lines = lines.len();
        let line_number = rng.gen_range(1..max);
        let prime: u64 = lines[line_number];
        prime
    }


    #[test]
    fn test_factoring_by_random_squares() {
        let mut rng = thread_rng();
        let loops = 3;
        let number_of_relations = 300;
        let number_of_primes = 4;
        let max = 400;

        for _ in 0..loops {
            let mut integer: u128 = 1;
            for _ in 0..number_of_primes {
                integer *= get_random_prime(max) as u128;
            }
            let factors = factorization_by_random_squares(&Integer::from(integer), number_of_relations, 0);

            let mut prod = Integer::ONE.clone();
            for factor in factors {
                prod = prod * factor;
            }

            assert_eq!(prod, Integer::from(integer));
        }
    }


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
