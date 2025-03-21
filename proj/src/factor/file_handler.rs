use std::{fs, io::{BufReader, BufRead, Write}};


pub fn largest_chosen_prime() -> u64 {
    let chosen_primes_path = "chosen-primes";
    if !fs::exists(chosen_primes_path).expect("Can't check existence of file.") {
        panic!("There must exist a generated set of small-primes");
    }

    let chosen_primes = fs::File::open(chosen_primes_path).expect("small-primes file should have been generated");
    let reader = BufReader::new(chosen_primes);
    let largest_prime = reader
        .lines().last()
        .expect("There should be lines in file")
        .expect("Line should be readable")
        .parse::<u64>()
        .expect("Last line should be an integer");

    largest_prime
}


pub fn choose_primes(n: usize) {
    let path = "chosen-primes";
    let small_primes_path = "small-primes";

    if !fs::exists(small_primes_path).expect("Can't check existence of file.") {
        panic!("There must exist a generated set of small-primes");
    }

    let mut file = fs::OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(true)
        .open(path)
        .unwrap();

    let small_primes = fs::File::open(small_primes_path).expect("small-primes file should have been generated");
    let reader = BufReader::new(small_primes);
    for prime in reader.lines() {
        let p = prime.expect("Line should exist.").parse::<u64>().expect("Should be an integer.");
        if p > n as u64 {
            break;
        }
        if let Err(e) = writeln!(file, "{}", p) {
            eprintln!("Error while reading file: {}", e);
        } 
    }
}
