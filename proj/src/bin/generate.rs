use rand::{thread_rng, Rng};
use rug::integer::IsPrime;
use rug::{Complete, Integer};
use std::ops::Range;
use std::fs;
use std::io::{BufReader, BufRead, Write};
use std::str::FromStr;


fn generate_primes() {
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


fn generate_non_primes() {
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


fn generate_small_primes(n: usize) {
    println!("Opening file ./small-primes");
    let mut file = fs::OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(true)
        .open("small-primes")
        .unwrap();

    let mut p: Integer = Integer::from(2);
    
    println!("Generating list of small primes...");
    while p < n {
        if let Err(e) = writeln!(file, "{}", &p) {
            eprintln!("Couldn't write to file: {}", e);
        }
        p.next_prime_mut();
    }
}


fn main() {
    generate_primes();
    generate_non_primes();
    generate_small_primes(100000);
}
