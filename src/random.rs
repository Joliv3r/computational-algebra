use rand::{thread_rng, Rng};
use rug::Integer;


pub fn randint_bits(bits: usize) -> Integer {
    let mut n = Integer::from(1);
    let mut rng = thread_rng();
    for _ in 1..bits {
        n <<= 1;
        n += rng.gen_range(0..=1);
    }
    n
}


pub fn randint_bits_odd(bits: usize) -> Integer {
    let mut n = randint_bits(bits-1);
    n <<= 1;
    n += 1;
    n
}


pub fn randint_digits(digits: usize) -> Integer {
    let mut rng = thread_rng();
    let mut n: Integer = Integer::from(rng.gen_range(1..=9));

    for _ in 0..digits {
        n *= 10;
        n += rng.gen_range(0..=9);
    }
    n
}
