use rug::{Complete, Integer};
use crate::algebraic_structure::{Element, HasRepresentation};


pub fn extended_euclidean_ordered(a: &Integer, b: &Integer) -> (Integer, Integer, Integer) {
    // There have been no attempt to optimize this function using cloning and references
    // efficiently.
    let mut a1: Integer = a.clone();
    let mut b1: Integer = b.clone();

    if b.is_zero() {
        return (a1, Integer::ONE.clone(), Integer::ZERO.clone())
    };

    let mut x2: Integer = Integer::ONE.clone();
    let mut x1: Integer = Integer::ZERO.clone();
    let mut y2: Integer = Integer::ZERO.clone();
    let mut y1: Integer = Integer::ONE.clone();


    while b1 > Integer::ZERO {
        // q = a1/b1;
        let (q, _) = (a1.div_rem_floor_ref(&b1)).complete();
        let r = a1 - &q*&b1;
        let x = x2 - &q*&x1;
        let y = y2 - &q*&y1;
        a1 = b1;
        b1 = r;
        x2 = x1;
        x1 = x;
        y2 = y1;
        y1 = y;
    };

    (a1, x2, y2)
}


pub fn extended_euclidean_to_integers<T: HasRepresentation + Clone>(a: &Element<T>, b: &Element<T>) -> (Integer, Integer, Integer) {
    let a_rep: &Integer = a.get_rep();
    let b_rep: &Integer = b.get_rep();
    if a_rep > b_rep {
        extended_euclidean_ordered(a_rep, b_rep)
    } else {
        extended_euclidean_ordered(b_rep, a_rep)
    }
}


pub fn pow_rug(a: &Integer, b: &Integer, n: &Integer) -> Integer {
    let mut product: Integer = Integer::ONE.clone();
    let mut base = a.clone() % n;
    let mut exponent = b.clone();

    while exponent != 0 {
        if exponent.get_bit(0) {
            product = (product * &base).modulo(n);
        }
        base = base.square().modulo(n);
        exponent = exponent >> 1;
    }
    product
}


pub fn naive_pow(a: &Integer, b: &Integer, n: &Integer) -> Integer {
    let mut product: Integer = Integer::ONE.clone();
    let s = b.to_u64().expect("The number is WAY too high to naively calculate.");
    for _ in 0..s {
        product *= a;
        product %= n;
    }
    product
}
