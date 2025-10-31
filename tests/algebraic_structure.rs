#[cfg(test)]
mod algebraic_structure {
    use beralg::algebraic_structure::Element;
    use beralg::algebraic_structure::finite_field::{MultiplicativeGroup, FiniteField};
    use rug::{Integer, Complete, rand::RandState};
    use std::sync::Arc;

    #[test]
    fn test_algebraic_structure_arithmetic() {
        let n: u32 = 200;
        let mut prime: Integer = Integer::from(2);
        let mut rng = RandState::new();

        for _ in 2..n {
            let p = prime.clone();
            let f: Arc<FiniteField> = Arc::new(FiniteField::new(p.clone()).unwrap());
            let a_rand: Integer = Integer::from(rng.bits(32));
            let b_rand: Integer = Integer::from(rng.bits(32));
            let a = Element::new(f.clone(), a_rand.clone());
            let b = Element::new(f.clone(), b_rand.clone());

            assert_eq!(a.get_rep(), &(&a_rand % &p).complete());
            assert_eq!(b.get_rep(), &(&b_rand % &p).complete());

            let added = (a.add_ref(&b)).get_rep().clone();
            let subtracted = (a.sub_ref(&b)).get_rep().clone();
            let multiplied = (a.mul_ref(&b)).get_rep().clone();

            let added_check: Integer = ((&a_rand % &p).complete() + (&b_rand % &p).complete()) % &p;
            let sub_check: Integer = ((&a_rand % &p).complete() - (&b_rand % &p).complete()).modulo(&p);
            let mul_check: Integer = ((a_rand % &p) * (b_rand % &p)) % &p;

            assert_eq!(added, added_check, "Failed for: {} + {}, and got {} instead of {}", a.get_rep(), b.get_rep(), added, added_check);
            assert_eq!(subtracted, sub_check, "Failed for: {} + {}, and got {} instead of {}", a.get_rep(), b.get_rep(), subtracted, sub_check);
            assert_eq!(multiplied, mul_check, "Failed for: {} * {}, and got {} instead of {}", a.get_rep(), b.get_rep(), multiplied, mul_check);

            let g: Arc<MultiplicativeGroup> = Arc::new(MultiplicativeGroup::from_finite_field(&f));

            // To not create zero
            let modification: Integer = Integer::from(rng.bits(32));
            let a_rand: Integer = (Integer::from(rng.bits(32)).modulo(&(&p-&Integer::ONE.clone()).complete())) + 1 + (&p * &modification);
            let b_rand: Integer = (Integer::from(rng.bits(32)).modulo(&(&p-&Integer::ONE.clone()).complete())) + 1 + (&p * modification);

            let a = Element::new(g.clone(), a_rand.clone());
            let b = Element::new(g.clone(), b_rand.clone());

            assert_eq!(a.get_rep(), &(&a_rand % &p).complete(), "Failed representation of {} in Z_{}*,  got {} instead of {}", &a_rand, &p, a.get_rep(), (&a_rand % &p).complete());
            assert_eq!(b.get_rep(), &(&b_rand % &p).complete(), "Failed representation of {} in Z_{}*,  got {} instead of {}", &b_rand, &p, b.get_rep(), (&b_rand % &p).complete());

            let multiplied = a.mul_ref(&b);
            let mul_check: Integer = ((&a_rand % &p).complete() * (&b_rand % &p).complete()) % &p;

            assert_eq!(multiplied.get_rep(), &mul_check);

            let divided = a.div_ref(&b);

            assert_eq!(divided.mul_ref(&b).get_rep(), &(&a_rand % &p).complete(), "Failed division {}/{} in Z_{}*", &a_rand, &b_rand, &p);


            prime.next_prime_mut();
        }
    }

    #[test]
    #[should_panic(expected = "Zero Division")]
    fn test_zero_division() {
        let mut rng = RandState::new();
        let p: Integer = Integer::from(rng.bits(32)).next_prime();
        let a_rand: Integer = Integer::from(rng.bits(32));
        let f: Arc<FiniteField> = Arc::new(FiniteField::new(p).unwrap());
        let a_elem: Element<FiniteField> = Element::new(f.clone(), a_rand);
        let b_elem: Element<FiniteField> = Element::new(f, Integer::ZERO);

        a_elem.div_ref(&b_elem);
    }

    #[test]
    fn test_extended_euclidean() {
        let n: u32 = 200;
        let mut prime: Integer = Integer::from(2);
        let mut rng = RandState::new();

        for _ in 2..n {
            let p = prime.clone();
            let f = Arc::new(FiniteField::new(p.clone()).unwrap());
            let a_rand = (Integer::from(rng.bits(32)).modulo(&(&p-&Integer::ONE.clone()).complete())) + 1 ;
            let a = Element::new(f.clone(), a_rand);

            let a_inv = a.mul_inv();
            assert_eq!(Integer::ONE, a_inv.mul_ref(&a).get_rep(), "We have a: {}, a_inv: {}, p: {}", a.get_rep(), a_inv.get_rep(), p);

            let g = Arc::new(MultiplicativeGroup::from_finite_field(&f));
            let a_rand = (Integer::from(rng.bits(32)).modulo(&(&p-&Integer::ONE.clone()).complete())) + 1 ;
            let a = Element::new(g.clone(), a_rand);
            let a_inv = a.mul_inv();
            assert_eq!(Integer::ONE, a_inv.mul_ref(&a).get_rep(), "We have a: {}, a_inv: {}, p: {}", a.get_rep(), a_inv.get_rep(), p);


            prime.next_prime_mut();

        }
    }


    #[test]
    fn test_exponentiation() {
        let n: u32 = 200;
        let mut rng = RandState::new();
        let mut prime: Integer = Integer::from(2);


        for _ in 2..n {
            let p = prime.clone();
            let f = Arc::new(FiniteField::new(p.clone()).unwrap());

            let a_rand = Integer::from(rng.bits(32));
            let a = Element::new(f.clone(), a_rand.clone());

            let x_rand = Integer::from(rng.bits(32));

            // let a_exp_check = a_rand.exp_residue(&x_rand, &p);
            let a_exp_check = a_rand.pow_mod_ref(&x_rand, &p).unwrap().complete();
            let a_exp = a.pow(&x_rand);

            assert_eq!(a_exp.get_rep(), &a_exp_check);

            prime.next_prime_mut();
        }
    }
}
