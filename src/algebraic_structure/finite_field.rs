use crate::algebraic_structure::{Element, HasAdd, HasMul, HasRepresentation, HasSub};
use crate::integers::integer_computations::{extended_euclidean_ordered, extended_euclidean_to_integers, pow_rug};
use rug::ops::SubFrom;
use rug::{integer::IsPrime, Complete, Integer};
use std::sync::Arc;


use super::HasDiv;

#[derive(Debug, Clone)]
pub struct FiniteField {
    // This struct will only consider finite fields isomorphic to Z_p for p prime.
    size: Integer,
}


impl HasRepresentation for FiniteField {
    fn make_representation(&self, repr: Integer) -> Integer {
        repr.modulo(self.mod_num())
    }
}


impl HasMul for FiniteField {
    fn mul(&self, a: &Element<FiniteField>, b: &Element<FiniteField>) -> Element<FiniteField> {
        Element::new(
            a.get_outer_structure(),
            (a.get_rep() * b.get_rep()).complete() % self.mod_num()
        )
    }

    fn pow(&self, a: &Element<Self>, b: &Integer) -> Element<Self> {
        Element::new(
            a.get_outer_structure(),
            pow_rug(a.get_rep(), b, self.mod_num())
        )
    }
}


impl HasAdd for FiniteField {
    fn add(&self, a: &Element<FiniteField>, b: &Element<FiniteField>) -> Element<FiniteField> {
        Element::new(
            a.get_outer_structure().clone(),
            (a.get_rep() + b.get_rep()).complete() % self.mod_num()
        )
    }
}


impl HasDiv for FiniteField {
    fn mul_inv(&self, a: &Element<Self>) -> Element<Self> {
        if a.get_rep().is_zero() {
            panic!("Zero Division");
        }
        let (_, _, y) = extended_euclidean_ordered(self.mod_num(), a.get_rep());
        Element::new(
            a.get_outer_structure(),
            y
        )
    }
}


impl HasSub for FiniteField {
    fn add_inv(&self, a: &Element<Self>) -> Element<Self> {
        let mut representation: Integer = a.get_rep().clone();
        representation.sub_from(self.mod_num());
        Element {
            outer_structure: a.get_outer_structure().clone(),
            representation,
        }
    }
}


impl FiniteField {
    pub fn new(size: Integer) -> Option<FiniteField> {
        if size.is_probably_prime(30) != IsPrime::No {
            Some(FiniteField {
                size,
            })
        } else {
            None
        }
    }

    pub fn one(self) -> Element<FiniteField> {
        Element {
            outer_structure: Arc::new(self),
            representation: Integer::ONE.clone(),
        }
    }

    pub fn zero(self) -> Element<FiniteField> {
        Element {
            outer_structure: Arc::new(self),
            representation: Integer::ZERO.clone(),
        }
    }

    pub fn get_size(&self) -> Integer {
        self.size.clone()
    }

    fn mod_num(&self) -> &Integer {
        &self.size
    }


    pub fn extended_euclidean(&self, a: &Element<FiniteField>, b: &Element<FiniteField>) -> (Element<FiniteField>, Element<FiniteField>, Element<FiniteField>) {
        let (d, x, y) = extended_euclidean_to_integers(a, b);
        (
            Element::new(
                a.get_outer_structure(),
                d
            ),
            Element::new(
                a.get_outer_structure(),
                x
            ),
            Element::new(
                a.get_outer_structure(),
                y
            ),
        )
    }
}


#[derive(Debug, Clone)]
pub struct MultiplicativeGroup {
    mod_num: Integer,
}


impl HasRepresentation for MultiplicativeGroup {
    // Makes representation for creating elements.
    // As 0 is not present in the mutliplicative group, it is assumed you meant identity and will
    // get 1 instead.
    fn make_representation(&self, repr: Integer) -> Integer {
        let representation = repr.modulo(self.mod_num());
        if representation == 0 {
            Integer::ONE.clone()
        } else {
            representation
        }
    }
}


impl HasMul for MultiplicativeGroup {
    fn mul(&self, a: &Element<MultiplicativeGroup>, b: &Element<MultiplicativeGroup>) -> Element<MultiplicativeGroup> {
        Element::new(
            a.get_outer_structure(),
            (a.get_rep() * b.get_rep()).complete() % self.mod_num()
        )
    }

    fn pow(&self, a: &Element<Self>, b: &Integer) -> Element<Self> {
        Element::new(
            a.get_outer_structure(),
            pow_rug(a.get_rep(), b, self.mod_num()),
        )
    }
}


impl MultiplicativeGroup {
    pub fn new(mod_num: Integer) -> MultiplicativeGroup {
        MultiplicativeGroup {
            mod_num,
        }
    }


    pub fn from_finite_field(finite_field: &FiniteField) -> MultiplicativeGroup {
        let mod_num = finite_field.mod_num().clone();
        MultiplicativeGroup {
            mod_num,
        }
    }


    pub fn mod_num(&self) -> &Integer {
        &self.mod_num
    }


    pub fn get_size(&self) -> Integer {
        (self.mod_num() - Integer::ONE).complete()
    }
}


impl HasDiv for MultiplicativeGroup {
    fn mul_inv(&self, a: &Element<Self>) -> Element<Self> {
        let (_, _, y) = extended_euclidean_ordered(self.mod_num(), a.get_rep());
        Element::new(
            a.get_outer_structure(),
            y
        )
    }
}
