use rug::Integer;
use std::sync::Arc;
pub mod finite_field;
pub mod z2;


pub trait HasRepresentation {
    fn make_representation(&self, repr: Integer) -> Integer;
}


pub trait HasMul: HasRepresentation + Clone {
    fn mul(&self, a: &Element<Self>, b: &Element<Self>) -> Element<Self>;
    fn pow(&self, a: &Element<Self>, b: &Integer) -> Element<Self>;
}


pub trait HasAdd: HasRepresentation + Clone {
    fn add(&self, a: &Element<Self>, b: &Element<Self>) -> Element<Self>;
}


pub trait HasSub: HasAdd {
    fn add_inv(&self, a: &Element<Self>) -> Element<Self>;
    fn sub_ref(&self, a: &Element<Self>, b: &Element<Self>) -> Element<Self> {
        self.add(a, &self.add_inv(b))
    }
}


pub trait HasDiv: HasMul {
    fn mul_inv(&self, a: &Element<Self>) -> Element<Self>;
    fn div(&self, a: &Element<Self>, b: &Element<Self>) -> Element<Self> {
        self.mul(a, &self.mul_inv(&b))
    }
}


#[derive(Debug, Clone)]
pub struct Element<T: HasRepresentation + Clone> {
    outer_structure: Arc<T>,
    representation: Integer,
}


impl<T: HasRepresentation + Clone> Element<T> {
    pub fn new(outer_structure: Arc<T>, repr: Integer) -> Element<T> {
        let representation: Integer = outer_structure.make_representation(repr);
        Element { outer_structure, representation }
    }


    pub fn get_outer_structure(&self) -> Arc<T> {
        self.outer_structure.clone()
    }


    pub fn get_rep(&self) -> &Integer {
        &self.representation
    }
}


impl<T: HasMul> Element<T> {
    pub fn mul_ref(&self, _rhs: &Element<T>) -> Element<T> {
        self.outer_structure.mul(&self, _rhs)
    }

    pub fn pow(&self, a: &Integer) -> Element<T> {
        self.get_outer_structure().pow(self, a)
    }
}


impl<T: HasAdd> Element<T> {
    pub fn add_ref(&self, _rhs: &Element<T>) -> Element<T> {
        self.outer_structure.add(&self, _rhs)
    }
}


impl<T: HasDiv> Element<T> {
    pub fn mul_inv(&self) -> Element<T> {
        self.outer_structure.mul_inv(&self)
    }

    pub fn div_ref(&self, b: &Element<T>) -> Element<T> {
        self.outer_structure.div(&self, b)
    }
}


impl<T: HasSub> Element<T> {
    pub fn add_inv(&self) -> Element<T> {
        self.outer_structure.add_inv(&self)
    }

    pub fn sub_ref(&self, b: &Element<T>) -> Element<T> {
        self.outer_structure.sub_ref(&self, b)
    }
}
