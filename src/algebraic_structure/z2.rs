use std::ops::{Add, AddAssign, Div, Mul, MulAssign, Not, SubAssign};
use num::traits::{Zero, One};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Z2(pub bool);

impl Add for Z2 {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Z2(self.0 ^ other.0)  // XOR for addition
    }
}

impl Add<&Z2> for Z2 {
    type Output = Self;
    fn add(self, other: &Self) -> Self {
        Z2(self.0 ^ other.0)
    }
}

impl Add<Z2> for &Z2 {
    type Output = Z2;
    fn add(self, other: Z2) -> Z2 {
        Z2(self.0 ^ other.0)
    }
}

impl Add<&Z2> for &Z2 {
    type Output = Z2;
    fn add(self, other: &Z2) -> Z2 {
        Z2(self.0 ^ other.0)
    }
}

impl AddAssign for Z2 {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl SubAssign for Z2 {
    fn sub_assign(&mut self, rhs: Self) {
        *self += rhs;
    }
}

impl Div for Z2 {
    type Output = Self;
    fn div(self, rhs: Self) -> Self::Output {
        assert!(rhs.0, "Division by zero is undefined in Z_2");
        self
    }
}

impl Mul for Z2 {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        Z2(self.0 & other.0)  // AND for multiplication
    }
}

impl Mul<&Z2> for Z2 {
    type Output = Self;
    fn mul(self, other: &Self) -> Self {
        Z2(self.0 & other.0)  // AND for multiplication
    }
}

impl Mul<Z2> for &Z2 {
    type Output = Z2;
    fn mul(self, other: Z2) -> Z2 {
        Z2(self.0 & other.0)  // AND for multiplication
    }
}

impl Mul<&Z2> for &Z2 {
    type Output = Z2;
    fn mul(self, other: &Z2) -> Z2 {
        Z2(self.0 & other.0)  // AND for multiplication
    }
}

impl MulAssign for Z2 {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

// Implement Zero and One traits for ndarray compatibility
impl Zero for Z2 {
    fn zero() -> Self {
        Z2(false)
    }
    fn is_zero(&self) -> bool {
        !self.0
    }
}

impl One for Z2 {
    fn one() -> Self {
        Z2(true)
    }
}

// Implement Not since ndarray might expect it
impl Not for Z2 {
    type Output = Self;
    fn not(self) -> Self {
        Z2(!self.0)
    }
}

// Display implementation for debugging
impl std::fmt::Display for Z2 {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", if self.0 { "1" } else { "0" })
    }
}
