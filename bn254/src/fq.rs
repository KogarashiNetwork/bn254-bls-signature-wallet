use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub};

use crate::limbs::{add, double, invert, little_fermat, mont, mul, neg, square, sub};

pub(crate) const MODULUS: [u64; 4] = [
    0x3c208c16d87cfd47,
    0x97816a916871ca8d,
    0xb85045b68181585d,
    0x30644e72e131a029,
];

/// R = 2^256 mod q
pub(crate) const R: [u64; 4] = [
    0xd35d438dc58f0d9d,
    0x0a78eb28f5c70b3d,
    0x666ea36f7879462c,
    0x0e0a77c19a07df2f,
];

/// INV = -(q^{-1} mod 2^64) mod 2^64
pub(crate) const INV: u64 = 0x87d20782e4866389;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Fq(pub(crate) [u64; 4]);

impl Fq {
    pub(crate) const fn zero() -> Self {
        Self([0; 4])
    }

    pub(crate) const fn is_zero(self) -> bool {
        self.0[0] == 0 && self.0[1] == 0 && self.0[2] == 0 && self.0[3] == 0
    }

    pub(crate) const fn one() -> Self {
        Self(R)
    }

    pub(crate) const fn double(self) -> Self {
        Self(double(self.0, MODULUS))
    }

    pub(crate) const fn square(self) -> Self {
        Self(square(self.0, MODULUS, INV))
    }

    pub(crate) fn invert(self) -> Option<Self> {
        match invert(self.0, little_fermat(MODULUS), R, MODULUS, INV) {
            Some(x) => Some(Self(x)),
            None => None,
        }
    }

    pub(crate) const fn to_mont_form(val: [u64; 4]) -> Self {
        Self(mont(
            [val[0], val[1], val[2], val[3], 0, 0, 0, 0],
            MODULUS,
            INV,
        ))
    }
}

impl Add for Fq {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self(add(self.0, rhs.0, MODULUS))
    }
}

impl AddAssign for Fq {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl Sub for Fq {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        Self(sub(self.0, rhs.0, MODULUS))
    }
}

impl Neg for Fq {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self(neg(self.0, MODULUS))
    }
}

impl Mul<Fq> for Fq {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        Self(mul(self.0, rhs.0, MODULUS, INV))
    }
}

impl MulAssign for Fq {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}
