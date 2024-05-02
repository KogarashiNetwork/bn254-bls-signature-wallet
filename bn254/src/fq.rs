use core::ops::Neg;

use crate::limbs::neg;

pub(crate) const MODULUS: [u64; 4] = [
    0x3c208c16d87cfd47,
    0x97816a916871ca8d,
    0xb85045b68181585d,
    0x30644e72e131a029,
];

pub(crate) const R: [u64; 4] = [
    0xd35d438dc58f0d9d,
    0x0a78eb28f5c70b3d,
    0x666ea36f7879462c,
    0x0e0a77c19a07df2f,
];

#[derive(Clone, Copy, Debug)]
pub(crate) struct Fq([u64; 4]);

impl Fq {
    pub(crate) fn zero() -> Self {
        Self([0; 4])
    }

    pub(crate) fn one() -> Self {
        Self(R)
    }
}

impl Neg for Fq {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self(neg(self.0, MODULUS))
    }
}
