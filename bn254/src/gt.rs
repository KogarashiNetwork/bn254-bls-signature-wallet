use core::ops::{Add, Neg};

use crate::fq12::Fq12;

#[derive(Debug, PartialEq, Eq)]
pub struct Gt(pub Fq12);

impl Gt {
    pub fn identity() -> Self {
        Self(Fq12::one())
    }

    pub fn generator() -> Self {
        Self(Fq12::generator())
    }
}

impl Add for Gt {
    type Output = Gt;

    fn add(self, rhs: Gt) -> Gt {
        Self(self.0 * rhs.0)
    }
}

impl Neg for Gt {
    type Output = Gt;

    fn neg(self) -> Gt {
        Gt(self.0.conjugate())
    }
}
