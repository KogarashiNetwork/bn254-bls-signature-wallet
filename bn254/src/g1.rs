use core::ops::Neg;

use crate::fq::Fq;

pub(crate) const G1_GENERATOR_X: Fq = Fq::one();
pub(crate) const G1_GENERATOR_Y: Fq = Fq::to_mont_form([2, 0, 0, 0]);

#[derive(Clone, Copy, Debug)]
pub struct G1Affine {
    pub(crate) x: Fq,
    pub(crate) y: Fq,
    is_infinity: bool,
}

impl G1Affine {
    pub(crate) fn is_identity(self) -> bool {
        self.is_infinity
    }

    pub const fn generator() -> Self {
        Self {
            x: G1_GENERATOR_X,
            y: G1_GENERATOR_Y,
            is_infinity: false,
        }
    }
}

impl Neg for G1Affine {
    type Output = Self;

    fn neg(self) -> Self {
        Self {
            x: self.x,
            y: -self.y,
            is_infinity: self.is_infinity,
        }
    }
}
