use crate::fq2::Fq2;

use core::ops::Neg;

pub(crate) struct G2Affine {
    x: Fq2,
    y: Fq2,
    is_infinity: bool,
}

impl G2Affine {
    fn is_identity(self) -> bool {
        self.is_infinity
    }
}

impl Neg for G2Affine {
    type Output = Self;

    fn neg(self) -> Self {
        Self {
            x: self.x,
            y: -self.y,
            is_infinity: self.is_infinity,
        }
    }
}

struct G2Projective {
    x: Fq2,
    y: Fq2,
    z: Fq2,
}

impl From<G2Affine> for G2Projective {
    fn from(affine: G2Affine) -> Self {
        if affine.is_identity() {
            Self {
                x: Fq2::zero(),
                y: Fq2::one(),
                z: Fq2::zero(),
            }
        } else {
            Self {
                x: affine.x,
                y: affine.y,
                z: Fq2::one(),
            }
        }
    }
}

struct PairingCoeff(Fq2, Fq2, Fq2);

pub(crate) struct G2PairingAffine {
    pub coeffs: Vec<PairingCoeff>,
    is_infinity: bool,
}

impl From<G2Affine> for G2PairingAffine {
    fn from(g2: G2Affine) -> Self {
        if g2.is_identity() {
            Self {
                coeffs: vec![],
                is_infinity: true,
            }
        } else {
            let mut coeffs = vec![];
            let mut g2_projective = G2Projective::from(g2);
            let neg = -g2;

            Self {
                coeffs,
                is_infinity: false,
            }
        }
    }
}
