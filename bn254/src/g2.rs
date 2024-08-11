use crate::fq2::Fq2;
use crate::pairing::{SIX_U_PLUS_2_NAF, XI_TO_Q_MINUS_1_OVER_2};
use crate::params::FROBENIUS_COEFF_FQ6_C1;

use core::ops::Neg;

#[derive(Clone, Copy, Debug)]
pub struct G2Affine {
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

impl G2Projective {
    pub(crate) fn double_eval(&mut self) -> PairingCoeff {
        // Adaptation of Algorithm 26, https://eprint.iacr.org/2010/354.pdf
        let tmp0 = self.x.square();
        let tmp1 = self.y.square();
        let tmp2 = tmp1.square();
        let tmp3 = (tmp1 + self.x).square() - tmp0 - tmp2;
        let tmp3 = tmp3.double();
        let tmp4 = tmp0.double() + tmp0;
        let tmp6 = self.x + tmp4;
        let tmp5 = tmp4.square();
        let zsquared = self.z.square();
        self.x = tmp5 - tmp3.double();
        self.z = (self.z + self.y).square() - tmp1 - zsquared;
        self.y = (tmp3 - self.x) * tmp4 - tmp2.double().double().double();
        let tmp3 = -(tmp4 * zsquared).double();
        let tmp6 = tmp6.square() - tmp0 - tmp5;
        let tmp1 = tmp1.double().double();
        let tmp6 = tmp6 - tmp1;
        let tmp0 = self.z * zsquared;
        let tmp0 = tmp0.double();

        PairingCoeff(tmp0, tmp3, tmp6)
    }

    pub(crate) fn add_eval(&mut self, rhs: G2Affine) -> PairingCoeff {
        // Adaptation of Algorithm 27, https://eprint.iacr.org/2010/354.pdf
        let zsquared = self.z.square();
        let ysquared = rhs.y.square();
        let t0 = zsquared * rhs.x;
        let t1 = ((rhs.y + self.z).square() - ysquared - zsquared) * zsquared;
        let t2 = t0 - self.x;
        let t3 = t2.square();
        let t4 = t3.double().double();
        let t5 = t4 * t2;
        let t6 = t1 - self.y.double();
        let t9 = t6 * rhs.x;
        let t7 = t4 * self.x;
        self.x = t6.square() - t5 - t7.double();
        self.z = (self.z + t2).square() - zsquared - t3;
        let t10 = rhs.y + self.z;
        let t8 = (t7 - self.x) * t6;
        let t0 = self.y * t5;
        self.y = t8 - t0.double();
        let t10 = t10.square() - ysquared;
        let ztsquared = self.z.square();
        let t10 = t10 - ztsquared;
        let t9 = t9.double() - t10;
        let t10 = self.z.double();
        let t1 = -t6.double();

        PairingCoeff(t10, t1, t9)
    }
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

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct PairingCoeff(pub(crate) Fq2, pub(crate) Fq2, pub(crate) Fq2);

pub struct G2PairingAffine {
    pub(crate) coeffs: Vec<PairingCoeff>,
    is_infinity: bool,
}

impl G2PairingAffine {
    pub fn is_identity(&self) -> bool {
        self.is_infinity
    }
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

            for i in (1..SIX_U_PLUS_2_NAF.len()).rev() {
                coeffs.push(g2_projective.double_eval());
                let x = SIX_U_PLUS_2_NAF[i - 1];
                match x {
                    1 => {
                        coeffs.push(g2_projective.add_eval(g2));
                    }
                    -1 => {
                        coeffs.push(g2_projective.add_eval(neg));
                    }
                    _ => continue,
                }
            }

            let mut q = g2;

            q.x.0[1] = -q.x.0[1];
            q.x *= FROBENIUS_COEFF_FQ6_C1[1];

            q.y.0[1] = -q.y.0[1];
            q.y *= XI_TO_Q_MINUS_1_OVER_2;

            coeffs.push(g2_projective.add_eval(q));

            let mut minusq2 = g2;
            minusq2.x *= FROBENIUS_COEFF_FQ6_C1[2];

            coeffs.push(g2_projective.add_eval(minusq2));

            Self {
                coeffs,
                is_infinity: false,
            }
        }
    }
}
