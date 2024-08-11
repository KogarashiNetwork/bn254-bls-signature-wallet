use crate::fq2::Fq2;
use crate::fq6::Fq6;
use crate::g1::G1Affine;
use crate::g2::PairingCoeff;

#[derive(Clone, Copy, Debug)]
pub(crate) struct Fq12([Fq6; 2]);

impl Fq12 {
    pub const fn one() -> Self {
        Self([Fq6::one(), Fq6::zero()])
    }

    pub fn square(self) -> Self {
        let re = self.0[0].square() - self.0[1].square();
        let im = (self.0[0] * self.0[1]).double();
        Self([re, im])
    }

    pub fn square_assign(&mut self) {
        *self = self.square()
    }

    // twisting isomorphism from E to E'
    pub(crate) fn untwist(self, coeffs: PairingCoeff, g1: G1Affine) -> Self {
        let mut c0 = coeffs.0;
        let mut c1 = coeffs.1;

        c0.0[0] *= g1.y;
        c0.0[1] *= g1.y;

        c1.0[0] *= g1.x;
        c1.0[1] *= g1.x;

        self.mul_by_034(c0, c1, coeffs.2)
    }

    pub fn mul_by_034(self, c0: Fq2, c3: Fq2, c4: Fq2) -> Self {
        let t0 = Fq6([
            self.0[0].0[0] * c0,
            self.0[0].0[1] * c0,
            self.0[0].0[2] * c0,
        ]);
        let mut t1 = self.0[1];
        t1 = t1.mul_by_01(c3, c4);
        let o = c0 + c3;
        let mut t2 = self.0[0] + self.0[1];
        t2 = t2.mul_by_01(o, c4);
        t2 -= t0;
        let b = t2 - t1;
        t1 = t1.mul_by_nonres();
        let a = t0 + t1;
        Self([a, b])
    }
}
