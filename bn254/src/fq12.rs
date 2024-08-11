use core::ops::{Mul, MulAssign};

use crate::fq::Fq;
use crate::fq2::Fq2;
use crate::fq6::Fq6;
use crate::g1::G1Affine;
use crate::g2::PairingCoeff;
use crate::gt::Gt;
use crate::params::{BN_X, FROBENIUS_COEFF_FQ12_C1};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Fq12([Fq6; 2]);

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

    pub fn invert(self) -> Option<Self> {
        (self.0[0].square() - self.0[1].square().mul_by_nonres())
            .invert()
            .map(|t| Self([self.0[0] * t, self.0[1] * -t]))
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

    pub(crate) fn mul_by_034(self, c0: Fq2, c3: Fq2, c4: Fq2) -> Self {
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

    pub fn conjugate(self) -> Self {
        Self([self.0[0], -self.0[1]])
    }

    pub fn frobenius_map(self) -> Self {
        let c0 = self.0[0].frobenius_map();
        let c1 =
            self.0[1].frobenius_map() * Fq6([FROBENIUS_COEFF_FQ12_C1[1], Fq2::zero(), Fq2::zero()]);

        Self([c0, c1])
    }

    fn frobenius_maps(self, power: usize) -> Self {
        let c0 = self.0[0].frobenius_maps(power);
        let c1 = self.0[1].frobenius_maps(power);
        let c1 = Fq6([
            c1.0[0] * FROBENIUS_COEFF_FQ12_C1[power % 12],
            c1.0[1] * FROBENIUS_COEFF_FQ12_C1[power % 12],
            c1.0[2] * FROBENIUS_COEFF_FQ12_C1[power % 12],
        ]);

        Self([c0, c1])
    }

    pub fn final_exp(self) -> Gt {
        fn fp4_square(a: Fq2, b: Fq2) -> (Fq2, Fq2) {
            let t0 = a.square();
            let t1 = b.square();
            let mut t2 = t1.mul_by_nonres();
            let c0 = t2 + t0;
            t2 = a + b;
            t2 = t2.square();
            t2 -= t0;
            let c1 = t2 - t1;

            (c0, c1)
        }
        // Adaptation of Algorithm 5.5.4, Guide to Pairing-Based Cryptography
        // Faster Squaring in the Cyclotomic Subgroup of Sixth Degree Extensions
        // https://eprint.iacr.org/2009/565.pdf
        #[must_use]
        fn cyclotomic_square(f: Fq12) -> Fq12 {
            let mut z0 = f.0[0].0[0];
            let mut z4 = f.0[0].0[1];
            let mut z3 = f.0[0].0[2];
            let mut z2 = f.0[1].0[0];
            let mut z1 = f.0[1].0[1];
            let mut z5 = f.0[1].0[2];

            let (t0, t1) = fp4_square(z0, z1);

            // For A
            z0 = t0 - z0;
            z0 = z0.double() + t0;

            z1 = t1 + z1;
            z1 = z1.double() + t1;

            let (mut t0, t1) = fp4_square(z2, z3);
            let (t2, t3) = fp4_square(z4, z5);

            // For C
            z4 = t0 - z4;
            z4 = z4.double() + t0;

            z5 = t1 + z5;
            z5 = z5.double() + t1;

            // For B
            t0 = t3.mul_by_nonres();
            z2 = t0 + z2;
            z2 = z2.double() + t0;

            z3 = t2 - z3;
            z3 = z3.double() + t2;

            Fq12([Fq6([z0, z4, z3]), Fq6([z2, z1, z5])])
        }

        #[must_use]
        fn cycolotomic_exp(f: Fq12) -> Fq12 {
            let mut res = Fq12::one();
            for is_one in (0..64).rev().map(|i| ((BN_X >> i) & 1) == 1) {
                res = cyclotomic_square(res);
                if is_one {
                    res *= f;
                }
            }
            res
        }

        let f = self;
        let f1 = f.conjugate();
        Gt(f.invert()
            .map(|mut f2| {
                f2 *= f1;
                let r = f2.frobenius_maps(2) * f2;

                let fp = r.frobenius_maps(1);
                let fp2 = r.frobenius_maps(2);
                let fp3 = fp2.frobenius_maps(1);

                let fu = cycolotomic_exp(r);
                let fu2 = cycolotomic_exp(fu);
                let fu3 = cycolotomic_exp(fu2);

                let y3 = fu.frobenius_maps(1).conjugate();

                let fu2p = fu2.frobenius_maps(1);
                let fu3p = fu3.frobenius_maps(1);

                let y2 = fu2.frobenius_maps(2);

                let y0 = fp * fp2 * fp3;
                let y1 = r.conjugate();
                let y5 = fu2.conjugate();

                let y4 = (fu * fu2p).conjugate();

                let mut y6 = cyclotomic_square((fu3 * fu3p).conjugate()) * y4 * y5;

                let mut t1 = y3 * y5 * y6;
                y6 *= y2;
                t1 = cyclotomic_square(cyclotomic_square(t1) * y6);

                let mut t0 = t1 * y1;
                t1 *= y0;
                t0 = cyclotomic_square(t0) * t1;
                t0
            })
            .unwrap())
    }

    pub(crate) const fn generator() -> Self {
        Fq12([
            Fq6([
                Fq2([
                    Fq([
                        0xc556f62b2a98671d,
                        0x23a59ac167bcf363,
                        0x5ef208445f5f6f37,
                        0x12adf27ccb29382a,
                    ]),
                    Fq([
                        0x2e02a64acbd60549,
                        0xd618018ea58e4add,
                        0x14d585f1a45ba647,
                        0x1832226987c434fc,
                    ]),
                ]),
                Fq2([
                    Fq([
                        0x2306e4312363b991,
                        0x465f6072d4023bf4,
                        0xa2ff062a4a77e736,
                        0x76ea6f18435864a,
                    ]),
                    Fq([
                        0x172d1f257a4d598e,
                        0xddf5bc7b7ffb5ac0,
                        0xae0b22c0bbb0f602,
                        0x1B158F3C2FAE9B18,
                    ]),
                ]),
                Fq2([
                    Fq([
                        0x5cf9cc917da86724,
                        0xc799dc487a0b2753,
                        0xdf2027bf1de17a7,
                        0x197cda6cc3e20636,
                    ]),
                    Fq([
                        0xf16c96d081754cdb,
                        0xce0394312bceeb55,
                        0x644e4dcf1f01ff0a,
                        0xcbea85ee0b236cc,
                    ]),
                ]),
            ]),
            Fq6([
                Fq2([
                    Fq([
                        0x1bb0ce0def1b82a1,
                        0x4c4c9fe1cadefa95,
                        0x746d9990cb12b27e,
                        0x13495c08e5d415c5,
                    ]),
                    Fq([
                        0x9458abcb56d24998,
                        0xb17540bd2a9e5adb,
                        0x9a9983c82e401a9f,
                        0x1614817a84c16291,
                    ]),
                ]),
                Fq2([
                    Fq([
                        0x8975b68a2bab1f9c,
                        0x2fdd826b796e0f35,
                        0x6a90a35fa03dfaa5,
                        0x1ffef4581607fc37,
                    ]),
                    Fq([
                        0x7002907c28ebfe11,
                        0x7b0591d3d080da67,
                        0xde7e5aa2181f138e,
                        0x210e437dfc43d951,
                    ]),
                ]),
                Fq2([
                    Fq([
                        0x988ae2485b36cf53,
                        0x5091cc0581334e54,
                        0xda7903229312ca0f,
                        0x2a2341538eaee95c,
                    ]),
                    Fq([
                        0xd34bab373157aa84,
                        0x3511ed44fd0d8598,
                        0x67e42a0bc2ced972,
                        0x2b8f1d5dfd20c55b,
                    ]),
                ]),
            ]),
        ])
    }
}

impl Mul<Fq12> for Fq12 {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        let aa = self.0[0] * rhs.0[0];
        let bb = self.0[1] * rhs.0[1];
        let o = rhs.0[0] + rhs.0[1];
        let c1 = self.0[1] + self.0[0];
        let c1 = c1 * o;
        let c1 = c1 - aa;
        let c1 = c1 - bb;
        let c0 = bb.mul_by_nonres();
        let c0 = c0 + aa;

        Self([c0, c1])
    }
}

impl MulAssign for Fq12 {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}
