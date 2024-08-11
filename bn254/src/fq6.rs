use crate::fq2::Fq2;
use crate::params::{FROBENIUS_COEFF_FQ6_C1, FROBENIUS_COEFF_FQ6_C2};
use core::ops::{Add, Mul, Sub, SubAssign};
use std::ops::Neg;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Fq6(pub(crate) [Fq2; 3]);

impl Fq6 {
    pub const fn zero() -> Self {
        Self([Fq2::zero(); 3])
    }

    pub const fn one() -> Self {
        Self([Fq2::one(), Fq2::zero(), Fq2::zero()])
    }

    pub(crate) fn double(self) -> Self {
        Self([self.0[0].double(), self.0[1].double(), self.0[2].double()])
    }

    pub fn square(self) -> Self {
        let s0 = self.0[0].square();
        let ab = self.0[0] * self.0[1];
        let s1 = ab.double();
        let mut s2 = self.0[0];
        s2 -= self.0[1];
        s2 += self.0[2];
        s2 = s2.square();
        let bc = self.0[1] * self.0[2];
        let s3 = bc.double();
        let s4 = self.0[2].square();

        let c0 = s3.mul_by_nonres() + s0;
        let c1 = s4.mul_by_nonres() + s1;
        let c2 = s1 + s2 + s3 - s0 - s4;

        Self([c0, c1, c2])
    }

    pub(crate) fn invert(self) -> Option<Self> {
        let c0 = (self.0[1] * self.0[2]).mul_by_nonres();
        let c0 = self.0[0].square() - c0;

        let c1 = self.0[2].square().mul_by_nonres();
        let c1 = c1 - (self.0[0] * self.0[1]);

        let c2 = self.0[1].square();
        let c2 = c2 - (self.0[0] * self.0[2]);

        let tmp = ((self.0[1] * c2) + (self.0[2] * c1)).mul_by_nonres();
        let tmp = tmp + (self.0[0] * c0);

        tmp.invert().map(|t| Self([t * c0, t * c1, t * c2]))
    }

    pub(crate) fn mul_by_01(&self, c0: Fq2, c1: Fq2) -> Self {
        let a_a = self.0[0] * c0;
        let b_b = self.0[1] * c1;
        let t1 = ((self.0[1] + self.0[2]) * c1 - b_b).mul_by_nonres() + a_a;
        let t2 = (c0 + c1) * (self.0[0] + self.0[1]) - a_a - b_b;
        let t3 = (self.0[0] + self.0[2]) * c0 - a_a + b_b;

        Self([t1, t2, t3])
    }

    pub(crate) fn mul_by_nonres(self) -> Self {
        Self([self.0[2].mul_by_nonres(), self.0[0], self.0[1]])
    }

    pub(crate) fn frobenius_map(&self) -> Self {
        let c0 = self.0[0].frobenius_map();
        let c1 = self.0[1].frobenius_map() * FROBENIUS_COEFF_FQ6_C1[1];
        let c2 = self.0[2].frobenius_map() * FROBENIUS_COEFF_FQ6_C2[1];

        Fq6([c0, c1, c2])
    }

    pub(crate) fn frobenius_maps(self, power: usize) -> Self {
        let c0 = self.0[0].frobenius_maps(power);
        let c1 = self.0[1].frobenius_maps(power) * FROBENIUS_COEFF_FQ6_C1[power % 6];
        let c2 = self.0[2].frobenius_maps(power) * FROBENIUS_COEFF_FQ6_C2[power % 6];

        Self([c0, c1, c2])
    }
}

impl Add for Fq6 {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Self([
            self.0[0] + rhs.0[0],
            self.0[1] + rhs.0[1],
            self.0[2] + rhs.0[2],
        ])
    }
}

impl Neg for Fq6 {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self([-self.0[0], -self.0[1], -self.0[2]])
    }
}

impl Sub for Fq6 {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        Self([
            self.0[0] - rhs.0[0],
            self.0[1] - rhs.0[1],
            self.0[2] - rhs.0[2],
        ])
    }
}

impl SubAssign for Fq6 {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl Mul<Fq6> for Fq6 {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        let a_a = self.0[0] * rhs.0[0];
        let b_b = self.0[1] * rhs.0[1];
        let c_c = self.0[2] * rhs.0[2];

        let mut t1 = rhs.0[1] + rhs.0[2];
        {
            let tmp = self.0[1] + self.0[2];

            t1 *= tmp;
            t1 -= b_b;
            t1 -= c_c;
            t1 = t1.mul_by_nonres();
            t1 += a_a;
        }

        let mut t3 = rhs.0[0] + rhs.0[2];
        {
            let tmp = self.0[0] + self.0[2];

            t3 *= tmp;
            t3 -= a_a;
            t3 += b_b;
            t3 -= c_c;
        }

        let mut t2 = rhs.0[0] + rhs.0[1];
        {
            let tmp = self.0[0] + self.0[1];

            t2 *= tmp;
            t2 -= a_a;
            t2 -= b_b;
            t2 += c_c.mul_by_nonres();
        }

        Self([t1, t2, t3])
    }
}
