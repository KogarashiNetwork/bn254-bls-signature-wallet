use crate::fq::Fq;
use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Fq2(pub(crate) [Fq; 2]);

impl Fq2 {
    pub(crate) const fn zero() -> Self {
        Self([Fq::zero(); 2])
    }

    pub(crate) const fn one() -> Self {
        Self([Fq::one(), Fq::zero()])
    }

    pub(crate) fn double(self) -> Self {
        Self([self.0[0].double(), self.0[1].double()])
    }

    pub(crate) fn square(self) -> Self {
        let re = self.0[0].square() - self.0[1].square();
        let im = (self.0[0] * self.0[1]).double();
        Self([re, im])
    }

    /// Multiply this element by quadratic nonresidue 9 + u.
    pub(crate) fn mul_by_nonres(self) -> Self {
        // (xi+y)(i+9) = (9x+y)i+(9y-x)
        let t0 = self.0[0];
        let t1 = self.0[1];
        // 8*x*i + 8*y
        let mut res = self.double().double().double();

        // 9*y - x
        res.0[0] += t0 - t1;
        // (9*x + y)i
        res.0[1] += t0 + t1;
        res
    }
}

impl Add for Fq2 {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self([self.0[0] + rhs.0[0], self.0[1] + rhs.0[1]])
    }
}

impl AddAssign for Fq2 {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl Neg for Fq2 {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self {
        Self([-self.0[0], -self.0[1]])
    }
}

impl Sub for Fq2 {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self([self.0[0] - rhs.0[0], self.0[1] - rhs.0[1]])
    }
}

impl SubAssign for Fq2 {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl Mul<Fq2> for Fq2 {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        let re = (self.0[0] * rhs.0[0]) - (self.0[1] * rhs.0[1]);
        let im = (self.0[0] * rhs.0[1]) + (self.0[1] * rhs.0[0]);
        Self([re, im])
    }
}

impl MulAssign for Fq2 {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}
