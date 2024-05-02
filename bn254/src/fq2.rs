use crate::fq::Fq;

pub(crate) struct Fq2([Fq; 2]);

impl Fq2 {
    pub(crate) fn zero() -> Self {
        Self([Fq::zero(); 2])
    }

    pub(crate) fn one() -> Self {
        Self([Fq::one(), Fq::zero()])
    }
}
