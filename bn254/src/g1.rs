use crate::fq::Fq;

#[derive(Clone, Copy, Debug)]
pub(crate) struct G1Affine {
    pub(crate) x: Fq,
    pub(crate) y: Fq,
    is_infinity: bool,
}

impl G1Affine {
    pub(crate) fn is_identity(self) -> bool {
        self.is_infinity
    }
}
