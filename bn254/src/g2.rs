use crate::fq2::Fq2;

pub(crate) struct G2Affine {
    x: Fq2,
    y: Fq2,
    is_infinity: bool,
}
