use crate::math::{sba, sbb};

pub(crate) const fn neg(a: [u64; 4], p: [u64; 4]) -> [u64; 4] {
    if (a[0] | a[1] | a[2] | a[3]) == 0 {
        a
    } else {
        let (l0, b) = sba(p[0], a[0]);
        let (l1, b) = sbb(p[1], a[1], b);
        let (l2, b) = sbb(p[2], a[2], b);
        let l3 = (p[3]).wrapping_sub(a[3]).wrapping_sub(b >> 63);

        [l0, l1, l2, l3]
    }
}
