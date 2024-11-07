use core::iter::Product;
use std::ops::{Mul, MulAssign};

use crate::{
    ConfigZZVecGoldilocks256, ConfigZZpX, ConfigZZpXGoldilocks256, NTTDomain, PolyGoldilock256, Polynomial, Goldilocks,
    PolynomialRing, ZZVec,
};

/// Ring over ZZ_q/(x^512+1)
pub type RingGoldilock256 = PolyGoldilock256;
/// Configuration for ring over ZZ_q/(x^512+1)
pub type ConfigRingGoldilocks256 = ConfigZZpXGoldilocks256;

// ========================
// multiplications
// ========================
impl Mul for RingGoldilock256 {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        self.mul(&rhs)
    }
}

impl<'b> Mul<&'b RingGoldilock256> for RingGoldilock256 {
    type Output = RingGoldilock256;

    #[inline]
    fn mul(self, rhs: &'b RingGoldilock256) -> RingGoldilock256 {
        let mut res = self;
        res.mul_assign(rhs);
        res
    }
}

impl MulAssign for RingGoldilock256 {
    #[inline]
    fn mul_assign(&mut self, rhs: RingGoldilock256) {
        self.mul_assign(&rhs)
    }
}

impl<'b> MulAssign<&'b RingGoldilock256> for RingGoldilock256 {
    #[inline]
    fn mul_assign(&mut self, rhs: &'b RingGoldilock256) {
        let a: ZZVec<ConfigZZVecGoldilocks256> = NTTDomain::forward_ntt(self);
        let b: ZZVec<ConfigZZVecGoldilocks256> = NTTDomain::forward_ntt(rhs);
        println!("a: {}", a);
        println!("b: {}", b);
        let c = a * b;
        println!("c: {}", c);
        *self = c.reverse_ntt();
    }
}

impl<T> Product<T> for RingGoldilock256
where
    T: core::borrow::Borrow<Self>,
{
    fn product<I: Iterator<Item = T>>(iter: I) -> Self {
        iter.fold(Self::one(), |acc, item| acc * item.borrow())
    }
}

impl PolynomialRing<ConfigRingGoldilocks256, ConfigZZpXGoldilocks256> for RingGoldilock256 {}

impl RingGoldilock256 {
    /// school book multiplication
    /// output = a(x) * b(x) mod x^N +1 mod MODULUS
    /// using school-book multiplications
    pub fn schoolbook_mul(a: &Self, b: &Self) -> Self {
        use crate::ConfigZZp;
        use crate::ConfigZZpX;

        let a = &a.coeffs;
        let b = &b.coeffs;
        let modulus = <ConfigRingGoldilocks256 as ConfigZZpX>::BaseConfig::MODULUS;
        const N: usize = <ConfigRingGoldilocks256 as ConfigZZpX>::DIM;

        let mut buf = [0u128; N << 1];
        let mut c = [0; N];
        for i in 0..N {
            for j in 0..N {
                buf[i + j] += (a[i].0 as u128 * b[j].0 as u128) % modulus as u128;
            }
        }

        for i in 0..N {
            c[i] = ((buf[i] + modulus as u128 - (buf[i + N] % modulus as u128)) % modulus as u128)
                as u64;
        }
        Self::from_primitive_types(&c)
    }

    /// Dot product of two vectors, and vectors elements are polynomials ring elements
    pub fn dot_product(a: &Vec<Self>, b: &Vec<Self>) -> Self {
        let mut res = Self::zero();
        for (i, coeff) in a.iter().enumerate() {
            res += Self::schoolbook_mul(coeff, &b[i]);
        }
        res
    }
}


#[test]
fn test_ring_mul() {
    use rand::SeedableRng;
    use rand_xorshift::XorShiftRng;

    let mut rng = XorShiftRng::from_seed([
        0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc,
        0xe5,
    ]);
    let a = RingGoldilock256::random(&mut rng, None);
    let b = RingGoldilock256::random(&mut rng, None);
    let c = RingGoldilock256::schoolbook_mul(&a, &b);
    let d = a * b;
    assert_eq!(c, d)
}


#[test]
fn test_ring_dot_product_example() {
    // Define the polynomials for s1
    // s1 = (1 + 2x, 3 + 4x)
    let coeffs_s1_0 = [
        Goldilocks::from(1u64), // Coefficient for x^0
        Goldilocks::from(2u64), // Coefficient for x^1
    ];
    let mut s1_0_coeffs = vec![Goldilocks::from(0u64); ConfigZZpXGoldilocks256::DIM];
    s1_0_coeffs[..coeffs_s1_0.len()].copy_from_slice(&coeffs_s1_0);
    let s1_0 = RingGoldilock256::from_coefficients_vec_unchecked(s1_0_coeffs);

    let coeffs_s1_1 = [
        Goldilocks::from(3u64), // Coefficient for x^0
        Goldilocks::from(4u64), // Coefficient for x^1
    ];
    let mut s1_1_coeffs = vec![Goldilocks::from(0u64); ConfigZZpXGoldilocks256::DIM];
    s1_1_coeffs[..coeffs_s1_1.len()].copy_from_slice(&coeffs_s1_1);
    let s1_1 = RingGoldilock256::from_coefficients_vec_unchecked(s1_1_coeffs);

    let s1 = vec![s1_0, s1_1];

    // Define the polynomials for s2
    // s2 = (5 + 6x, 7 + 8x)
    let coeffs_s2_0 = [
        Goldilocks::from(5u64), // Coefficient for x^0
        Goldilocks::from(6u64), // Coefficient for x^1
    ];
    let mut s2_0_coeffs = vec![Goldilocks::from(0u64); ConfigZZpXGoldilocks256::DIM];
    s2_0_coeffs[..coeffs_s2_0.len()].copy_from_slice(&coeffs_s2_0);
    let s2_0 = RingGoldilock256::from_coefficients_vec_unchecked(s2_0_coeffs);

    let coeffs_s2_1 = [
        Goldilocks::from(7u64), // Coefficient for x^0
        Goldilocks::from(8u64), // Coefficient for x^1
    ];
    let mut s2_1_coeffs = vec![Goldilocks::from(0u64); ConfigZZpXGoldilocks256::DIM];
    s2_1_coeffs[..coeffs_s2_1.len()].copy_from_slice(&coeffs_s2_1);
    let s2_1 = RingGoldilock256::from_coefficients_vec_unchecked(s2_1_coeffs);

    let s2 = vec![s2_0, s2_1];

    // Compute the dot product <s1, s2>
    let result = RingGoldilock256::dot_product(&s1, &s2);

    // Manually compute the expected result
    // (1 + 2x)(5 + 6x) = 5 + 6x + 10x + 12x^2 = 5 + 16x + 12x^2
    // (3 + 4x)(7 + 8x) = 21 + 24x + 28x + 32x^2 = 21 + 52x + 32x^2
    // <s1, s2> = (5 + 16x + 12x^2) + (21 + 52x + 32x^2)
    //          = 26 + 68x + 44x^2

    let coeffs_expected = [
        Goldilocks::from(26u64), // Coefficient for x^0
        Goldilocks::from(68u64), // Coefficient for x^1
        Goldilocks::from(44u64), // Coefficient for x^2
    ];
    let mut expected_coeffs = vec![Goldilocks::from(0u64); ConfigZZpXGoldilocks256::DIM];
    expected_coeffs[..coeffs_expected.len()].copy_from_slice(&coeffs_expected);
    let expected = RingGoldilock256::from_coefficients_vec_unchecked(expected_coeffs);

    // Assert that the computed result matches the expected result
    assert_eq!(result, expected);
}
