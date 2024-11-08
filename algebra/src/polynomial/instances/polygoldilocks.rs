use crate::{ConfigZZpGoldilocks, ConfigZZpX, ZZpX};

/// Configuration for ZZ[x]/(x^512+1) mod 12289
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub struct ConfigZZpXGoldilocks256;

impl ConfigZZpX for ConfigZZpXGoldilocks256 {
    /// Config for the base field
    type BaseConfig = ConfigZZpGoldilocks;
    /// Number of coefficients in a poly
    const DIM: usize = 256;
}

/// Polynomial with coefficient from ZZ_q where q=12289.
pub type PolyGoldilock256 = ZZpX<ConfigZZpXGoldilocks256>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_poly() {
        use crate::Goldilocks;
        // ConfigZZpXGoldilocks256::DIM == numbers of coefficients
        let coeffs = (0..ConfigZZpXGoldilocks256::DIM)
            .map(|x| Goldilocks::from(x as u64))
            .collect::<Vec<_>>()
            .try_into()
            .unwrap();
        let poly = PolyGoldilock256 { coeffs };
        println!("poly {}", poly);
        println!("poly {}", poly.clone() + poly);
    }

    #[test]
    fn test_infinity_norm() {
        use crate::Goldilocks;
        use crate::Polynomial;
        let coeffs = (0..ConfigZZpXGoldilocks256::DIM)
            .map(|x| Goldilocks::from(x as u64))
            .collect::<Vec<_>>()
            .try_into()
            .unwrap();
        let poly = PolyGoldilock256 { coeffs };
        let expected = (0..ConfigZZpXGoldilocks256::DIM)
            .map(|x| x as u64)
            .max()
            .unwrap();
        let result = poly.infinity_norm();
        assert_eq!(expected, result.try_into().unwrap());
    }

    #[test]
    fn test_l2_norm() {
        use crate::Goldilocks;
        use crate::Polynomial;
        use num::integer::Roots;

        let coeffs = (0..ConfigZZpXGoldilocks256::DIM)
            .map(|x| Goldilocks::from(x as u64))
            .collect::<Vec<_>>()
            .try_into()
            .unwrap();
        let poly = PolyGoldilock256 { coeffs };
        let expected = (0..ConfigZZpXGoldilocks256::DIM)
            .map(|x| (x * x) as u64)
            .sum::<u64>()
            .sqrt();

        println!("expected {:?}", expected);
        let result = poly.l2_norm();
        assert_eq!(expected, result.try_into().unwrap());
    }
}
