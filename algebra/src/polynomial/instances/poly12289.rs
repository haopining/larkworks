use crate::{ConfigZZp12289, ConfigZZpX, ZZpX};

/// Configuration for ZZ[x]/(x^512+1) mod 12289
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub struct ConfigZZpX12289_512;

impl ConfigZZpX for ConfigZZpX12289_512 {
    /// Config for the base field
    type BaseConfig = ConfigZZp12289;
    /// Number of coefficients in a poly
    const DIM: usize = 512;
}

/// Polynomial with coefficient from ZZ_q where q=12289.
pub type Poly12289_512 = ZZpX<ConfigZZpX12289_512>;

#[test]
fn test_poly() {
    use crate::F12289;
    let coeffs = (0..ConfigZZpX12289_512::DIM)
        .map(|x| F12289::from(x as u64))
        .collect::<Vec<_>>()
        .try_into()
        .unwrap();
    let poly = Poly12289_512 { coeffs };
    println!("poly {}", poly);
    println!("poly {}", poly.clone() + poly);
}
