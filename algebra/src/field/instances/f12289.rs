use crate::{field::zz_p::ZZp, ConfigZZp};

/// Configuration parameter for ZZ mod 12289
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub struct ConfigZZp12289;

impl ConfigZZp for ConfigZZp12289 {
    type PrimitiveType = u16;
    type ProductType = u32;
    const MODULUS: Self::PrimitiveType = 12289;
    /// The place where the multiplication algorithm is actually implemented.
    fn mul_internal(a: &Self::PrimitiveType, b: &Self::PrimitiveType) -> Self::PrimitiveType {
        (*a as Self::ProductType * *b as Self::ProductType % Self::MODULUS as Self::ProductType)
            as Self::PrimitiveType
    }

    /// The place where the addition algorithm is actually implemented.
    fn add_internal(a: &Self::PrimitiveType, b: &Self::PrimitiveType) -> Self::PrimitiveType {
        let mut tmp = a + b;
        if tmp >= Self::MODULUS {
            tmp -= Self::MODULUS
        }
        tmp
    }

    /// The place where the subtraction algorithm is actually implemented.
    fn sub_internal(a: &Self::PrimitiveType, b: &Self::PrimitiveType) -> Self::PrimitiveType {
        if a >= b {
            a - b
        } else {
            a + Self::MODULUS - b
        }
    }

    fn eq_internal(a: &Self::PrimitiveType, b: &Self::PrimitiveType) -> bool {
        a % Self::MODULUS == b % Self::MODULUS
    }
}

/// ZZ mod 12289
pub type F12289 = ZZp<ConfigZZp12289>;

#[cfg(test)]
mod tests {
    use super::F12289;
    use crate::tests::field::random_field_tests;

    #[test]
    fn test_integer() {
        random_field_tests::<F12289>("F12289".to_string());
    }
}
