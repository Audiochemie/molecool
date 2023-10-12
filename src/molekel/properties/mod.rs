use std::ops::Div;

use super::Molecule;
use itertools::Itertools;
use nalgebra::DMatrix;
use num::{traits::AsPrimitive, Float};

pub struct Hessian<T>(pub DMatrix<T>);

pub fn mass_weight_hessian<T>(m: &Molecule<T>, h: &mut Hessian<T>)
where
    T: Float
        + std::fmt::Debug
        + std::str::FromStr
        + nalgebra::ClosedMul
        + nalgebra::ClosedAdd
        + nalgebra::ClosedDiv
        + num::cast::AsPrimitive<T>
        + 'static,
    <T as std::str::FromStr>::Err: std::fmt::Debug,
    f64: num::traits::AsPrimitive<T>,
    f32: num::traits::AsPrimitive<T>,
{
    let n_atms = m.building_atoms.len();
    if h.0.shape() != (n_atms * 3, n_atms * 3) {
        panic!(
            "Shape of hessian {:?} does not match (NAtms={} * 3, NAtms * 3) ",
            h.0.shape(),
            n_atms
        )
    }
    let mut mass_scaling: Vec<f64> = Vec::with_capacity(n_atms * 3 * n_atms * 3);
    let cart_prod = m.building_atoms.iter().cartesian_product(0..3);
    let cart_prod = cart_prod.clone().cartesian_product(cart_prod);
    for ((at1, _c1), (at2, _c2)) in cart_prod {
        mass_scaling.push(1.0.div((at1.atomic_mass * at2.atomic_mass).sqrt()));
    }
    for (j, ele) in h.0.iter_mut().enumerate() {
        *ele *= mass_scaling[j].as_()
    }
}
