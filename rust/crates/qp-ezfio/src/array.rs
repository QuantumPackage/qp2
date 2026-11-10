//! Rank-`n` arrays with Fortran (column-major) element order.

use crate::error::{Error, Result};

/// An EZFIO array: dimensions plus the flattened elements.
///
/// EZFIO stores arrays as a rank followed by the extents followed by the
/// elements in **Fortran order** — the *first* dimension varies fastest, exactly
/// like `reshape(dat, dims)` in Fortran and like the reference OCaml client's
/// `flattened_ezfio`. `data.len()` is therefore always
/// `dims.iter().product::<usize>()`.
///
/// Indices in the accessors are **1-based** to match the Fortran/IRPF90 side
/// (e.g. `psi_det` is indexed `(1:n_int*bit_kind/8, 1:2, 1:n_det)`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Array<T> {
    dims: Vec<usize>,
    data: Vec<T>,
}

impl<T> Array<T> {
    /// Build an array from explicit dimensions and flattened Fortran-order data.
    ///
    /// Fails with [`Error::InconsistentArray`] if `data.len()` is not exactly
    /// the product of `dims` (the invariant every reader and writer relies on).
    pub fn new(dims: Vec<usize>, data: Vec<T>) -> Result<Self> {
        let expected = dims
            .iter()
            .try_fold(1usize, |acc, &d| acc.checked_mul(d))
            .unwrap_or(usize::MAX);
        if expected != data.len() {
            return Err(Error::InconsistentArray {
                dims,
                expected,
                found: data.len(),
            });
        }
        Ok(Self { dims, data })
    }

    /// Build a rank-1 array from a vector.
    pub fn from_vec(data: Vec<T>) -> Self {
        let dims = vec![data.len()];
        Self { dims, data }
    }

    /// Dimensions, outermost index last (Fortran order).
    pub fn dims(&self) -> &[usize] {
        &self.dims
    }

    /// Rank: the number of dimensions. Always at least 1 for a value that came
    /// from an EZFIO file.
    pub fn rank(&self) -> usize {
        self.dims.len()
    }

    /// Total number of elements.
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// True when the array holds no element (some EZFIO arrays have a zero
    /// extent, e.g. `pseudo/pseudo_dz_k` before pseudopotentials are set).
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// The flattened elements, first dimension varying fastest.
    pub fn data(&self) -> &[T] {
        &self.data
    }

    /// Mutable access to the flattened elements.
    pub fn data_mut(&mut self) -> &mut [T] {
        &mut self.data
    }

    /// Consume the array and return the flattened elements.
    pub fn into_data(self) -> Vec<T> {
        self.data
    }

    /// Elements in Fortran order.
    pub fn iter(&self) -> std::slice::Iter<'_, T> {
        self.data.iter()
    }

    /// Linear (0-based) index of a 1-based multi-index.
    ///
    /// `indices.len()` must equal the rank and every index must be within
    /// `1..=dims[i]`.
    ///
    /// Fails with [`Error::IndexOutOfBounds`] on a rank or a value mismatch.
    /// The linear index cannot overflow: [`Array::new`] has already checked
    /// that the product of the dimensions fits in a `usize`, so every valid
    /// multi-index maps below `data.len()`.
    pub fn index_of(&self, indices: &[usize]) -> Result<usize> {
        let out_of_bounds = || Error::IndexOutOfBounds {
            dims: self.dims.clone(),
            indices: indices.to_vec(),
        };
        if indices.len() != self.dims.len() {
            return Err(out_of_bounds());
        }
        let mut linear = 0usize;
        let mut stride = 1usize;
        for (&index, &dim) in indices.iter().zip(self.dims.iter()) {
            if index == 0 || index > dim {
                return Err(out_of_bounds());
            }
            linear += (index - 1) * stride;
            stride *= dim;
        }
        Ok(linear)
    }

    /// Element at a 1-based multi-index, or `None` if out of bounds.
    pub fn at(&self, indices: &[usize]) -> Option<&T> {
        self.index_of(indices).ok().map(|i| &self.data[i])
    }
}

impl<T> From<Vec<T>> for Array<T> {
    fn from(data: Vec<T>) -> Self {
        Self::from_vec(data)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fortran_order_indexing() {
        // dims (2,3): index 1 varies fastest.
        let a = Array::new(vec![2, 3], (1..=6).collect::<Vec<i32>>()).unwrap();
        assert_eq!(a.rank(), 2);
        assert_eq!(a.len(), 6);
        assert_eq!(a.at(&[1, 1]), Some(&1));
        assert_eq!(a.at(&[2, 1]), Some(&2));
        assert_eq!(a.at(&[1, 2]), Some(&3));
        assert_eq!(a.at(&[2, 3]), Some(&6));
        assert_eq!(a.at(&[3, 1]), None);
        assert_eq!(a.at(&[1, 4]), None);
        assert_eq!(a.at(&[1]), None);
    }

    #[test]
    fn inconsistent_dims_are_rejected() {
        let err = Array::new(vec![2, 3], vec![1, 2, 3]).unwrap_err();
        assert!(matches!(err, Error::InconsistentArray { .. }));
    }

    #[test]
    fn zero_extent_is_allowed() {
        let a = Array::new(vec![3, 0], Vec::<f64>::new()).unwrap();
        assert!(a.is_empty());
        assert_eq!(a.dims(), &[3, 0]);
    }
}
