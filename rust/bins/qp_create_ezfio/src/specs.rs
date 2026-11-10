//! Pure value computations of the EZFIO write.
//!
//! Every function here is free of `qp_types` / `qp_ezfio` types so each one
//! is unit-testable in isolation; [`crate::write`] extracts the raw numbers
//! from the parsed data structures and calls these. The computations
//! reproduce `ocaml/qp_create_ezfio.ml` `write_file` line for line:
//! the OCaml folds start at the same bases (`0`), the divisions are OCaml
//! integer divisions (floor for the positive values involved), and the
//! array flattenings are the OCaml `concat_2d` / `concat_3d`
//! (row-major over the built matrices, which is column-major over the
//! stored `dim` — first index fastest).

/// `Ezfio.set_ezfio_files_ezfio_convention`: `if qmc then 20210101 else
/// 20250211`.
#[must_use]
pub const fn convention(qmc: bool) -> i64 {
    if qmc { 20210101 } else { 20250211 }
}

/// The `"<path> already exists"` payload of the `Sys.file_exists` check
/// (`qp_create_ezfio.ml:273`).
#[must_use]
pub fn exists_message(path: &str) -> String {
    format!("{path} already exists")
}

/// `n_elec_to_remove` folded over the pseudopotentials, applied to the
/// electron counts:
///
/// ```text
/// elec_alpha' = alpha − n / 2
/// elec_beta'  = beta  − (n − n / 2)
/// ```
///
/// OCaml integer division on the (non-negative) values involved, so the
/// odd electron is removed from beta.
#[must_use]
pub fn adjust_electrons(alpha: i64, beta: i64, n_remove: i64) -> (i64, i64) {
    (alpha - n_remove / 2, beta - (n_remove - n_remove / 2))
}

/// `klocmax`, `kmax` and `lmax` of the pseudo section
/// (`qp_create_ezfio.ml:368–404`).
///
/// * `local_lens[j]` — the number of local primitives of nucleus `j`;
/// * `projs[j]` — the `proj` values of nucleus `j`'s non-local primitives,
///   in file order.
///
/// The OCaml folds start at `0`, so with no local components `klocmax` is
/// `0` and with no non-local components `lmax` is `0` (one empty block).
/// `kmax` is the maximum, over every block `l ≤ lmax`, of the number of
/// non-local primitives with `proj = l` on the busiest nucleus.
#[must_use]
pub fn pseudo_maxima(local_lens: &[usize], projs: &[&[u64]]) -> (usize, usize, usize) {
    let klocmax = local_lens.iter().copied().max().unwrap_or(0);
    let lmax = projs
        .iter()
        .flat_map(|list| list.iter().copied())
        .max()
        .unwrap_or(0) as usize;
    let mut kmax = 0_usize;
    for l in 0..=lmax {
        let l = l as u64;
        for list in projs {
            let count = list.iter().filter(|p| **p == l).count();
            if count > kmax {
                kmax = count;
            }
        }
    }
    (klocmax, kmax, lmax)
}

/// Zero-pad the per-AO `(coefficient, exponent)` lists to
/// `ao_prim_num_max` columns (the OCaml `get_n` returns `0.` past the end).
///
/// Returns the padded lists and the pad width.
#[must_use]
pub fn pad_prim_lists(pairs: &[&[(f64, f64)]]) -> (Vec<Vec<(f64, f64)>>, usize) {
    let max = pairs.iter().map(|p| p.len()).max().unwrap_or(0);
    let padded = pairs
        .iter()
        .map(|pair| {
            let mut out: Vec<(f64, f64)> = pair.to_vec();
            out.resize_with(max, || (0.0, 0.0));
            out
        })
        .collect();
    (padded, max)
}

/// `nucl_coord` data: all `x`, then all `y`, then all `z` (the OCaml
/// `xs @ ys @ zs`), stored under `dim [nucl_num; 3]`, first index fastest.
#[must_use]
pub fn flatten_coords(coords: &[(f64, f64, f64)]) -> Vec<f64> {
    let mut data = Vec::with_capacity(3 * coords.len());
    for (x, _, _) in coords {
        data.push(*x);
    }
    for (_, y, _) in coords {
        data.push(*y);
    }
    for (_, _, z) in coords {
        data.push(*z);
    }
    data
}

/// `shell_index`: for every primitive, the 1-based index of the shell it
/// belongs to (length `prim_num`).
#[must_use]
pub fn shell_indices(shell_prims: &[usize]) -> Vec<i64> {
    let total: usize = shell_prims.iter().sum();
    let mut data = Vec::with_capacity(total);
    for (shell, prim) in shell_prims.iter().enumerate() {
        for _ in 0..*prim {
            data.push((shell + 1) as i64);
        }
    }
    data
}

/// `nucleus_shell_num`: per nucleus (1-based in `shell_nuclei`), the number
/// of shells on that nucleus (length `nucl_num`).
#[must_use]
pub fn shells_per_nucleus(nucl_num: usize, shell_nuclei: &[i64]) -> Vec<i64> {
    let mut counts = vec![0_i64; nucl_num];
    for &nucleus in shell_nuclei {
        let index = (nucleus - 1) as usize;
        if index < counts.len() {
            counts[index] += 1;
        }
    }
    counts
}

/// Flatten the local pseudo matrices the way OCaml `concat_2d` does: the
/// built matrices are `[klocmax; nucl_num]` (row `i` = primitive slot,
/// column `j` = nucleus), flattened row-major, and stored under
/// `dim [nucl_num; klocmax]` — first index fastest, i.e. data position
/// `j + nucl_num * i`.
#[must_use]
pub fn local_matrix_data(
    nucl_num: usize,
    klocmax: usize,
    get: impl Fn(usize, usize) -> f64,
) -> Vec<f64> {
    let mut data = Vec::with_capacity(nucl_num * klocmax);
    for i in 0..klocmax {
        for j in 0..nucl_num {
            data.push(get(i, j));
        }
    }
    data
}

/// Flatten the non-local pseudo matrices the way OCaml `concat_3d` does:
/// the built matrices are `[lmax+1; kmax; nucl_num]` (block `l = proj`,
/// slot `i` within the block, nucleus `j`), flattened row-major over
/// `(l, i, j)`, and stored under `dim [nucl_num; kmax; lmax+1]` — data
/// position `j + nucl_num * i + nucl_num * kmax * l`.
#[must_use]
pub fn nonlocal_matrix_data(
    nucl_num: usize,
    kmax: usize,
    lmax: usize,
    get: impl Fn(usize, usize, usize) -> f64,
) -> Vec<f64> {
    let mut data = Vec::with_capacity(nucl_num * kmax * (lmax + 1));
    for l in 0..=lmax {
        for i in 0..kmax {
            for j in 0..nucl_num {
                data.push(get(l, i, j));
            }
        }
    }
    data
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_convention_follows_the_qmc_flag() {
        assert_eq!(convention(false), 20250211);
        assert_eq!(convention(true), 20210101);
    }

    #[test]
    fn the_exists_message_keeps_the_path_verbatim() {
        assert_eq!(
            exists_message("/tmp/x.ezfio"),
            "/tmp/x.ezfio already exists"
        );
        assert_eq!(exists_message("rel.ezfio"), "rel.ezfio already exists");
    }

    #[test]
    fn the_electron_adjustment_sends_the_odd_one_to_beta() {
        assert_eq!(adjust_electrons(2, 2, 0), (2, 2));
        assert_eq!(adjust_electrons(5, 5, 2), (4, 4));
        assert_eq!(adjust_electrons(5, 5, 3), (4, 3));
        assert_eq!(adjust_electrons(6, 4, 5), (4, 1));
    }

    #[test]
    fn the_pseudo_maxima_fold_from_zero() {
        // No pseudo at all: klocmax 0, lmax 0 (one empty block), kmax 0.
        assert_eq!(pseudo_maxima(&[], &[&[]]), (0, 0, 0));
        assert_eq!(pseudo_maxima(&[0, 0], &[&[], &[]]), (0, 0, 0));
        // Two local slots on the busiest nucleus.
        assert_eq!(pseudo_maxima(&[1, 2, 0], &[&[], &[], &[]]), (2, 0, 0));
        // projs [0; 0; 2] on one nucleus: lmax 2, kmax 2 (block 0 twice).
        assert_eq!(pseudo_maxima(&[0], &[&[0, 0, 2]]), (0, 2, 2));
        // The busiest nucleus wins per block, across nuclei.
        assert_eq!(pseudo_maxima(&[0], &[&[0, 0], &[1]]), (0, 2, 1));
    }

    #[test]
    fn the_padding_reaches_the_widest_ao() {
        let a = [(1.0, 2.0)];
        let b = [(3.0, 4.0), (5.0, 6.0)];
        let (padded, max) = pad_prim_lists(&[&a, &b, &b]);
        assert_eq!(max, 2);
        assert_eq!(padded[0], vec![(1.0, 2.0), (0.0, 0.0)]);
        assert_eq!(padded[1], vec![(3.0, 4.0), (5.0, 6.0)]);
    }

    #[test]
    fn the_coord_flattening_is_x_then_y_then_z() {
        assert_eq!(
            flatten_coords(&[(1.0, 11.0, 111.0), (2.0, 12.0, 112.0)]),
            vec![1.0, 2.0, 11.0, 12.0, 111.0, 112.0]
        );
        assert!(flatten_coords(&[]).is_empty());
    }

    #[test]
    fn the_shell_indices_are_one_based_per_primitive() {
        assert_eq!(shell_indices(&[1, 3, 0]), vec![1, 2, 2, 2]);
        assert!(shell_indices(&[]).is_empty());
    }

    #[test]
    fn the_shell_counts_cover_every_nucleus() {
        assert_eq!(shells_per_nucleus(3, &[1, 1, 2, 3, 1]), vec![3, 1, 1]);
        assert_eq!(shells_per_nucleus(2, &[]), vec![0, 0]);
    }

    #[test]
    fn the_local_flattening_matches_column_major_storage() {
        // Built [2; 3]: rows i = 0,1 ; columns j = 0,1,2.
        let data = local_matrix_data(3, 2, |i, j| (i * 10 + j) as f64);
        // (j fastest): 0,1,2, 10,11,12.
        assert_eq!(data, vec![0.0, 1.0, 2.0, 10.0, 11.0, 12.0]);
        assert!(local_matrix_data(3, 0, |_, _| 0.0).is_empty());
    }

    #[test]
    fn the_nonlocal_flattening_matches_column_major_storage() {
        // Built [2; 2; 2]: (l, i, j) -> l*4 + i*2 + j.
        let data = nonlocal_matrix_data(2, 2, 1, |l, i, j| (l * 4 + i * 2 + j) as f64);
        // stored (j fastest, then i, then l): 0,1, 2,3, 4,5, 6,7.
        assert_eq!(data, vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]);
    }
}
