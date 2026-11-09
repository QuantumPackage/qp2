//! The EZFIO write: `ocaml/qp_create_ezfio.ml` `write_file` (lines
//! 282–695), transcribed.
//!
//! The sequence and the array `dim`s follow the reference line for line
//! (see `write-sequence.md`, parent-digested 2026-11-05, from a full read
//! of the source). Values are extracted from the parsed data structures
//! and passed through the pure computations in [`crate::specs`].
//!
//! # Errors
//!
//! One [`ToolError`] per failure (wave S2; the reference's caught vs
//! uncaught split had the same exit code and is collapsed):
//!
//! * the `Sys.file_exists` failure (`"<path> already exists"`) is
//!   reported by [`check_output_available`];
//! * a `mkdir` failure of the directory creation is the plain
//!   `cannot create <path>: <io error>`;
//! * a malformed pseudopotential block keeps the crate's plain text;
//! * any other failure of the write sequence removes the partial
//!   directory (best effort) and re-raises.

use std::fs;
use std::io;
use std::io::Write;
use std::path::Path;

use qp_ezfio::{Array, ByteLayout, EzFio};
use qp_types::basis::{self, Basis, LongBasis};
use qp_types::element;
use qp_types::error::ToolError;
use qp_types::pseudo::Pseudo;

use crate::Options;
use crate::TOOL;
use crate::basis_read::{AtomBasis, PseudoTable};
use crate::geometry::Geometry;
use crate::specs;

/// `if Sys.file_exists ezfio_file then failwith (ezfio_file ^ " already
/// exists")` (`qp_create_ezfio.ml:273`).
///
/// Positioned between the pseudo-table build and the per-nucleus basis
/// reads, so an existing output directory wins over a basis-read failure.
///
/// # Errors
///
/// The plain `<path> already exists` message when the path exists (a file
/// or a directory, like `Sys.file_exists`).
pub fn check_output_available(ezfio_file: &Path) -> Result<(), ToolError> {
    if ezfio_file.symlink_metadata().is_ok() {
        return Err(ToolError::new(
            TOOL,
            specs::exists_message(&ezfio_file.to_string_lossy()),
        ));
    }
    Ok(())
}

/// The whole action after the inputs: create the directory, write the
/// sections, run the trailing `qp_edit -c` (failure ignored) and print the
/// directory.
///
/// # Errors
///
/// See the module docs. Any failure after the directory exists removes it
/// (best effort) before re-raising, like the OCaml
/// `try write_file () with ex -> rmdir; raise ex`.
pub fn write(
    options: &Options,
    geometry: &Geometry,
    bases: &[AtomBasis],
    pseudo_table: &PseudoTable,
) -> Result<(), ToolError> {
    // The cross-module invariant the writer relies on: one `AtomBasis` per
    // nucleus, and `AtomBasis::nucleus` is the 1-based position in
    // `geometry.nuclei` (ghosts included).
    debug_assert!(
        bases.len() == geometry.nuclei.len(),
        "invariant violated: one AtomBasis per nucleus"
    );
    for (offset, atom_basis) in bases.iter().enumerate() {
        debug_assert!(
            usize::try_from(atom_basis.nucleus).ok() == Some(offset + 1),
            "invariant violated: AtomBasis::nucleus is not the 1-based position"
        );
    }
    let ezfio_file = &options.output;
    // The reference's `Ezfio.set_file` (ezfio.ml:426) runs `mkdir
    // ezfio_file` (ezfio.ml:69), a *single-level* `Unix.mkdir`: a missing
    // parent directory is the plain `cannot create` failure (wave S2).
    if let Some(parent) = ezfio_file.parent()
        && !parent.as_os_str().is_empty()
        && !parent.exists()
    {
        // Raw errno 2 (`ENOENT`), the same value on Linux and macOS.
        let error = io::Error::from_raw_os_error(2);
        return Err(ToolError::new(
            TOOL,
            format!("cannot create {}: {error}", ezfio_file.display()),
        ));
    }
    let ezfio = EzFio::create(ezfio_file).map_err(|error| {
        // The pre-check above covers the measured ENOENT case; the remaining
        // create failures (e.g. a parent that is a file) keep their own
        // display (the reference rendered every one of them as ENOENT — a
        // documented deviation).
        ToolError::new(
            TOOL,
            format!("cannot create {}: {error}", ezfio_file.display()),
        )
    })?;
    let ezfio = ezfio.with_layout(ByteLayout::Ocaml);
    // `set_file` also runs `mkdir (Filename.concat filename "ezfio")` —
    // the `ezfio` group exists with its `.version` marker even when it
    // holds no key.
    if let Err(error) = ezfio.create_group("ezfio") {
        let _ = fs::remove_dir_all(ezfio_file);
        return Err(app(error));
    }
    // `set_file` then stamps the group with machine state
    // (ezfio.ml:431-434). `creation` is not cosmetic: the Fortran
    // consumers of a created tree (e.g. `scf`, `ezfio.irp.f:30`) probe
    // for it and `stop 1` without it.
    if let Err(error) = write_stamps(ezfio_file) {
        let _ = fs::remove_dir_all(ezfio_file);
        return Err(error);
    }
    match write_file(&ezfio, options, geometry, bases, pseudo_table) {
        Ok(()) => {
            // `ignore @@ Sys.command ("qp_edit -c " ^ ezfio_file)`: the
            // program comes from `QP_EDIT` (test seam), its exit status is
            // ignored, and its stderr passes through.
            run_qp_edit(ezfio_file);
            let mut out = std::io::stdout().lock();
            let _ = out.write_all(format!("{}\n", ezfio_file.display()).as_bytes());
            let _ = out.flush();
            Ok(())
        }
        Err(error) => {
            let _ = fs::remove_dir_all(ezfio_file);
            Err(error)
        }
    }
}

/// `write_file ()`: the sections, in the reference's order.
fn write_file(
    ezfio: &EzFio,
    options: &Options,
    geometry: &Geometry,
    bases: &[AtomBasis],
    pseudo_table: &PseudoTable,
) -> Result<(), ToolError> {
    let nuclei = &geometry.nuclei;
    let nucl_num = nuclei.len();

    // -- Create + convention (282–285).
    ezfio
        .set::<i64>(
            "ezfio_files",
            "ezfio_convention",
            &specs::convention(options.qmc),
        )
        .map_err(app)?;
    ezfio
        .set::<bool>("basis", "ao_normalized", &true)
        .map_err(app)?;

    // -- Pseudo (287–485).
    let pseudo: Vec<Pseudo> = nuclei
        .iter()
        .map(|atom| -> Result<Pseudo, ToolError> {
            let symbol = atom.element;
            let pseudo = match pseudo_table.text(atom.element) {
                Some(text) => Pseudo::read_element(text, symbol)
                    .map_err(|error| ToolError::new(TOOL, error.to_string()))?,
                None => Pseudo::empty(symbol),
            };
            Ok(pseudo)
        })
        .collect::<Result<_, ToolError>>()?;

    let z_core: Vec<f64> = pseudo.iter().map(|p| removed(p) as f64).collect();
    let nucl_charge_remove = Array::new(vec![nucl_num], z_core).map_err(app)?;
    ezfio
        .set_array::<f64>("pseudo", "nucl_charge_remove", &nucl_charge_remove)
        .map_err(app)?;

    // Molecule adjustment (306–330): the core electrons are subtracted
    // from the electron counts (the odd one goes to beta) and from each
    // nucleus's charge.
    let n_remove: i64 = pseudo.iter().map(removed).sum();
    let (elec_alpha, elec_beta) = specs::adjust_electrons(
        i64::from(geometry.molecule.elec_alpha),
        i64::from(geometry.molecule.elec_beta),
        n_remove,
    );
    let charges: Vec<f64> = nuclei
        .iter()
        .zip(&pseudo)
        .map(|(atom, p)| atom.charge.to_float() - removed(p) as f64)
        .collect();

    ezfio
        .set::<i64>("electrons", "elec_alpha_num", &elec_alpha)
        .map_err(app)?;
    ezfio
        .set::<i64>("electrons", "elec_beta_num", &elec_beta)
        .map_err(app)?;

    // Nuclei (343–358): labels, adjusted charges, coordinates (bohr).
    ezfio
        .set::<i64>("nuclei", "nucl_num", &(nucl_num as i64))
        .map_err(app)?;
    let labels: Vec<String> = nuclei
        .iter()
        .map(|atom| element::to_string(atom.element).to_owned())
        .collect();
    let label_array = Array::new(vec![nucl_num], labels).map_err(app)?;
    ezfio
        .set_array::<String>("nuclei", "nucl_label", &label_array)
        .map_err(app)?;
    let charge_array = Array::new(vec![nucl_num], charges).map_err(app)?;
    ezfio
        .set_array::<f64>("nuclei", "nucl_charge", &charge_array)
        .map_err(app)?;
    // The reference never quantizes at the write site; its trailing
    // `qp_edit -c` (master `qp_create_ezfio.ml:694`, a shell-out whose exit
    // is ignored, ported as `run_qp_edit`) applied the parsed master-format
    // values, and since w12 that write-back quantizes the nuclei section.
    // The write path stays byte-identical to the reference's write path
    // (contract C3); the quantization is the tail's job (w12 review,
    // 2026-11-07: the spec's write-site grid application is retracted).
    let coords: Vec<(f64, f64, f64)> = nuclei
        .iter()
        .map(|atom| (atom.coord.x, atom.coord.y, atom.coord.z))
        .collect();
    let coord_array = Array::new(vec![nucl_num, 3], specs::flatten_coords(&coords)).map_err(app)?;
    ezfio
        .set_array::<f64>("nuclei", "nucl_coord", &coord_array)
        .map_err(app)?;

    // do_pseudo: the `-p` option itself, not the parsed potentials.
    ezfio
        .set::<bool>("pseudo", "do_pseudo", &options.pseudo.is_some())
        .map_err(app)?;

    // Pseudo maxima + the local / non-local matrices (368–485).
    let local_lens: Vec<usize> = pseudo.iter().map(|p| p.local.len()).collect();
    let projs: Vec<Vec<u64>> = pseudo
        .iter()
        .map(|p| p.non_local.iter().map(|(prim, _)| prim.proj()).collect())
        .collect();
    let projs_refs: Vec<&[u64]> = projs.iter().map(Vec::as_slice).collect();
    let (klocmax, kmax, lmax) = specs::pseudo_maxima(&local_lens, &projs_refs);

    ezfio
        .set::<i64>("pseudo", "pseudo_klocmax", &(klocmax as i64))
        .map_err(app)?;
    ezfio
        .set::<i64>("pseudo", "pseudo_kmax", &(kmax as i64))
        .map_err(app)?;
    ezfio
        .set::<i64>("pseudo", "pseudo_lmax", &(lmax as i64))
        .map_err(app)?;

    // Local: matrix [klocmax; nucl_num], row i = slot, column j = nucleus.
    let v_k = specs::local_matrix_data(nucl_num, klocmax, |i, j| {
        pseudo[j].local.get(i).map(|(_, coef)| *coef).unwrap_or(0.0)
    });
    let dz_k = specs::local_matrix_data(nucl_num, klocmax, |i, j| {
        pseudo[j]
            .local
            .get(i)
            .map(|(prim, _)| prim.expo())
            .unwrap_or(0.0)
    });
    let n_k: Vec<i64> = specs::local_matrix_data(nucl_num, klocmax, |i, j| {
        pseudo[j]
            .local
            .get(i)
            .map(|(prim, _)| prim.r_power())
            .unwrap_or(0) as f64
    })
    .into_iter()
    .map(|x| x as i64)
    .collect();
    for (key, data) in [("pseudo_v_k", v_k), ("pseudo_dz_k", dz_k)] {
        let array = Array::new(vec![nucl_num, klocmax], data).map_err(app)?;
        ezfio.set_array::<f64>("pseudo", key, &array).map_err(app)?;
    }
    {
        let array = Array::new(vec![nucl_num, klocmax], n_k).map_err(app)?;
        ezfio
            .set_array::<i64>("pseudo", "pseudo_n_k", &array)
            .map_err(app)?;
    }

    // Non-local: matrix [lmax+1; kmax; nucl_num], index (l = proj, i = slot
    // within the block, j = nucleus). The OCaml fill walks each nucleus's
    // `non_local` in file order, so the slot is the position of the entry
    // inside the subsequence with `proj = l`.
    let by_block: Vec<Vec<Vec<(f64, f64, i64)>>> = pseudo
        .iter()
        .map(|p| {
            (0..=lmax)
                .map(|l| {
                    p.non_local
                        .iter()
                        .filter(|(prim, _)| prim.proj() as usize == l)
                        .map(|(prim, coef)| (*coef, prim.expo(), prim.r_power()))
                        .collect()
                })
                .collect()
        })
        .collect();
    let get_nonlocal = |l: usize, i: usize, j: usize| -> Option<(f64, f64, i64)> {
        by_block[j][l].get(i).copied()
    };
    let v_kl = specs::nonlocal_matrix_data(nucl_num, kmax, lmax, |l, i, j| {
        get_nonlocal(l, i, j).map(|(c, _, _)| c).unwrap_or(0.0)
    });
    let dz_kl = specs::nonlocal_matrix_data(nucl_num, kmax, lmax, |l, i, j| {
        get_nonlocal(l, i, j).map(|(_, e, _)| e).unwrap_or(0.0)
    });
    let n_kl: Vec<i64> = specs::nonlocal_matrix_data(nucl_num, kmax, lmax, |l, i, j| {
        get_nonlocal(l, i, j).map(|(_, _, r)| r).unwrap_or(0) as f64
    })
    .into_iter()
    .map(|x| x as i64)
    .collect();
    for (key, data) in [("pseudo_v_kl", v_kl), ("pseudo_dz_kl", dz_kl)] {
        let array = Array::new(vec![nucl_num, kmax, lmax + 1], data).map_err(app)?;
        ezfio.set_array::<f64>("pseudo", key, &array).map_err(app)?;
    }
    {
        let array = Array::new(vec![nucl_num, kmax, lmax + 1], n_kl).map_err(app)?;
        ezfio
            .set_array::<i64>("pseudo", "pseudo_n_kl", &array)
            .map_err(app)?;
    }

    // -- Basis set (491–678).
    // `basis` = the per-nucleus reads concatenated in nucleus order
    // (the reference's `get_basis` result); every shell is already tagged
    // with its 1-based nucleus.
    let short: Basis = bases
        .iter()
        .flat_map(|atom_basis| atom_basis.basis.iter().cloned())
        .collect();
    let long = LongBasis::of_basis(&short);
    let ao_num = long.len();
    if ao_num > 0 {
        write_basis(ezfio, options, nucl_num, &short, &long)?;
    }
    Ok(())
}

/// The basis section (596–678), inside the reference's `if ao_num > 0`.
fn write_basis(
    ezfio: &EzFio,
    options: &Options,
    nucl_num: usize,
    short: &Basis,
    long: &LongBasis,
) -> Result<(), ToolError> {
    let ao_num = long.len();
    ezfio
        .set::<i64>("ao_basis", "ao_num", &(ao_num as i64))
        .map_err(app)?;
    ezfio
        .set::<String>("ao_basis", "ao_basis", &options.basis)
        .map_err(app)?;
    ezfio
        .set::<String>("basis", "basis", &options.basis)
        .map_err(app)?;

    let ao_prim_num: Vec<i64> = long
        .0
        .iter()
        .map(|(_, gto, _)| gto.lc.len() as i64)
        .collect();
    let ao_nucl: Vec<i64> = long.0.iter().map(|(_, _, n)| *n as i64).collect();
    let ao_power: Vec<i64> = specs::flatten_coords(
        &long
            .0
            .iter()
            .map(|(xyz, _, _)| (xyz.x as f64, xyz.y as f64, xyz.z as f64))
            .collect::<Vec<_>>(),
    )
    .into_iter()
    .map(|x| x as i64)
    .collect();

    // Zero-padded coefficient / exponent columns (the OCaml
    // `create_expo_coef`): data position `n * ao_num + j` is AO `j`'s
    // primitive `n` — i.e. `dim [ao_num; ao_prim_num_max]`, first index
    // fastest.
    let per_ao: Vec<Vec<(f64, f64)>> = long
        .0
        .iter()
        .map(|(_, gto, _)| {
            gto.lc
                .iter()
                .map(|(prim, coef)| (*coef, prim.expo))
                .collect()
        })
        .collect();
    let refs: Vec<&[(f64, f64)]> = per_ao.iter().map(Vec::as_slice).collect();
    let (padded, ao_prim_num_max) = specs::pad_prim_lists(&refs);
    let mut ao_coef = Vec::with_capacity(ao_num * ao_prim_num_max);
    let mut ao_expo = Vec::with_capacity(ao_num * ao_prim_num_max);
    // The AO index is the fastest, the primitive slot the slowest (the
    // OCaml `build` appends one `get_n n` block per slot); clippy's
    // iterator suggestion would transpose the stored order.
    #[allow(clippy::needless_range_loop)]
    for n in 0..ao_prim_num_max {
        for j in 0..ao_num {
            ao_coef.push(padded[j][n].0);
            ao_expo.push(padded[j][n].1);
        }
    }

    // Shell level (596–661): over the short basis, in shell order.
    let shell_num = short.len();
    let lc: Vec<&Vec<(qp_types::GaussianPrimitive, f64)>> =
        short.iter().map(|(gto, _)| &gto.lc).collect();
    // `Angmom.to_l` (`S` → 0 … `L` → 9).
    let ang_mom: Vec<i64> = lc
        .iter()
        .map(|l| {
            let sym = l.first().expect("a shell has at least one primitive").0.sym;
            i64::try_from(sym.to_l()).expect("to_l is 0..=9")
        })
        .collect();
    let mut expo: Vec<f64> = Vec::new();
    let mut coef: Vec<f64> = Vec::new();
    for l in &lc {
        for (prim, c) in l.iter() {
            expo.push(prim.expo);
            coef.push(*c);
        }
    }
    let shell_prim_num: Vec<i64> = lc.iter().map(|l| l.len() as i64).collect();
    let shell_idx = specs::shell_indices(
        &shell_prim_num
            .iter()
            .map(|&n| n as usize)
            .collect::<Vec<_>>(),
    );
    let prim_num = coef.len();
    let basis_nucleus_index: Vec<i64> = short.iter().map(|(_, n)| *n as i64).collect();
    let nucleus_shell_num = specs::shells_per_nucleus(nucl_num, &basis_nucleus_index);

    ezfio
        .set::<String>("basis", "typ", &"Gaussian".to_owned())
        .map_err(app)?;
    ezfio
        .set::<i64>("basis", "shell_num", &(shell_num as i64))
        .map_err(app)?;
    ezfio
        .set::<i64>("basis", "prim_num", &(prim_num as i64))
        .map_err(app)?;
    set_1d(ezfio, "basis", "shell_prim_num", &shell_prim_num)?;
    set_1d(ezfio, "basis", "shell_ang_mom", &ang_mom)?;
    set_1d(ezfio, "basis", "shell_index", &shell_idx)?;
    set_1d(ezfio, "basis", "basis_nucleus_index", &basis_nucleus_index)?;
    set_1d(ezfio, "basis", "nucleus_shell_num", &nucleus_shell_num)?;
    set_1d_f64(ezfio, "basis", "prim_coef", &coef)?;
    set_1d_f64(ezfio, "basis", "prim_expo", &expo)?;

    // AO-level arrays (664–674).
    set_1d(ezfio, "ao_basis", "ao_prim_num", &ao_prim_num)?;
    set_1d(ezfio, "ao_basis", "ao_nucl", &ao_nucl)?;
    {
        let array = Array::new(vec![ao_num, 3], ao_power).map_err(app)?;
        ezfio
            .set_array::<i64>("ao_basis", "ao_power", &array)
            .map_err(app)?;
    }
    {
        let array = Array::new(vec![ao_num, ao_prim_num_max], ao_coef).map_err(app)?;
        ezfio
            .set_array::<f64>("ao_basis", "ao_coef", &array)
            .map_err(app)?;
    }
    {
        let array = Array::new(vec![ao_num, ao_prim_num_max], ao_expo).map_err(app)?;
        ezfio
            .set_array::<f64>("ao_basis", "ao_expo", &array)
            .map_err(app)?;
    }
    ezfio
        .set::<bool>("ao_basis", "ao_cartesian", &options.cartesian)
        .map_err(app)?;

    // The md5 stamp (676–678): the reference reads the `ao_basis` group
    // back and stamps `Basis.to_md5` over the short basis it rebuilds —
    // here computed on the in-memory short basis (the `%24.15e` EZFIO
    // round trip is a no-op for file-derived values; the 27-row corpus
    // verified this in P3). `ao_extra_md5` is NOT stamped here: that is
    // `qp_edit`'s job.
    let md5 = basis::to_md5(short);
    ezfio
        .set::<String>("ao_basis", "ao_md5", &md5)
        .map_err(app)?;
    Ok(())
}

/// `ignore @@ Sys.command ("<QP_EDIT or qp_edit> -c " ^ dir)`: through
/// `/bin/sh -c`, so a missing program prints the shell's own
/// `command not found` and the exit status is discarded.
fn run_qp_edit(ezfio_file: &Path) {
    let program = std::env::var("QP_EDIT").unwrap_or_else(|_| "qp_edit".to_owned());
    let _ = std::process::Command::new("/bin/sh")
        .arg("-c")
        .arg(format!("{program} -c {}", ezfio_file.display()))
        .status();
}

/// `qp_ezfio` errors are port bugs, not reference failure classes.
fn app(error: qp_ezfio::Error) -> ToolError {
    ToolError::new(TOOL, format!("internal EZFIO error: {error}"))
}

/// The three `set_file` stamps (ezfio.ml:431-434), written after the
/// `ezfio` group exists:
///
/// ```sh
/// LANG= date > <file>/ezfio/creation
/// echo $USER > <file>/ezfio/user
/// echo <library> > <file>/ezfio/library
/// ```
///
/// The reference folds them into one shell command and fails with a
/// `Failure` naming the file when it fails; the port runs `date`
/// directly (the reference's sprintf'd command would misparse a path
/// with spaces) and reports the plain one-line failure. `creation` must
/// exist: the Fortran consumers of a created tree (e.g. `scf`,
/// `ezfio.irp.f:30`) probe for it and `stop 1` without it. All three
/// are VOLATILE in the golden harness, so their content is machine
/// state.
fn write_stamps(ezfio_file: &Path) -> Result<(), ToolError> {
    let failure = || {
        ToolError::new(
            TOOL,
            format!("cannot write the ezfio stamps in {}", ezfio_file.display()),
        )
    };
    let dir = ezfio_file.join("ezfio");
    // `LANG= date`: the local time in the C locale.
    let date = std::process::Command::new("date")
        .env_remove("LANG")
        .output()
        .map_err(|_| failure())?;
    if !date.status.success() {
        return Err(failure());
    }
    // `echo $USER`: empty line when the variable is unset.
    let user = std::env::var("USER").unwrap_or_default();
    // The reference stamps its configure-time ezfio library path
    // (`ezfio.ml:28`). A static pure-Rust build has no such path, and
    // nothing in the reference tool layer reads the file back; the
    // built-in implementation's name stands in (documented deviation).
    const LIBRARY: &str = "qp-ezfio";
    if fs::write(dir.join("creation"), &date.stdout).is_err()
        || fs::write(dir.join("user"), format!("{user}\n").as_bytes()).is_err()
        || fs::write(dir.join("library"), format!("{LIBRARY}\n").as_bytes()).is_err()
    {
        return Err(failure());
    }
    Ok(())
}

/// Rank-1 `i64` array write (dim = data length).
fn set_1d(ezfio: &EzFio, group: &str, key: &str, data: &[i64]) -> Result<(), ToolError> {
    let array = Array::new(vec![data.len()], data.to_vec()).map_err(app)?;
    ezfio.set_array::<i64>(group, key, &array).map_err(app)
}

/// Rank-1 `f64` array write (dim = data length).
fn set_1d_f64(ezfio: &EzFio, group: &str, key: &str, data: &[f64]) -> Result<(), ToolError> {
    let array = Array::new(vec![data.len()], data.to_vec()).map_err(app)?;
    ezfio.set_array::<f64>(group, key, &array).map_err(app)
}

/// `Positive_int.to_int x.Pseudo.n_elec`: the removed electron count as a
/// signed integer (a count, realistically ≤ 118).
fn removed(p: &Pseudo) -> i64 {
    i64::try_from(p.n_elec).unwrap_or(i64::MAX)
}
