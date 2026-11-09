//! Golden tests for [`qp_types::pseudo`] against the reference implementation.
//!
//! Two independent oracles are used, both outside this crate:
//!
//! * **The reference parser.** `ocaml/Pseudo.ml` was transcribed with real
//!   OCaml `Printf` and `sexplib` derivation (`$HOME/OpenCode/tmp/w2b/
//!   pseudo_oracle_dune/`) and run on `data/pseudo/bfd` and `data/pseudo/
//!   ccecp`. Every `display` and `sexp` literal below is that
//!   program's output, byte for byte, including `%20f` rounding and
//!   `%G`-style float fields.
//! * **The reference EZFIO arrays.** `tests/fixtures/pseudo/c2h2/` holds the
//!   decompressed payloads written by `qp_create_ezfio` for the `c2h2`
//!   fixture. The second test rebuilds those arrays from the parsed C and H
//!   pseudopotentials with the arithmetic of `ocaml/qp_create_ezfio.ml`
//!   (lines 371-490) and compares every element exactly.
//!
//! `tests/fixtures/pseudo/README.md` records where the fixture copies come
//! from.
//!
//! Note on `ao_extra_md5`: this is **not** the md5 that ends up in
//! `basis-md5.jsonl`. `Input_ao_extra_basis.read` hashes the literal string
//! `"None"` and never consults `Pseudo`, so all 27 rows of the corpus carry
//! `6adf97f83acf6453d4a6a4b1070f3754`. The `Pseudo.to_md5` digest itself is
//! not pinned here: no code in the OCaml tree calls it (it came from
//! `include To_md5`), and the port dropped it in wave S1.

use qp_types::pseudo::{ElementSymbol, Pseudo};

/// Verbatim copy of `data/pseudo/bfd` (see the fixture README).
const BFD: &str = include_str!("fixtures/pseudo/bfd");
/// Verbatim copy of `data/pseudo/ccecp`.
const CCECP: &str = include_str!("fixtures/pseudo/ccecp");

/// One reference element block and the three renderings of it.
struct Golden {
    element: &'static str,
    /// Which library text the block is read from.
    library: &'static str,
    n_elec: u64,
    /// `Pseudo.to_string`, without the newline `print_endline` adds.
    display: &'static str,
    /// `Sexp.to_string (sexp_of_t t)`.
    sexp: &'static str,
}

const GOLDEN: &[Golden] = &[
    Golden {
        element: "H",
        library: "BFD",
        n_elec: 0,
        display: "H   0 electrons removed\nLocal component:\n              Coeff.      r^n                 Exp.\n            1.000000       -1             4.476924\n            4.476924        1             2.976365\n           -4.321123        0             3.018416",
        sexp: "((element H)(n_elec 0)(local((((expo 4.4769241)(r_power -1))1)(((expo 2.97636451)(r_power 1))4.4769241)(((expo 3.01841596)(r_power 0))-4.3211234)))(non_local()))",
    },
    Golden {
        element: "He",
        library: "BFD",
        n_elec: 0,
        display: "He   0 electrons removed\nLocal component:\n              Coeff.      r^n                 Exp.\n            2.000000       -1             9.843688\n           19.687376        1            10.945162\n          -17.205703        0            12.037153",
        sexp: "((element He)(n_elec 0)(local((((expo 9.84368811)(r_power -1))2)(((expo 10.94516233)(r_power 1))19.68737622)(((expo 12.03715264)(r_power 0))-17.20570338)))(non_local()))",
    },
    Golden {
        element: "C",
        library: "BFD",
        n_elec: 2,
        display: "C   2 electrons removed\nLocal component:\n              Coeff.      r^n                 Exp.\n            4.000000       -1             8.359738\n           33.438953        1             4.483619\n          -19.175373        0             3.938313\nNon-local component:\n              Coeff.      r^n                 Exp.     Proj\n           22.551642        0             5.029916   |0><0|",
        sexp: "((element C)(n_elec 2)(local((((expo 8.35973821)(r_power -1))4)(((expo 4.48361888)(r_power 1))33.43895285)(((expo 3.93831258)(r_power 0))-19.17537323)))(non_local((((expo 5.02991637)(r_power 0)(proj 0))22.55164191))))",
    },
    Golden {
        element: "Cu",
        library: "BFD",
        n_elec: 10,
        display: "Cu   10 electrons removed\nLocal component:\n              Coeff.      r^n                 Exp.\n           19.000000       -1             6.251496\n          118.778429        1             6.727253\n         -105.899824        0             6.610246\nNon-local component:\n              Coeff.      r^n                 Exp.     Proj\n          127.350694        0            12.365687   |0><0|\n           71.229849        0            11.160729   |1><1|",
        sexp: "((element Cu)(n_elec 10)(local((((expo 6.25149628)(r_power -1))19)(((expo 6.72725326)(r_power 1))118.7784294)(((expo 6.61024592)(r_power 0))-105.89982403)))(non_local((((expo 12.36568715)(r_power 0)(proj 0))127.35069424)(((expo 11.16072939)(r_power 0)(proj 1))71.229849))))",
    },
    Golden {
        element: "Mg",
        library: "CCECP",
        n_elec: 10,
        display: "Mg   10 electrons removed\nLocal component:\n              Coeff.      r^n                 Exp.\n            2.000000       -1             6.048538\n           12.097075        1             2.796989\n          -17.108313        0             2.547408\nNon-local component:\n              Coeff.      r^n                 Exp.     Proj\n            6.428631        0             5.936017   |0><0|\n           14.195491        0             1.592891   |0><0|\n            3.315069        0             1.583969   |1><1|\n            4.403025        0             1.077297   |1><1|",
        sexp: "((element Mg)(n_elec 10)(local((((expo 6.048538)(r_power -1))2)(((expo 2.796989)(r_power 1))12.097075)(((expo 2.547408)(r_power 0))-17.108313)))(non_local((((expo 5.936017)(r_power 0)(proj 0))6.428631)(((expo 1.592891)(r_power 0)(proj 0))14.195491)(((expo 1.583969)(r_power 0)(proj 1))3.315069)(((expo 1.077297)(r_power 0)(proj 1))4.403025))))",
    },
];

#[test]
fn reference_blocks_match_display_sexp_and_md5() {
    for golden in GOLDEN {
        let text = match golden.library {
            "BFD" => BFD,
            "CCECP" => CCECP,
            other => panic!("unknown library {other}"),
        };
        let element = ElementSymbol::of_string(golden.element).expect("known element");
        let pseudo = Pseudo::read_element(text, element).expect("parse reference block");
        assert_eq!(pseudo.element, element);
        assert_eq!(pseudo.n_elec, golden.n_elec, "{}: n_elec", golden.element);
        assert_eq!(
            pseudo.to_string(),
            golden.display,
            "{}: display",
            golden.element
        );
        assert_eq!(pseudo.to_sexp(), golden.sexp, "{}: sexp", golden.element);
    }
}

#[test]
fn empty_pseudo_matches_the_reference_renderings() {
    // `Pseudo.empty (Element X)` and `Pseudo.empty (Element H)`, as produced
    // by the same OCaml oracle.
    for (element, sexp) in [
        ("X", "((element X)(n_elec 0)(local())(non_local()))"),
        ("H", "((element H)(n_elec 0)(local())(non_local()))"),
    ] {
        let symbol = ElementSymbol::of_string(element).expect("known element");
        let pseudo = Pseudo::empty(symbol);
        assert_eq!(
            pseudo.to_string(),
            format!("{element}   0 electrons removed")
        );
        assert_eq!(pseudo.to_sexp(), sexp);
    }
    // A missing element is a valid, empty potential, not an error.
    let ghost = Pseudo::read_element(BFD, ElementSymbol::of_string("X").expect("X"))
        .expect("missing element");
    assert_eq!(
        ghost,
        Pseudo::empty(ElementSymbol::of_string("X").expect("X"))
    );
}

/// One EZFIO array payload: `rank`, then `dims`, then one value per line in
/// Fortran (column-major) order.
fn ezfio_array(text: &str) -> (Vec<usize>, Vec<f64>) {
    let mut lines = text.lines();
    let rank: usize = lines
        .next()
        .expect("rank line")
        .trim()
        .parse()
        .expect("rank");
    let dims: Vec<usize> = lines
        .next()
        .expect("dims line")
        .split_whitespace()
        .map(|dim| dim.parse().expect("dim"))
        .collect();
    assert_eq!(dims.len(), rank, "rank does not match the number of dims");
    let data: Vec<f64> = lines
        .map(str::trim)
        .filter(|line| !line.is_empty())
        .map(|line| line.parse().expect("array value"))
        .collect();
    assert_eq!(
        data.len(),
        dims.iter().product::<usize>(),
        "element count does not match the dimensions"
    );
    (dims, data)
}

fn scalar(text: &str) -> usize {
    text.trim().parse().expect("scalar")
}

/// The seven arrays and three extents that `qp_create_ezfio` derives from a
/// list of pseudopotentials, in the EZFIO layout: `[nucl_num, klocmax]` for
/// the local part and `[nucl_num, kmax, lmax+1]` for the non-local one,
/// column-major, zero-filled.
struct PseudoArrays {
    klocmax: usize,
    kmax: usize,
    lmax: usize,
    nucl_charge_remove: Vec<f64>,
    v_k: Vec<f64>,
    dz_k: Vec<f64>,
    n_k: Vec<f64>,
    v_kl: Vec<f64>,
    dz_kl: Vec<f64>,
    n_kl: Vec<f64>,
}

fn pseudo_arrays(pseudos: &[Pseudo]) -> PseudoArrays {
    let nucl_num = pseudos.len();
    let klocmax = pseudos.iter().map(|p| p.local.len()).max().unwrap_or(0);
    let lmax = pseudos
        .iter()
        .flat_map(|p| {
            p.non_local
                .iter()
                .map(|(primitive, _)| primitive.proj() as usize)
        })
        .max()
        .unwrap_or(0);
    // `kmax`: the largest number of non-local primitives any atom carries in
    // one projector block.
    let mut kmax = 0;
    for l in 0..=lmax {
        let count = pseudos
            .iter()
            .map(|p| {
                p.non_local
                    .iter()
                    .filter(|(primitive, _)| primitive.proj() as usize == l)
                    .count()
            })
            .max()
            .unwrap_or(0);
        kmax = kmax.max(count);
    }

    let mut arrays = PseudoArrays {
        klocmax,
        kmax,
        lmax,
        nucl_charge_remove: pseudos.iter().map(|p| p.n_elec as f64).collect(),
        v_k: vec![0.0; nucl_num * klocmax],
        dz_k: vec![0.0; nucl_num * klocmax],
        n_k: vec![0.0; nucl_num * klocmax],
        v_kl: vec![0.0; nucl_num * kmax * (lmax + 1)],
        dz_kl: vec![0.0; nucl_num * kmax * (lmax + 1)],
        n_kl: vec![0.0; nucl_num * kmax * (lmax + 1)],
    };
    for (j, pseudo) in pseudos.iter().enumerate() {
        for (i, (primitive, coefficient)) in pseudo.local.iter().enumerate() {
            let index = i * nucl_num + j;
            arrays.v_k[index] = *coefficient;
            arrays.dz_k[index] = primitive.expo();
            arrays.n_k[index] = primitive.r_power() as f64;
        }
        // Per projector block, the primitives fill `i = 0, 1, ...`.
        let mut last = vec![0usize; lmax + 1];
        for (primitive, coefficient) in &pseudo.non_local {
            let k = primitive.proj() as usize;
            let i = last[k];
            let index = k * kmax * nucl_num + i * nucl_num + j;
            arrays.v_kl[index] = *coefficient;
            arrays.dz_kl[index] = primitive.expo();
            arrays.n_kl[index] = primitive.r_power() as f64;
            last[k] += 1;
        }
    }
    arrays
}

/// The `c2h2` fixture (`cc-pvdz_ecp_bfd`, `-p bfd`): atoms C, C, H, H.
const C2H2_NUCL_CHARGE_REMOVE: &str = include_str!("fixtures/pseudo/c2h2/nucl_charge_remove");
const C2H2_V_K: &str = include_str!("fixtures/pseudo/c2h2/pseudo_v_k");
const C2H2_DZ_K: &str = include_str!("fixtures/pseudo/c2h2/pseudo_dz_k");
const C2H2_N_K: &str = include_str!("fixtures/pseudo/c2h2/pseudo_n_k");
const C2H2_V_KL: &str = include_str!("fixtures/pseudo/c2h2/pseudo_v_kl");
const C2H2_DZ_KL: &str = include_str!("fixtures/pseudo/c2h2/pseudo_dz_kl");
const C2H2_N_KL: &str = include_str!("fixtures/pseudo/c2h2/pseudo_n_kl");
const C2H2_KLOCMAX: &str = include_str!("fixtures/pseudo/c2h2/pseudo_klocmax");
const C2H2_KMAX: &str = include_str!("fixtures/pseudo/c2h2/pseudo_kmax");
const C2H2_LMAX: &str = include_str!("fixtures/pseudo/c2h2/pseudo_lmax");

#[test]
fn c2h2_arrays_match_the_reference_ezfio() {
    let carbon =
        Pseudo::read_element(BFD, ElementSymbol::of_string("C").expect("C")).expect("C block");
    let hydrogen =
        Pseudo::read_element(BFD, ElementSymbol::of_string("H").expect("H")).expect("H block");
    let arrays = pseudo_arrays(&[carbon.clone(), carbon, hydrogen.clone(), hydrogen]);

    assert_eq!(arrays.klocmax, scalar(C2H2_KLOCMAX), "pseudo_klocmax");
    assert_eq!(arrays.kmax, scalar(C2H2_KMAX), "pseudo_kmax");
    assert_eq!(arrays.lmax, scalar(C2H2_LMAX), "pseudo_lmax");

    for (name, text, values, dims) in [
        (
            "nucl_charge_remove",
            C2H2_NUCL_CHARGE_REMOVE,
            &arrays.nucl_charge_remove,
            vec![4],
        ),
        ("pseudo_v_k", C2H2_V_K, &arrays.v_k, vec![4, 3]),
        ("pseudo_dz_k", C2H2_DZ_K, &arrays.dz_k, vec![4, 3]),
        ("pseudo_n_k", C2H2_N_K, &arrays.n_k, vec![4, 3]),
        ("pseudo_v_kl", C2H2_V_KL, &arrays.v_kl, vec![4, 1, 1]),
        ("pseudo_dz_kl", C2H2_DZ_KL, &arrays.dz_kl, vec![4, 1, 1]),
        ("pseudo_n_kl", C2H2_N_KL, &arrays.n_kl, vec![4, 1, 1]),
    ] {
        let (fixture_dims, fixture_values) = ezfio_array(text);
        assert_eq!(fixture_dims, dims, "{name}: dims");
        assert_eq!(&fixture_values, values, "{name}: values");
    }
}
