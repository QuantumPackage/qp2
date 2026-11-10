# QP2 Rust tool layer

Rust port of the QP2 OCaml tool layer (the 6 `qp_*` executables +
`element_create_db`, the EZfIO client, and the ZMQ task server). The
Fortran/IRPF90 core is untouched; the EZFIO on-disk format is unchanged.

The frozen behavioural contracts are in [`CONTRACTS.md`](CONTRACTS.md); the
design notes are in [`DESIGN-qp_edit-config.md`](DESIGN-qp_edit-config.md).
Every user-visible difference from the OCaml tools, with the reason, is in
[`DEVIATIONS.md`](DEVIATIONS.md).

## Crates

| Crate | Role | Phase |
|---|---|---|
| `qp-ezfio` | QP2 tool layer over `ezfio-core` (the ezfio-rs on-disk format implementation) + typed registry | P1 |
| `qp-input` | `qp_edit` config schema + loader + generic engine + 7 handlers | P2/P4 |
| `qp-types` | Physics-adjacent data types (Atom, Molecule, Gto, Basis, …) | P3 |
| `qp-zmq` | Message protocol + TaskServer | P4 |
| `qp-cli` | Program lookup, `$QP_ROOT` resolution | P4/P5 |

## Build

```
cd rust
cargo build          # dev
cargo test           # unit + golden tests (fixtures in tests/golden, P0)
cargo clippy --all-targets -- -D warnings
```

Style: edition 2024, rustfmt, every public item documented, no `unsafe`
without a written justification. Toolchain pinned in `rust-toolchain.toml`.
