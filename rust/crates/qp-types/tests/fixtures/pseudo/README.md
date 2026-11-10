# Pseudo fixtures

Everything here is a **copy** of reference data; nothing is edited. The tests
read only this directory.

## `bfd`, `ccecp`

Verbatim copies of the reference pseudopotential library files
`data/pseudo/bfd` and `data/pseudo/ccecp` at commit `be4b71c1`, taken from the
same tree as `rust/tests/golden/fixtures/`:

| file    | sha256 |
| ------- | ------ |
| `bfd`   | `b22f417637d685011e14dd0875eb3ce9940aae1cde7e88f4fbd925b66de6dafb` |
| `ccecp` | `f5deba29602f9bbc75e32af92ed1599c490c38762e436144cacbbe697e952dc8` |

`ocaml/tests/test_pseudo.ml` reads the Cu block of `bfd`; the Rust
`tests/test_group_b_pseudo_md5.rs` pins the same element plus He, C (`bfd`) and
Mg (`ccecp`), and the C/H blocks against the `c2h2` EZFIO arrays below.

## `c2h2/`

Decompressed payloads of `rust/tests/golden/fixtures/c2h2/pseudo/*.gz`, the
reference EZFIO arrays written by `qp_create_ezfio` for the `c2h2` fixture
(`cc-pvdz_ecp_bfd` + `-p bfd`, atoms C, C, H, H). The golden corpus stores the
gzip container; the payload is the oracle (plan D6), so the test reads the
decompressed text and does not depend on `flate2`.

`qp-ezfio`'s textual array format is:

```
  <rank>
<dim_1> <dim_2> ... <dim_rank>
<value>          (one per line, Fortran/column-major order)
```

The test decodes that layout itself (`read_ezfio_array` in
`tests/test_group_b_pseudo_md5.rs`, ten lines) rather than depending on the
`qp-ezfio` crate.

Provenance hashes of the source `.gz` files:

```
79a8a43aef984976fc0929ca0378f5dae5a81da9787242d3bbcd22e68ee71db1  nucl_charge_remove.gz
a6efa06b8335083936575944a4e6cb83930f786188bcaaf240879c024ca58b5b  pseudo_dz_k.gz
36c6a938d60b6d324cc52047a0658418a70af9f7838b46f20ec3c12f83cb6bf7  pseudo_dz_kl.gz
86e3c1716d23e1a472a808edca1a4bf98890255c521dce4fbd045bedb41ac1a1  pseudo_n_k.gz
9c00ad6729a64786c69557b65e5c10f2935c07b927dc5f97fd9e51bd5c63a08a  pseudo_n_kl.gz
ec65d867a3f2957ef57826ad577a56995798d65a96fd1e30e62ebcd23c97b135  pseudo_v_k.gz
c3bd6af71b3248a041a3706531f59d9e95dedb7ea8b01ee3fca86c14d3bc1cbe  pseudo_v_kl.gz
```

Regenerate with:

```sh
for f in nucl_charge_remove pseudo_dz_k pseudo_v_k pseudo_n_k \
         pseudo_dz_kl pseudo_v_kl pseudo_n_kl; do
  gzip -dc rust/tests/golden/fixtures/c2h2/pseudo/$f.gz > \
      rust/crates/qp-types/tests/fixtures/pseudo/c2h2/$f
done
```
