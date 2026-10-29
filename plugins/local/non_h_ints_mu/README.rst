=============
non_h_ints_mu
=============

The non-Hermitian part of the mu-transcorrelated Hamiltonian.

This module computes the non-Hermitian potential of the mu-TC
Hamiltonian on the |AO| and the BI-ORTHO |MO| basis. The operator is
defined in Eq. 33 of JCP 154, 084119 (2021).

The two providers are:

* `ao_non_hermit_term_chemist`: the non-Hermitian part of the
  two-electron TC Hamiltonian in the |AO| basis
* `mo_non_hermit_term_chemist`: the same operator in the BI-ORTHO
  |MO| basis

The integrals are built from the Jastrow derivative routines
(`jastrow_psi.irp.f`, `jast_deriv*.irp.f`), the mean field
:math:`\mu(r)` and its gradient (`mu_of_r.irp.f`), and the one- and
two-body J-mu integral kernels, in the bump (`j_bump.irp.f`) and
Gaussian (`deb_j_gauss.irp.f`) forms. The `deb_*.irp.f` files are
debug programs for the individual pieces.
