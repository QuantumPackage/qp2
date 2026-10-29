===============
selectors_utils
===============

Quantities derived from the selectors.

From the selectors of the wave function, this module computes:

* the excitation degree of each selector with respect to the |HF|
  reference, :c:data:`exc_degree_per_selectors`, and the list
  :c:data:`double_index_selectors` and the number
  :c:data:`n_double_selectors` of the double excitations among them
* the correlation energy per selector,
  :c:data:`e_corr_per_selectors`: for a double excitation
  :math:`D_i`, it is
  :math:`\langle D_i | H | \text{HF} \rangle c(D_i)/c(\text{HF})`
* the second order estimate of the correlation energy,
  :c:data:`e_corr_second_order`, and its variant restricted to the
  double excitations, :c:data:`e_corr_double_only`
* the diagonal elements of the Hamiltonian matrix for each selector,
  :c:data:`psi_selectors_diag_h_mat`
