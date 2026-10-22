=====
utils
=====

General purpose utilities.

* units conversions: :c:func:`ha_to_ev`, :c:func:`au_to_d`, and the
  constants :c:data:`light_speed`, :c:data:`planck_cte`
* mathematics: factorials and binomials (:c:func:`fact`,
  :c:func:`logfact`, :c:func:`binom`), the error function and its
  complex variants, the Boys functions (:c:func:`zboysfun` and its
  variants)
* polynomial and Gaussian products, the workhorse of the integral
  computations: :c:func:`multiply_poly` and its variants,
  :c:func:`give_explicit_poly_and_gaussian` and its variants,
  :c:func:`recentered_poly2` and its variants,
  :c:func:`overlap_gaussian_x` and its variants
* one dimensional numerical integration: :c:func:`crint_1`,
  :c:func:`crint_2` and their variants, :c:func:`v2_over_x`
* linear algebra: :c:func:`lapack_diag` and its variants,
  :c:func:`svd` and its variants, the orthonormalizations
  :c:func:`ortho_lowdin`, :c:func:`ortho_qr`,
  :c:func:`ortho_canonical`, :c:func:`pivoted_cholesky`,
  :c:func:`randomized_svd`, :c:func:`exp_matrix`,
  :c:func:`get_inverse` and its variants
* maps for the sparse storage of the integrals
  (:c:func:`map_load_from_disk`, :c:func:`map_save_to_disk`)
* the Lebedev grids for the angular integration
  (:c:data:`phi_angular_integration_lebedev`,
  :c:data:`theta_angular_integration_lebedev`,
  :c:data:`weights_angular_integration_lebedev`)
* memory management (:c:func:`check_mem`, :c:func:`total_memory`,
  :c:data:`memory_of_double` and its variants) and sorting
  (:c:func:`set_order` and its variants)
