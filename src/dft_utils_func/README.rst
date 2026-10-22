===============
dft_utils_func
===============

Exchange and correlation functionals of |DFT| and |RSDFT|.

* the LDA exchange and correlation, and their short range
  variants: `ecmd_lda.irp.f` and `routines_exc_sr_lda.irp.f`
  (:c:data:`ex_lda`, :c:data:`ec_lda`, :c:data:`ex_lda_sr`,
  :c:data:`ec_lda_sr`, and the range separated correlation
  :c:data:`ecorrlr` / :c:data:`vcorrlr`)
* the PBE and the SR-PBE functionals: `ecmd_pbe_general.irp.f` and
  `routines_exc_sr_pbe.irp.f` (:c:data:`ex_pbe_sr`,
  :c:data:`ec_pbe_sr`)
* the meta-GGA |PBE| on-top correlation, with its implementation in
  terms of the uniform electron gas: `ecmd_pbe_on_top.irp.f`,
  `ecmd_pbe_ueg.irp.f` and `on_top_from_ueg.irp.f`
  (:c:data:`ecmd_pbe_ueg_at_r`, :c:data:`g0_ueg_mu`,
  :c:data:`mu_correction_of_on_top`)
* the range separation functions of the interaction:
  `mu_of_r_dft.irp.f` (:c:data:`mu_rs_c`,
  :c:data:`mu_grad_rho_func`) and the error function and screened
  Coulomb functions in `utils.irp.f` (:c:data:`g0`, :c:data:`g0f`,
  :c:data:`g0s`, :c:data:`g0d` and their derivatives)
* the transformations between the spin resolved and the total
  density, for the screened Coulomb interaction:
  `rho_ab_to_rho_tot.irp.f` (:c:data:`rho_oc_to_rho_ab`,
  :c:data:`v_rho_oc_to_v_rho_ab` and their variants)
