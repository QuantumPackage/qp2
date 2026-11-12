! The core PT2 drivers, pt2_stochastic and omp_selection, no longer
! call a plugin hook: they are plain routines in module cipsi_utils
! and cannot know about the transcorrelated data.  This subroutine
! does what the old distributed hook (the slave-side PT2 provider)
! did, and each call site in this plugin calls it just before
! pt2_stochastic or omp_selection.
subroutine provide_for_tc_pt2
  PROVIDE psi_selectors_coef_transp_tc psi_det_sorted_tc psi_det_sorted_tc_order
  PROVIDE dominant_dets_of_cfgs N_dominant_dets_of_cfgs
end
