! The core PT2 driver, pt2_stochastic, no longer calls a plugin hook: it
! is a plain routine in module cipsi and cannot know about the
! transcorrelated data.  This subroutine does what the old
! provide_for_zmq_pt2 did, and each of the five call sites in this plugin
! calls it just before pt2_stochastic.
subroutine provide_for_tc_pt2
  PROVIDE psi_selectors_coef_transp_tc psi_det_sorted_tc psi_det_sorted_tc_order
end
