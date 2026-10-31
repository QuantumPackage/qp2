subroutine provide_for_selection_slave
  BEGIN_DOC
  ! Provides, beyond the values listed in the caller, the
  ! values used by the CIPSI selection slaves: the
  ! norm-sorted wave function (psi_det_sorted,
  ! psi_det_sorted_order) and the transposed coefficients
  ! of the selected determinants
  ! (psi_selectors_coef_transp).
  END_DOC
  PROVIDE psi_det_sorted_order
  PROVIDE psi_selectors_coef_transp psi_det_sorted 
end

