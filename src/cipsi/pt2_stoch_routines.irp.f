subroutine provide_for_zmq_pt2
  BEGIN_DOC
  ! Provides, beyond the values listed in the caller, the
  ! values used by the stochastic PT2 slaves: the
  ! norm-sorted wave function (psi_det_sorted,
  ! psi_det_sorted_order), the transposed coefficients of
  ! the selected determinants (psi_selectors_coef_transp),
  ! and the diagonal of H in the determinant basis (psi_det_hii).
  END_DOC
  PROVIDE psi_selectors_coef_transp psi_det_sorted psi_det_sorted_order psi_det_hii 
end
