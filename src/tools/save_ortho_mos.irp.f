program save_ortho_mos
  implicit none
  BEGIN_DOC
  ! Save orthonormalized MOs in the EZFIO.
  !
  ! This program reads the current MOs, computes the corresponding overlap matrix in the MO basis
  !
  ! and perform a Lowdin orthonormalization : :math:`MO_{new} = S^{-1/2} MO_{guess}`.
  !
  ! Thanks to the Lowdin orthonormalization, the new MOs are the most similar to the guess MOs.
  END_DOC
  call orthonormalize_mos
  mo_label = 'Orthonormalized'
  SOFT_TOUCH mo_label
  call save_mos
end
