subroutine extract_ref
  implicit none
  BEGIN_DOC
  ! Replaces the total wave function by the normalized projection on the reference
  END_DOC

  integer                        :: i,j,k
  do k=1,N_states
    do j=1,N_det_ref
      psi_coef(j,k) = psi_ref_coef_normalized(j,k)
    enddo
  enddo

  do j=1,N_det_ref
    do k=1,N_int
      psi_det(k,1,j) = psi_ref(k,1,j)
      psi_det(k,2,j) = psi_ref(k,2,j)
    enddo
  enddo
  N_det = N_det_ref

  call save_wavefunction

end
