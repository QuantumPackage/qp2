subroutine extract_cas
  implicit none
  BEGIN_DOC
  ! Replaces the total wave function by the normalized projection on the CAS.
  END_DOC

  integer                        :: i,j,k
  do k=1,N_states
    do j=1,N_det_generators
      psi_coef(j,k) = psi_coef_generators(j,k)
    enddo
  enddo

  do j=1,N_det_generators
    do k=1,N_int
      psi_det(k,1,j) = psi_det_generators(k,1,j)
      psi_det(k,2,j) = psi_det_generators(k,2,j)
    enddo
  enddo
  N_det = N_det_generators

  SOFT_TOUCH N_det psi_det psi_coef
end
