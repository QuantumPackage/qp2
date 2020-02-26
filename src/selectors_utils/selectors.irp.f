use bitmasks

BEGIN_PROVIDER [ integer, psi_selectors_size ]
 implicit none
 psi_selectors_size = psi_det_size
END_PROVIDER

BEGIN_PROVIDER [ double precision, psi_selectors_coef_transp, (N_states,psi_selectors_size) ]
  implicit none
  BEGIN_DOC
  ! Transposed psi_selectors
  END_DOC
  integer                        :: i,k

  do i=1,N_det_selectors
    do k=1,N_states
      psi_selectors_coef_transp(k,i) = psi_selectors_coef(i,k)
    enddo
  enddo
END_PROVIDER

BEGIN_PROVIDER [ double precision, psi_selectors_diag_h_mat, (psi_selectors_size) ]
  implicit none
  BEGIN_DOC
  ! Diagonal elements of the H matrix for each selectors
  END_DOC
  integer                        :: i
  double precision :: diag_H_mat_elem
  do i = 1, N_det_selectors
   psi_selectors_diag_h_mat(i) = diag_H_mat_elem(psi_selectors(1,1,i),N_int)
  enddo
END_PROVIDER


BEGIN_PROVIDER [ complex*16, psi_selectors_coef_transp_complex, (N_states,psi_selectors_size) ]
  implicit none
  BEGIN_DOC
  ! Transposed psi_selectors
  END_DOC
  integer                        :: i,k

  do i=1,N_det_selectors
    do k=1,N_states
      psi_selectors_coef_transp_complex(k,i) = psi_selectors_coef_complex(i,k)
    enddo
  enddo
END_PROVIDER

