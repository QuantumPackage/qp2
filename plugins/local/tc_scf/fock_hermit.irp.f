
! ---

BEGIN_PROVIDER [ double precision, good_hermit_tc_fock_mat, (mo_num, mo_num)]

  BEGIN_DOC
! good_hermit_tc_fock_mat = Hermitian Upper triangular Fock matrix 
!
! The converged eigenvectors of such matrix yield to orthonormal vectors satisfying the left Brillouin theorem
  END_DOC
  implicit none
  integer :: i, j

  good_hermit_tc_fock_mat = Fock_matrix_tc_mo_tot
  do j = 1, mo_num
    do i = 1, j-1
      good_hermit_tc_fock_mat(i,j) = Fock_matrix_tc_mo_tot(j,i) 
    enddo
  enddo

END_PROVIDER 

BEGIN_PROVIDER [ double precision, hermit_average_tc_fock_mat, (mo_num, mo_num)]

  BEGIN_DOC
! hermit_average_tc_fock_mat = (F + F^\dagger)/2
  END_DOC
  implicit none
  integer :: i, j

  hermit_average_tc_fock_mat = Fock_matrix_tc_mo_tot
  do j = 1, mo_num
    do i = 1, mo_num
      hermit_average_tc_fock_mat(i,j) = 0.5d0 * (Fock_matrix_tc_mo_tot(j,i) + Fock_matrix_tc_mo_tot(i,j))
    enddo
  enddo

END_PROVIDER 


! ---
BEGIN_PROVIDER [ double precision, grad_hermit]
 implicit none
 BEGIN_DOC
 ! square of gradient of the energy
 END_DOC
 if(symetric_fock_tc)then
  grad_hermit = grad_hermit_average_tc_fock_mat
 else
  grad_hermit = grad_good_hermit_tc_fock_mat
 endif

END_PROVIDER 

BEGIN_PROVIDER [ double precision, grad_good_hermit_tc_fock_mat]
  implicit none
  BEGIN_DOC
  ! grad_good_hermit_tc_fock_mat = norm of gradients of the upper triangular TC fock
  END_DOC
  integer :: i, j
  grad_good_hermit_tc_fock_mat = 0.d0
  do i = 1, elec_alpha_num
    do j = elec_alpha_num+1, mo_num
      grad_good_hermit_tc_fock_mat += dabs(good_hermit_tc_fock_mat(i,j))
    enddo
  enddo
END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, grad_hermit_average_tc_fock_mat]
  implicit none
  BEGIN_DOC
  ! grad_hermit_average_tc_fock_mat = norm of gradients of the upper triangular TC fock
  END_DOC
  integer :: i, j
  grad_hermit_average_tc_fock_mat = 0.d0
  do i = 1, elec_alpha_num
    do j = elec_alpha_num+1, mo_num
      grad_hermit_average_tc_fock_mat += dabs(hermit_average_tc_fock_mat(i,j))
    enddo
  enddo
END_PROVIDER 


! ---

subroutine save_good_hermit_tc_eigvectors()

  implicit none
  integer        :: sign
  character*(64) :: label
  logical        :: output

  sign = 1
  label = "Canonical"
  output = .False.
  
  if(symetric_fock_tc)then
   call mo_as_eigvectors_of_mo_matrix(hermit_average_tc_fock_mat, mo_num, mo_num, label, sign, output)
  else
   call mo_as_eigvectors_of_mo_matrix(good_hermit_tc_fock_mat, mo_num, mo_num, label, sign, output)
  endif
end subroutine save_good_hermit_tc_eigvectors

! ---

