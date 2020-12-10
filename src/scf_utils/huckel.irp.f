subroutine huckel_guess
  implicit none
  BEGIN_DOC
! Build the MOs using the extended Huckel model
  END_DOC
  integer                        :: i,j
  double precision               :: accu
  double precision               :: c
  character*(64)                 :: label
  double precision, allocatable  :: A(:,:)
  label = "Guess"
  c = 0.5d0 * 1.75d0

  allocate (A(ao_num, ao_num))
  A = 0.d0
  do j=1,ao_num
    do i=1,ao_num
      A(i,j) = c * ao_overlap(i,j) * &
         (ao_one_e_integrals_diag(i) + ao_one_e_integrals_diag(j)  ) 
    enddo
    A(j,j) = ao_one_e_integrals_diag(j) + ao_two_e_integral_alpha(j,j) 
  enddo

  Fock_matrix_ao_alpha(1:ao_num,1:ao_num) = A(1:ao_num,1:ao_num)
  Fock_matrix_ao_beta (1:ao_num,1:ao_num) = A(1:ao_num,1:ao_num)

  TOUCH Fock_matrix_ao_alpha Fock_matrix_ao_beta
  mo_coef = eigenvectors_fock_matrix_mo
  call restore_symmetry(ao_num,mo_num,mo_coef,size(mo_coef,1),1.d-10)
  call orthonormalize_mos
  SOFT_TOUCH mo_coef
  call save_mos
  deallocate(A)

end
