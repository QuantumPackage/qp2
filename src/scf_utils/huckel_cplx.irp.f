subroutine huckel_guess_complex
  implicit none
  BEGIN_DOC
! Build the MOs using the extended Huckel model
  END_DOC
  integer                        :: i,j
  double precision               :: accu
  double precision               :: c
  character*(64)                 :: label
  complex*16, allocatable  :: A(:,:)
  label = "Guess"
  c = 0.5d0 * 1.75d0

  allocate (A(ao_num, ao_num))
  A = 0.d0
  do j=1,ao_num
    do i=1,ao_num
      A(i,j) = c * ao_overlap_complex(i,j) * (ao_one_e_integrals_diag_complex(i) + ao_one_e_integrals_diag_complex(j))
    enddo
    A(j,j) = ao_one_e_integrals_diag_complex(j) + dble(ao_two_e_integral_alpha_complex(j,j))
    if (dabs(dimag(ao_two_e_integral_alpha_complex(j,j))) .gt. 1.0d-10) then
      stop 'diagonal elements of ao_two_e_integral_alpha should be real'
    endif
  enddo

!  Fock_matrix_ao_alpha(1:ao_num,1:ao_num) = A(1:ao_num,1:ao_num)
!  Fock_matrix_ao_beta (1:ao_num,1:ao_num) = A(1:ao_num,1:ao_num)
  call zlacpy('X', ao_num, ao_num, A, size(A,1), &
         Fock_matrix_ao_alpha_complex, size(Fock_matrix_ao_alpha_complex,1))
  call zlacpy('X', ao_num, ao_num, A, size(A,1), &
         Fock_matrix_ao_beta_complex,  size(Fock_matrix_ao_beta_complex, 1))
  

!  TOUCH mo_coef

  TOUCH Fock_matrix_ao_alpha_complex Fock_matrix_ao_beta_complex
  mo_coef_complex = eigenvectors_fock_matrix_mo_complex
  SOFT_TOUCH mo_coef_complex
  call save_mos
  deallocate(A)

end
!============================================!
!                                            !
!                    kpts                    !
!                                            !
!============================================!
subroutine huckel_guess_kpts
  implicit none
  BEGIN_DOC
! Build the MOs using the extended Huckel model
  END_DOC
  integer                        :: i,j,k
  double precision               :: accu
  double precision               :: c
  character*(64)                 :: label
  complex*16, allocatable  :: A(:,:)
  label = "Guess"
  c = 0.5d0 * 1.75d0

  allocate (A(ao_num, ao_num))
  do k=1,kpt_num
    A = (0.d0,0.d0)
    do j=1,ao_num_per_kpt
      do i=1,ao_num_per_kpt
        A(i,j) = c * ao_overlap_kpts(i,j,k) * (ao_one_e_integrals_diag_kpts(i,k) + ao_one_e_integrals_diag_kpts(j,k))
      enddo
      A(j,j) = ao_one_e_integrals_diag_kpts(j,k) + dble(ao_two_e_integral_alpha_kpts(j,j,k))
      if (dabs(dimag(ao_two_e_integral_alpha_kpts(j,j,k))) .gt. 1.0d-10) then
        stop 'diagonal elements of ao_two_e_integral_alpha should be real'
      endif
    enddo

!    Fock_matrix_ao_alpha(1:ao_num,1:ao_num) = A(1:ao_num,1:ao_num)
!    Fock_matrix_ao_beta (1:ao_num,1:ao_num) = A(1:ao_num,1:ao_num)
    call zlacpy('X', ao_num_per_kpt, ao_num_per_kpt, A, size(A,1), &
           Fock_matrix_ao_alpha_kpts(:,:,k), size(Fock_matrix_ao_alpha_kpts,1))
    call zlacpy('X', ao_num_per_kpt, ao_num_per_kpt, A, size(A,1), &
           Fock_matrix_ao_beta_kpts(:,:,k),  size(Fock_matrix_ao_beta_kpts, 1))
  enddo  

!  TOUCH mo_coef

  !TOUCH fock_matrix_ao_alpha_complex fock_matrix_ao_beta_kpts
  TOUCH fock_matrix_ao_alpha_kpts fock_matrix_ao_beta_kpts
  mo_coef_kpts = eigenvectors_fock_matrix_mo_kpts
  SOFT_TOUCH mo_coef_kpts
  call save_mos
  deallocate(A)

end
subroutine huckel_guess_kpts_real
  implicit none
  BEGIN_DOC
! Build the MOs using the extended Huckel model
  END_DOC
  integer                        :: i,j,k
  double precision               :: accu
  double precision               :: c
  character*(64)                 :: label
  !complex*16, allocatable  :: A(:,:)
  double precision, allocatable  :: A(:,:)
  label = "Guess"
  c = 0.5d0 * 1.75d0

  allocate (A(ao_num_per_kpt, ao_num_per_kpt))
  do k=1,kpt_num
    A = 0.d0
    do j=1,ao_num_per_kpt
      do i=1,ao_num_per_kpt
        A(i,j) = c * ao_overlap_kpts_real(i,j,k) * (ao_one_e_integrals_diag_kpts(i,k) + ao_one_e_integrals_diag_kpts(j,k))
      enddo
      A(j,j) = ao_one_e_integrals_diag_kpts(j,k) + dble(ao_two_e_integral_alpha_kpts(j,j,k))
      if (dabs(dimag(ao_two_e_integral_alpha_kpts(j,j,k))) .gt. 1.0d-10) then
        stop 'diagonal elements of ao_two_e_integral_alpha should be real'
      endif
    enddo

!    Fock_matrix_ao_alpha(1:ao_num,1:ao_num) = A(1:ao_num,1:ao_num)
!    Fock_matrix_ao_beta (1:ao_num,1:ao_num) = A(1:ao_num,1:ao_num)
    call zlacp2('X', ao_num_per_kpt, ao_num_per_kpt, A, size(A,1), &
           Fock_matrix_ao_alpha_kpts(:,:,k), size(Fock_matrix_ao_alpha_kpts,1))
    call zlacp2('X', ao_num_per_kpt, ao_num_per_kpt, A, size(A,1), &
           Fock_matrix_ao_beta_kpts(:,:,k),  size(Fock_matrix_ao_beta_kpts, 1))
    !call zlacpy('X', ao_num_per_kpt, ao_num_per_kpt, A, size(A,1), &
    !       Fock_matrix_ao_alpha_kpts(:,:,k), size(Fock_matrix_ao_alpha_kpts,1))
    !call zlacpy('X', ao_num_per_kpt, ao_num_per_kpt, A, size(A,1), &
    !       Fock_matrix_ao_beta_kpts(:,:,k),  size(Fock_matrix_ao_beta_kpts, 1))
  enddo  

!  TOUCH mo_coef

  !TOUCH fock_matrix_ao_alpha_complex fock_matrix_ao_beta_kpts
  TOUCH fock_matrix_ao_alpha_kpts fock_matrix_ao_beta_kpts
  mo_coef_kpts = eigenvectors_fock_matrix_mo_kpts_real
  SOFT_TOUCH mo_coef_kpts
  call save_mos
  deallocate(A)

end
