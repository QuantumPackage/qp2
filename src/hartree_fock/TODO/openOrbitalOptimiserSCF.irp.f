subroutine c_print_callback(handle)
  use OpenOrbitalOptimizer
  use iso_c_binding
  implicit none
  type(c_ptr), intent(in), value :: handle
  integer(c_int64_t) :: iteration_SCF, dim_DIIS
  real(c_double) :: energy_SCF, Delta_energy_SCF, max_error_DIIS
  iteration_SCF = ooo_get_int64(handle, 'iter')
  energy_SCF = ooo_get_double(handle, 'E')
  Delta_energy_SCF = ooo_get_double(handle, 'dE')
  max_error_DIIS = ooo_get_double(handle, 'diis_max_error')
  dim_DIIS = ooo_get_int64(handle, 'history_size')
  write(6,'(I4, 1X, F16.10, 1X, F16.10, 1X, F16.10, 1X, F16.10, 1X, I3)')  &
      iteration_SCF, energy_SCF, Delta_energy_SCF, max_error_DIIS, 0.d0, dim_DIIS
end

BEGIN_PROVIDER [ double precision, partial_occupation, (mo_num) ]
 implicit none
 BEGIN_DOC
 ! Partial occupation numbers for optimal damping algorithm
 END_DOC
 partial_occupation = 0.d0
 integer  :: i
 do i=1,elec_beta_num
   partial_occupation(i) = 2.d0
 enddo
 do i=elec_beta_num+1, elec_alpha_num
   partial_occupation(i) = 1.d0
 enddo
END_PROVIDER


BEGIN_PROVIDER [ double precision, S_mo_coef_guess, (ao_num, mo_num) ]
 implicit none
 BEGIN_DOC
 ! Product S.C where S is the overlap matrix in the AO basis and C the mo_coef_begin_iteration matrix.
 END_DOC

 call dgemm('N','N', ao_num, mo_num, ao_num,                   &
     1.d0, ao_overlap,size(ao_overlap,1),      &
     mo_coef_begin_iteration, size(mo_coef_begin_iteration,1),                                     &
     0.d0, S_mo_coef_guess, size(S_mo_coef_guess,1))

END_PROVIDER


subroutine ao_to_mo_guess(A_ao,LDA_ao,A_mo,LDA_mo)
  implicit none
  BEGIN_DOC
  ! Transform A from the |AO| basis to the |MO| basis
  !
  ! $C^\dagger.A_{ao}.C$
  END_DOC
  integer, intent(in)            :: LDA_ao,LDA_mo
  double precision, intent(in)   :: A_ao(LDA_ao,ao_num)
  double precision, intent(out)  :: A_mo(LDA_mo,mo_num)
  double precision, allocatable  :: T(:,:)

  allocate ( T(ao_num,mo_num) )
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: T

  call dgemm('N','N', ao_num, mo_num, ao_num,                    &
      1.d0, A_ao,LDA_ao,                                             &
      mo_coef_begin_iteration, size(mo_coef_begin_iteration,1),                                      &
      0.d0, T, size(T,1))

  call dgemm('T','N', mo_num, mo_num, ao_num,                &
      1.d0, mo_coef_begin_iteration ,size(mo_coef_begin_iteration,1),                                 &
      T, ao_num,                                                     &
      0.d0, A_mo, size(A_mo,1))

  call restore_symmetry(mo_num,mo_num,A_mo,size(A_mo,1),1.d-15)
  deallocate(T)
end


subroutine c_fock_builder(Norb, C, n, F, Etot)
  use iso_c_binding
  implicit none

  integer(c_int64_t), intent(in), value :: Norb
  real(c_double), intent(in)            :: n(Norb)
  real(c_double), intent(in)            :: C(Norb,Norb)
  real(c_double), intent(out)           :: F(Norb,Norb)
  real(c_double), intent(out)           :: Etot


  if (Norb /= mo_num) then
    stop 'In c_fock_builder, Norb /= mo_num'
  endif

  integer :: i
  do i=1,Norb
    partial_occupation(i) = n(i)
  enddo

  call dgemm('N','N', ao_num, mo_num, mo_num, 1.d0, &
     mo_coef_begin_iteration, size(mo_coef_begin_iteration,1), &
     C, size(C,1), &
     0.d0, mo_coef, size(mo_coef,1))
  TOUCH mo_coef partial_occupation
  call ao_to_mo_guess(Fock_matrix_ao,size(Fock_matrix_ao,1),F,size(F,1))
  Etot = hf_energy
end

subroutine OpenOrbitalOptimiserSCF
  use iso_c_binding
  use OpenOrbitalOptimizer
  implicit none
  double precision :: energy
  type(C_FUNPTR) :: cb, cc
  external :: c_fock_builder, c_print_callback
  double precision, allocatable :: C(:,:)
  integer :: i

  call write_time(6)

  print*,'Energy of the guess = ',SCF_energy
  write(6,'(A4, 1X, A16, 1X, A16, 1X, A16, 1X, A16)')  &
    '====','================','================','================','================'
  write(6,'(A4, 1X, A16, 1X, A16, 1X, A16, 1X, A16)')  &
    '  N ', 'Energy  ', 'Energy diff  ',  'DIIS error  ', 'Level shift   '
  write(6,'(A4, 1X, A16, 1X, A16, 1X, A16, 1X, A16)')  &
    '====','================','================','================','================'

  call initialize_mo_coef_begin_iteration
  allocate(C(mo_num,mo_num))
  C = 0.d0
  do i=1,mo_num
    C(i,i) = 1.d0
  enddo

  cb = c_funloc(c_fock_builder)
  cc = c_funloc(c_print_callback)
  energy = rhf_solve_nosym(cb, cc, mo_num*1_8, elec_num*1_8, C, partial_occupation, dsqrt(thresh_SCF))

  write(6,'(A4, 1X, A16, 1X, A16, 1X, A16, 1X, A16)')  &
    '====','================','================','================','================'
  write(6,*)

  mo_label = 'Canonical'

  call save_mos
  call write_double(6, energy, 'SCF energy')
  call write_time(6)
  SOFT_TOUCH mo_coef partial_occupation
end

