
! ---

BEGIN_PROVIDER [ double precision, fock_vartc_eigvec_mo, (mo_num, mo_num)]

  implicit none

  integer                       :: i, j
  integer                       :: liwork, lwork, n, info
  integer,          allocatable :: iwork(:)
  double precision, allocatable :: work(:), F(:,:), F_save(:,:)
  double precision, allocatable :: diag(:)

  PROVIDE mo_r_coef
  PROVIDE Fock_matrix_vartc_mo_tot

  allocate( F(mo_num,mo_num), F_save(mo_num,mo_num)  )
  allocate (diag(mo_num) )

  do j = 1, mo_num
    do i = 1, mo_num
      F(i,j) = Fock_matrix_vartc_mo_tot(i,j)
    enddo
  enddo

  ! Insert level shift here
  do i = elec_beta_num+1, elec_alpha_num
    F(i,i) += 0.5d0 * level_shift_tcscf
  enddo
  do i = elec_alpha_num+1, mo_num
    F(i,i) += level_shift_tcscf
  enddo

  n = mo_num
  lwork = 1+6*n + 2*n*n
  liwork = 3 + 5*n

  allocate(work(lwork))
  allocate(iwork(liwork) )

  lwork = -1
  liwork = -1

  F_save = F
  call dsyevd('V', 'U', mo_num, F, size(F, 1), diag, work, lwork, iwork, liwork, info)

  if (info /= 0) then
    print *,  irp_here//' DSYEVD failed : ', info
    stop 1
  endif
  lwork = int(work(1))
  liwork = iwork(1)
  deallocate(iwork)
  deallocate(work)

  allocate(work(lwork))
  allocate(iwork(liwork) )
  call dsyevd('V', 'U', mo_num, F, size(F, 1), diag, work, lwork, iwork, liwork, info)
  deallocate(iwork)

  if (info /= 0) then
    F = F_save
    call dsyev('V', 'L', mo_num, F, size(F, 1), diag, work, lwork, info)

    if (info /= 0) then
      print *,  irp_here//' DSYEV failed : ', info
      stop 1
    endif
  endif

  do i = 1, mo_num
    do j = 1, mo_num
      fock_vartc_eigvec_mo(j,i) = F(j,i)
    enddo
  enddo

  deallocate(work, F, F_save, diag)

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, fock_vartc_eigvec_ao, (ao_num, mo_num)]

  implicit none

  PROVIDE mo_r_coef

  call dgemm( 'N', 'N', ao_num, mo_num, mo_num, 1.d0                                             &
            , mo_r_coef, size(mo_r_coef, 1), fock_vartc_eigvec_mo, size(fock_vartc_eigvec_mo, 1) &
            , 0.d0, fock_vartc_eigvec_ao, size(fock_vartc_eigvec_ao, 1))

END_PROVIDER

! ---

