program tc_rpa

  BEGIN_DOC
  !
  !
  !
  END_DOC

  my_grid_becke = .True.
  PROVIDE tc_grid1_a tc_grid1_r
  my_n_pt_r_grid = tc_grid1_r
  my_n_pt_a_grid = tc_grid1_a
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  if(j1b_type .ge. 100) then
    my_extra_grid_becke  = .True.
    PROVIDE tc_grid2_a tc_grid2_r
    my_n_pt_r_extra_grid = tc_grid2_r
    my_n_pt_a_extra_grid = tc_grid2_a
    touch my_extra_grid_becke my_n_pt_r_extra_grid my_n_pt_a_extra_grid

    call write_int(6, my_n_pt_r_extra_grid, 'radial  internal grid over')
    call write_int(6, my_n_pt_a_extra_grid, 'angular internal grid over')
  endif


  call main()

end

! ---

subroutine main()

  implicit none
  integer                       :: i, j, n
  integer                       :: n_good, n_real_eigv
  double precision              :: thr_cpx, thr_d, thr_nd
  double precision              :: accu_d, accu_nd
  integer,          allocatable :: list_good(:), iorder(:)
  double precision, allocatable :: WR(:), WI(:), VL(:,:), VR(:,:)
  double precision, allocatable :: Omega_p(:), Reigvec_p(:,:), Leigvec_p(:,:)
  double precision, allocatable :: Omega_m(:), Reigvec_m(:,:), Leigvec_m(:,:)
  double precision, allocatable :: S(:,:)

  PROVIDE M_RPA

  print *, ' '
  print *, ' Computing left/right eigenvectors for TC-RPA ...'
  print *, ' '
 

  n = 2 * nS_exc

  thr_cpx = 1d-7
  thr_d   = 1d-07
  thr_nd  = 1d-07


  allocate(WR(n), WI(n), VL(n,n), VR(n,n))  
  call lapack_diag_non_sym(n, M_RPA, WR, WI, VL, VR)
  FREE M_RPA

  print *, ' excitation energies:'
  do i = 1, nS_exc
    write(*, '(I3, X, 1000(F16.10,X))') i, WR(i), WI(i)
    if(dabs(WI(i)) .gt. thr_cpx) then
      print *, ' WARNING ! IMAGINARY EIGENVALUES !!!'
      write(*, '(1000(F16.10,X))') WR(i), WI(i+1)
    endif
  enddo

  print *, ' '
  print *, ' desexcitation energies:'
  do i = nS_exc+1, n
    write(*, '(I3, X, 1000(F16.10,X))') i, WR(i), WI(i)
    if(dabs(WI(i)) .gt. thr_cpx) then
      print *, ' WARNING ! IMAGINARY EIGENVALUES !!!'
      write(*, '(1000(F16.10,X))') WR(i), WI(i+1)
    endif
  enddo


  ! track & sort the real eigenvalues 

  n_good = 0
  do i = 1, nS_exc
    if(dabs(WI(i)) .lt. thr_cpx) then
      if(dabs(WI(nS_exc+i)) .lt. thr_cpx) then
        n_good += 1
      endif
    endif
  enddo
  n_real_eigv = n_good 

  print *, ' '
  print *, ' nb of real eigenvalues  = ', n_real_eigv
  print *, ' total nb of eigenvalues = ', nS_exc

  allocate(Omega_p(n_real_eigv), Reigvec_p(n,n_real_eigv), Leigvec_p(n,n_real_eigv))
  allocate(Omega_m(n_real_eigv), Reigvec_m(n,n_real_eigv), Leigvec_m(n,n_real_eigv))

  n_good = 0
  do i = 1, nS_exc
    if(dabs(WI(i)) .lt. thr_cpx) then
      if(dabs(WI(nS_exc+i)) .lt. thr_cpx) then
        n_good += 1

        Omega_p(n_good) = WR(i)
        do j = 1, n
          Reigvec_p(j,n_good) = VR(j,n_good)
          Leigvec_p(j,n_good) = VL(j,n_good)
        enddo

        Omega_m(n_good) = WR(nS_exc+i)
        do j = 1, n
          Reigvec_m(j,n_good) = VR(j,nS_exc+n_good)
          Leigvec_m(j,n_good) = VL(j,nS_exc+n_good)
        enddo
      endif
    endif
  enddo

  deallocate(WR, WI, VL, VR)  


  ! check bi-orthogonality

  ! first block

  allocate(S(n_real_eigv,n_real_eigv))

  call check_biorthog(n, n_real_eigv, Leigvec_p, Reigvec_p, accu_d, accu_nd, S, thr_d, thr_nd, .false.)
  print *, ' accu_d  = ', accu_d
  print *, ' accu_nd = ', accu_nd

  if((accu_nd .lt. thr_nd) .and. (dabs(accu_d-dble(n_real_eigv))/dble(n_real_eigv) .lt. thr_d)) then
    print *, ' RPA first-block eigenvectors are normalized and bi-orthogonalized'
  else
    print *, ' RPA first-block eigenvectors are neither normalized nor bi-orthogonalized'

    call reorder_degen_eigvec(n, Omega_p, Leigvec_p, Reigvec_p)
    call impose_biorthog_degen_eigvec(n, Omega_p, Leigvec_p, Reigvec_p)

    call check_biorthog(n, n_real_eigv, Leigvec_p, Reigvec_p, accu_d, accu_nd, S, thr_d, thr_nd, .false.)
    if( (accu_nd .lt. thr_nd) .and. (dabs(accu_d-dble(n_real_eigv)) .gt. thr_d) ) then
      call check_biorthog_binormalize(n, n_real_eigv, Leigvec_p, Reigvec_p, thr_d, thr_nd, .true.)
    endif
    call check_biorthog(n, n_real_eigv, Leigvec_p, Reigvec_p, accu_d, accu_nd, S, thr_d, thr_nd, .true.)
  endif


  ! second block

  call check_biorthog(n, n_real_eigv, Leigvec_m, Reigvec_m, accu_d, accu_nd, S, thr_d, thr_nd, .false.)
  print *, ' accu_d  = ', accu_d
  print *, ' accu_nd = ', accu_nd

  if((accu_nd .lt. thr_nd) .and. (dabs(accu_d-dble(n_real_eigv))/dble(n_real_eigv) .lt. thr_d)) then
    print *, ' RPA first-block eigenvectors are normalized and bi-orthogonalized'
  else
    print *, ' RPA first-block eigenvectors are neither normalized nor bi-orthogonalized'

    call reorder_degen_eigvec(n, Omega_m, Leigvec_m, Reigvec_m)
    call impose_biorthog_degen_eigvec(n, Omega_m, Leigvec_m, Reigvec_m)

    call check_biorthog(n, n_real_eigv, Leigvec_m, Reigvec_m, accu_d, accu_nd, S, thr_d, thr_nd, .false.)
    if( (accu_nd .lt. thr_nd) .and. (dabs(accu_d-dble(n_real_eigv)) .gt. thr_d) ) then
      call check_biorthog_binormalize(n, n_real_eigv, Leigvec_m, Reigvec_m, thr_d, thr_nd, .true.)
    endif
    call check_biorthog(n, n_real_eigv, Leigvec_m, Reigvec_m, accu_d, accu_nd, S, thr_d, thr_nd, .true.)
  endif

  deallocate(S)

  return

end

! ---

