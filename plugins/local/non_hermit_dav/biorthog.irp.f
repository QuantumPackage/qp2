
! ---

subroutine non_hrmt_bieig(n, A, thr_d, thr_nd, leigvec, reigvec, n_real_eigv, eigval)

  BEGIN_DOC
  ! 
  ! routine which returns the sorted REAL EIGENVALUES ONLY and corresponding LEFT/RIGHT eigenvetors 
  ! of a non hermitian matrix A(n,n)
  !
  ! n_real_eigv is the number of real eigenvalues, which might be smaller than the dimension "n" 
  !
  END_DOC

  implicit none
  integer,          intent(in)  :: n
  double precision, intent(in)  :: A(n,n)
  double precision, intent(in)  :: thr_d, thr_nd
  integer,          intent(out) :: n_real_eigv
  double precision, intent(out) :: reigvec(n,n), leigvec(n,n), eigval(n)

  integer                       :: i, j,k 
  integer                       :: n_good
  double precision              :: thr, thr_cut, thr_diag, thr_norm
  double precision              :: accu_d, accu_nd

  integer,          allocatable :: list_good(:), iorder(:), deg_num(:)
  double precision, allocatable :: WR(:), WI(:), VL(:,:), VR(:,:)
  double precision, allocatable :: S(:,:)
  double precision, allocatable  :: phi_1_tilde(:),phi_2_tilde(:),chi_1_tilde(:),chi_2_tilde(:)

  allocate(phi_1_tilde(n),phi_2_tilde(n),chi_1_tilde(n),chi_2_tilde(n))

  allocate(WR(n), WI(n), VL(n,n), VR(n,n)) 

  call lapack_diag_non_sym(n, A, WR, WI, VL, VR)

  thr_diag = 1d-06
  thr_norm = 1d+10

  ! ---

  ! track & sort the real eigenvalues 

  n_good = 0
  thr    = Im_thresh_tcscf
  do i = 1, n
    if(dabs(WI(i)) .lt. thr) then
      n_good += 1
    else
      print*, 'Found an imaginary component to eigenvalue on i = ', i
      print*, 'Re(i) + Im(i)', WR(i), WI(i)
    endif
  enddo

  if(n_good.ne.n) then
    print*,'there are some imaginary eigenvalues '
    thr_diag = 1d-03
    n_good = n
  endif

  allocate(list_good(n_good), iorder(n_good))

  n_good = 0
  do i = 1, n
    n_good += 1
    list_good(n_good) = i
    eigval(n_good) = WR(i)
  enddo

  deallocate( WR, WI )

  n_real_eigv = n_good 
  do i = 1, n_good
    iorder(i) = i
  enddo
  call dsort(eigval, iorder, n_good)
      
  reigvec(:,:) = 0.d0 
  leigvec(:,:) = 0.d0 
  do i = 1, n_real_eigv
    do j = 1, n
      reigvec(j,i) = VR(j,list_good(iorder(i)))
      leigvec(j,i) = VL(j,list_good(iorder(i)))
    enddo
  enddo

  deallocate( list_good, iorder )
  deallocate( VL, VR )

  ASSERT(n==n_real_eigv)

  ! ---

  ! check bi-orthogonality

  thr_diag = 10.d0
  thr_norm = 1d+10

  allocate( S(n_real_eigv,n_real_eigv) )
  call check_biorthog(n, n_real_eigv, leigvec, reigvec, accu_d, accu_nd, S, thr_d, thr_nd, .false.)

  if( (accu_nd .lt. thr_nd) .and. (dabs(accu_d-dble(n_real_eigv))/dble(n_real_eigv) .lt. thr_d) ) then

    print *, ' lapack vectors are normalized and bi-orthogonalized'
    deallocate(S)
    return

  ! accu_nd is modified after adding the normalization
  elseif( (accu_nd .lt. thr_nd) .and. (dabs(accu_d-dble(n_real_eigv))/dble(n_real_eigv) .gt. thr_d) ) then

    print *, ' lapack vectors are not normalized but bi-orthogonalized'
    call check_biorthog_binormalize(n, n_real_eigv, leigvec, reigvec, thr_d, thr_nd, .true.)

    call check_biorthog(n, n_real_eigv, leigvec, reigvec, accu_d, accu_nd, S, thr_d, thr_nd, .true.)
    call check_EIGVEC(n, n, A, eigval, leigvec, reigvec, thr_diag, thr_norm, .true.)

    deallocate(S)
    return

  else

    print *, ' lapack vectors are not normalized neither bi-orthogonalized'

    allocate(deg_num(n))
    call reorder_degen_eigvec(n, deg_num, eigval, leigvec, reigvec)
    call impose_biorthog_degen_eigvec(n, deg_num, eigval, leigvec, reigvec)
    deallocate(deg_num)

    call check_biorthog(n, n_real_eigv, leigvec, reigvec, accu_d, accu_nd, S, thr_d, thr_nd, .false.)
    if( (accu_nd .lt. thr_nd) .and. (dabs(accu_d-dble(n_real_eigv)) .gt. thr_d) ) then
      call check_biorthog_binormalize(n, n_real_eigv, leigvec, reigvec, thr_d, thr_nd, .true.)
    endif
    call check_biorthog(n, n_real_eigv, leigvec, reigvec, accu_d, accu_nd, S, thr_d, thr_nd, .true.)

    deallocate(S)

  endif

  return

end

! ---

subroutine check_bi_ortho(reigvec, leigvec, n, S, accu_nd)

  BEGIN_DOC
  ! retunrs the overlap matrix S = Leigvec^T Reigvec 
  !
  ! and the square root of the sum of the squared off-diagonal elements of S
  END_DOC

  implicit none
  integer,          intent(in)  :: n
  double precision, intent(in)  :: reigvec(n,n), leigvec(n,n)
  double precision, intent(out) :: S(n,n), accu_nd

  integer :: i,j

  ! S = VL x VR
  call dgemm( 'T', 'N', n, n, n, 1.d0 &
            , leigvec, size(leigvec, 1), reigvec, size(reigvec, 1)  &
            , 0.d0, S, size(S, 1) )

  accu_nd = 0.d0
  do i = 1, n
    do j = 1, n
      if(i.ne.j) then
        accu_nd = accu_nd + S(j,i) * S(j,i)
      endif
    enddo
  enddo
  accu_nd = dsqrt(accu_nd)

end


