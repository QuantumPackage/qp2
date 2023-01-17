subroutine non_hrmt_diag_split_degen(n, A, leigvec, reigvec, n_real_eigv, eigval)

  BEGIN_DOC
  !
  ! routine which returns the sorted REAL EIGENVALUES ONLY and corresponding LEFT/RIGHT eigenvetors 
  !
  ! of a non hermitian matrix A(n,n)
  !
  ! n_real_eigv is the number of real eigenvalues, which might be smaller than the dimension "n" 
  !
  END_DOC

  implicit none

  integer,          intent(in)  :: n
  double precision, intent(in)  :: A(n,n)
  integer,          intent(out) :: n_real_eigv
  double precision, intent(out) :: reigvec(n,n), leigvec(n,n), eigval(n)
  double precision, allocatable :: reigvec_tmp(:,:), leigvec_tmp(:,:)

  integer                       :: i, j, n_degen,k , iteration
  integer                       :: n_good
  double precision              :: shift,shift_current
  double precision              :: r,thr
  integer,          allocatable :: list_good(:), iorder_origin(:),iorder(:)
  double precision, allocatable :: WR(:), WI(:), Vl(:,:), VR(:,:),S(:,:)
  double precision, allocatable :: Aw(:,:),diag_elem(:),A_save(:,:)
  double precision, allocatable :: im_part(:),re_part(:)


  print*,'Computing the left/right eigenvectors ...'
  print*,'Using the degeneracy splitting algorithm'


  ! pre-processing the matrix :: sorting by diagonal elements
  allocate(reigvec_tmp(n,n), leigvec_tmp(n,n))
  allocate(diag_elem(n),iorder_origin(n),A_save(n,n))
  do i = 1, n
   iorder_origin(i) = i
   diag_elem(i) = A(i,i)
  enddo
  call dsort(diag_elem, iorder_origin, n)
  do i = 1, n
   do j = 1, n
    A_save(j,i) = A(iorder_origin(j),iorder_origin(i))
   enddo
  enddo

  shift = 1.d-15
  shift_current = shift
  iteration = 1 
  logical :: good_ortho
  good_ortho = .False.
  do while(n_real_eigv.ne.n.or. .not.good_ortho)
   if(shift.gt.1.d-3)then
    print*,'shift > 1.d-3 !!'
    print*,'Your matrix intrinsically contains complex eigenvalues'
    stop
   endif
   print*,'***** iteration = ',iteration
   print*,'shift = ',shift
   allocate(WR(n), WI(n), VL(n,n), VR(n,n), Aw(n,n))
   Aw = A_save
   do i = 1, n
    do j = 1, n
     if(dabs(Aw(j,i)).lt.shift)then
      Aw(j,i) = 0.d0
     endif
    enddo
   enddo
   call lapack_diag_non_sym(n,Aw,WR,WI,VL,VR)
   allocate(im_part(n),iorder(n))
   do i = 1, n
    im_part(i) = -dabs(WI(i))
    iorder(i) = i
   enddo
   call dsort(im_part, iorder, n)

   shift_current = max(10.d0 * dabs(im_part(1)),shift)
   print*,'Largest imaginary part found in eigenvalues = ',im_part(1)
   print*,'Splitting the degeneracies by ',shift_current
   Aw = A_save
   call split_matrix_degen(Aw,n,shift_current)
   deallocate( im_part, iorder )
   call lapack_diag_non_sym(n,Aw,WR,WI,VL,VR)
   ! You track the real eigenvalues 
   n_good = 0
   do i = 1, n
     if(dabs(WI(i)).lt.1.d-20)then
       n_good += 1
     else
       print*,'Found an imaginary component to eigenvalue'
       print*,'Re(i) + Im(i)',WR(i),WI(i)
     endif
   enddo
   allocate( list_good(n_good), iorder(n_good) )
   n_good = 0
   do i = 1, n
     if(dabs(WI(i)).lt.1.d-20)then
       n_good += 1
       list_good(n_good) = i
       eigval(n_good) = WR(i)
     endif
   enddo
   deallocate( WR, WI )
 
   n_real_eigv = n_good 
   do i = 1, n_good
     iorder(i) = i
   enddo
 
   ! You sort the real eigenvalues 
   call dsort(eigval, iorder, n_good)
 
   reigvec(:,:) = 0.d0 
   leigvec(:,:) = 0.d0 
   do i = 1, n_real_eigv
     do j = 1, n
       reigvec_tmp(j,i) = VR(j,list_good(iorder(i)))
       leigvec_tmp(j,i) = Vl(j,list_good(iorder(i)))
     enddo
   enddo

   if(n_real_eigv == n)then
    allocate(S(n,n))
    call check_bi_ortho(reigvec_tmp,leigvec_tmp,n,S,accu_nd)
    print*,'accu_nd = ',accu_nd
    double precision :: accu_nd
    good_ortho = accu_nd .lt. 1.d-10
    deallocate(S)
   endif
 
   deallocate( list_good, iorder )
   deallocate( VL, VR, Aw)
   shift *= 10.d0
   iteration += 1
  enddo
  do i = 1, n
   do j = 1, n
    reigvec(iorder_origin(j),i) = reigvec_tmp(j,i)
    leigvec(iorder_origin(j),i) = leigvec_tmp(j,i)
   enddo
  enddo

end subroutine non_hrmt_diag_split_degen

! ---

subroutine non_hrmt_real_diag_new(n, A, leigvec, reigvec, n_real_eigv, eigval)

  BEGIN_DOC
  !
  ! routine which returns the sorted REAL EIGENVALUES ONLY and corresponding LEFT/RIGHT eigenvetors 
  !
  ! of a non hermitian matrix A(n,n)
  !
  ! n_real_eigv is the number of real eigenvalues, which might be smaller than the dimension "n" 
  !
  END_DOC

  implicit none

  integer,          intent(in)  :: n
  double precision, intent(in)  :: A(n,n)
  integer,          intent(out) :: n_real_eigv
  double precision, intent(out) :: reigvec(n,n), leigvec(n,n), eigval(n)

  integer                       :: i, j
  integer                       :: n_good
  double precision              :: shift,shift_current
  double precision              :: r,thr
  integer,          allocatable :: list_good(:), iorder(:)
  double precision, allocatable :: WR(:), WI(:), Vl(:,:), VR(:,:)
  double precision, allocatable :: Aw(:,:)
  double precision, allocatable :: im_part(:)


  print*,'Computing the left/right eigenvectors ...'

  ! Eigvalue(n) = WR(n) + i * WI(n)
  shift = 1.d-10
  do while(n_real_eigv.ne.n.or.shift.gt.1.d-3)
   allocate(WR(n), WI(n), VL(n,n), VR(n,n), Aw(n,n))
   Aw = A
   call lapack_diag_non_sym(n,Aw,WR,WI,VL,VR)
   allocate(im_part(n), iorder(n))
   do i = 1, n
    im_part(i) = -dabs(WI(i))
    iorder(i) = i
   enddo
   shift_current = max(10.d0 * dabs(im_part(1)),shift)
   print*,'adding random number of magnitude ',shift_current
   Aw = A
   do i = 1, n
     call RANDOM_NUMBER(r)
     Aw(i,i) += shift_current * r
   enddo
   deallocate( im_part, iorder )
   call lapack_diag_non_sym(n,Aw,WR,WI,VL,VR)
 
   ! You track the real eigenvalues 
   thr = 1.d-10
   n_good = 0
   do i = 1, n
     if(dabs(WI(i)).lt.thr)then
       n_good += 1
     else
       print*,'Found an imaginary component to eigenvalue'
       print*,'Re(i) + Im(i)',WR(i),WI(i)
     endif
   enddo
 
   allocate( list_good(n_good), iorder(n_good) )
   n_good = 0
   do i = 1, n
     if(dabs(WI(i)).lt.thr)then
       n_good += 1
       list_good(n_good) = i
       eigval(n_good) = WR(i)
     endif
   enddo
 
   deallocate( WR, WI )
 
   n_real_eigv = n_good 
   do i = 1, n_good
     iorder(i) = i
   enddo
 
   ! You sort the real eigenvalues 
   call dsort(eigval, iorder, n_good)
 
   reigvec(:,:) = 0.d0 
   leigvec(:,:) = 0.d0 
   do i = 1, n_real_eigv
     do j = 1, n
       reigvec(j,i) = VR(j,list_good(iorder(i)))
       leigvec(j,i) = Vl(j,list_good(iorder(i)))
     enddo
   enddo
 
   deallocate( list_good, iorder )
   deallocate( VL, VR, Aw)
   shift *= 10.d0
  enddo
  if(shift.gt.1.d-3)then
   print*,'shift > 1.d-3 !!'
   print*,'Your matrix intrinsically contains complex eigenvalues'
  endif

end subroutine non_hrmt_real_diag_new

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

  integer                       :: i, j
  integer                       :: n_good
  double precision              :: thr, thr_cut, thr_diag, thr_norm
  double precision              :: accu_d, accu_nd

  integer,          allocatable :: list_good(:), iorder(:)
  double precision, allocatable :: WR(:), WI(:), VL(:,:), VR(:,:)
  double precision, allocatable :: S(:,:)


  ! -------------------------------------------------------------------------------------
  !

  !print *, ' '
  !print *, ' Computing the left/right eigenvectors ...'
  !print *, ' '

  allocate(WR(n), WI(n), VL(n,n), VR(n,n)) 
  
  !print *, ' fock matrix'
  !do i = 1, n
  !  write(*, '(1000(F16.10,X))') A(i,:)
  !enddo

  !thr_cut = 1.d-15
  !call cancel_small_elmts(A, n, thr_cut)

  !call lapack_diag_non_sym_right(n, A, WR, WI, VR)
  call lapack_diag_non_sym(n, A, WR, WI, VL, VR)
  !call lapack_diag_non_sym_new(n, A, WR, WI, VL, VR)

  !print *, ' '
  !print *, ' eigenvalues'
  !do i = 1, n
  !  write(*, '(1000(F16.10,X))') WR(i), WI(i)
  !enddo
  !print *, ' right eigenvect bef' 
  !do i = 1, n
  !  write(*, '(1000(F16.10,X))') VR(:,i)
  !enddo
  !print *, ' left eigenvect bef'
  !do i = 1, n
  !  write(*, '(1000(F16.10,X))') VL(:,i)
  !enddo

  thr_diag = 1d-06
  thr_norm = 1d+10
  call check_EIGVEC(n, n, A, WR, VL, VR, thr_diag, thr_norm, .false.)

  !
  ! -------------------------------------------------------------------------------------

  ! ---

  ! -------------------------------------------------------------------------------------
  !                  track & sort the real eigenvalues 

  n_good = 0
  !thr    = 100d0
  thr    = Im_thresh_tcscf
  do i = 1, n
    !print*, 'Re(i) + Im(i)', WR(i), WI(i)
    if(dabs(WI(i)) .lt. thr) then
      n_good += 1
    else
      print*, 'Found an imaginary component to eigenvalue on i = ', i
      print*, 'Re(i) + Im(i)', WR(i), WI(i)
    endif
  enddo

  if(n_good.ne.n)then
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

  !print *, ' eigenvalues'
  !do i = 1, n
  !  write(*, '(1000(F16.10,X))') eigval(i)
  !enddo
  !print *, ' right eigenvect aft ord' 
  !do i = 1, n
  !  write(*, '(1000(F16.10,X))') reigvec(:,i)
  !enddo
  !print *, ' left eigenvect aft ord'
  !do i = 1, n
  !  write(*, '(1000(F16.10,X))') leigvec(:,i)
  !enddo

  !
  ! -------------------------------------------------------------------------------------

  ! ---

  ! -------------------------------------------------------------------------------------
  !                               check bi-orthogonality

  thr_diag = 10.d0
  thr_norm = 1d+10

  allocate( S(n_real_eigv,n_real_eigv) )
  call check_biorthog(n, n_real_eigv, leigvec, reigvec, accu_d, accu_nd, S, thr_d, thr_nd, .false.)

  if( (accu_nd .lt. thr_nd) .and. (dabs(accu_d-dble(n_real_eigv))/dble(n_real_eigv) .lt. thr_d) ) then

    !print *, ' lapack vectors are normalized and bi-orthogonalized'
    deallocate(S)
    return

  ! accu_nd is modified after adding the normalization
  !elseif( (accu_nd .lt. thr_nd) .and. (dabs(accu_d-dble(n_real_eigv))/dble(n_real_eigv) .gt. thr_d) ) then

  !  print *, ' lapack vectors are not normalized but bi-orthogonalized'
  !  call check_biorthog_binormalize(n, n_real_eigv, leigvec, reigvec, thr_d, thr_nd, .true.)

  !  call check_EIGVEC(n, n, A, eigval, leigvec, reigvec, thr_diag, thr_norm, .true.)

  !  deallocate(S)
  !  return

  else

    !print *, ' lapack vectors are not normalized neither bi-orthogonalized'

    ! ---

!   call impose_orthog_degen_eigvec(n, eigval, reigvec)
!   call impose_orthog_degen_eigvec(n, eigval, leigvec)

    call impose_biorthog_degen_eigvec(n, eigval, leigvec, reigvec)


    !call impose_orthog_biorthog_degen_eigvec(n, thr_d, thr_nd, eigval, leigvec, reigvec)

    !call impose_unique_biorthog_degen_eigvec(n, eigval, mo_coef, ao_overlap, leigvec, reigvec)

    ! ---

    call check_biorthog(n, n_real_eigv, leigvec, reigvec, accu_d, accu_nd, S, thr_d, thr_nd, .false.)
    if( (accu_nd .lt. thr_nd) .and. (dabs(accu_d-dble(n_real_eigv)) .gt. thr_d) ) then
      call check_biorthog_binormalize(n, n_real_eigv, leigvec, reigvec, thr_d, thr_nd, .true.)
    endif
    call check_biorthog(n, n_real_eigv, leigvec, reigvec, accu_d, accu_nd, S, thr_d, thr_nd, .true.)

    !call impose_biorthog_qr(n, n_real_eigv, thr_d, thr_nd, leigvec, reigvec)
    !call impose_biorthog_lu(n, n_real_eigv, thr_d, thr_nd, leigvec, reigvec)

    ! ---

    call check_EIGVEC(n, n, A, eigval, leigvec, reigvec, thr_diag, thr_norm, .true.)

    deallocate(S)

  endif

  !
  ! -------------------------------------------------------------------------------------

  return

end subroutine non_hrmt_bieig

! ---

subroutine non_hrmt_bieig_random_diag(n, A, leigvec, reigvec, n_real_eigv, eigval)

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
  integer,          intent(out) :: n_real_eigv
  double precision, intent(out) :: reigvec(n,n), leigvec(n,n), eigval(n)

  integer                       :: i, j
  integer                       :: n_good
  double precision              :: thr
  double precision              :: accu_nd

  integer,          allocatable :: list_good(:), iorder(:)
  double precision, allocatable :: Aw(:,:)
  double precision, allocatable :: WR(:), WI(:), VL(:,:), VR(:,:)
  double precision, allocatable :: S(:,:)
  double precision :: r


  ! -------------------------------------------------------------------------------------
  !

  print *, 'Computing the left/right eigenvectors ...'
  allocate( WR(n), WI(n), VL(n,n), VR(n,n), Aw(n,n) )

  Aw(:,:) = A(:,:)
  call lapack_diag_non_sym_new(n, Aw, WR, WI, VL, VR)

  thr    = 1.d-12
  double precision, allocatable :: im_part(:)
  n_good = 0
  do i = 1, n
    if( dabs(WI(i)).lt.thr ) then
      n_good += 1
    else
      print*, 'Found an imaginary component to eigenvalue on i = ', i
      print*, 'Re(i) + Im(i)', WR(i), WI(i)
    endif
  enddo
  print*,'n_good = ',n_good
  if(n_good .lt. n)then
   print*,'Removing degeneracies to remove imaginary parts'
   allocate(im_part(n),iorder(n))
   r = 0.d0
   do i = 1, n
     im_part(i) = -dabs(WI(i))
     iorder(i) = i
   enddo
   call dsort(im_part,iorder,n) 
   thr = 10.d0 * dabs(im_part(1))
   print*,'adding random numbers on the diagonal of magnitude ',thr
   Aw(:,:) = A(:,:)
   do i = 1, n
     call RANDOM_NUMBER(r)
     print*,'r = ',r*thr
     Aw(i,i) += thr * r
   enddo
   print*,'Rediagonalizing the matrix with random numbers'
   call lapack_diag_non_sym_new(n, Aw, WR, WI, VL, VR)
   deallocate(im_part,iorder)
  endif
  deallocate( Aw )

  !
  ! -------------------------------------------------------------------------------------

  ! ---

  ! -------------------------------------------------------------------------------------
  !                  track & sort the real eigenvalues 

  n_good = 0
  thr    = 1.d-5
  do i = 1, n
    if( dabs(WI(i)).lt.thr ) then
      n_good += 1
    else
      print*, 'Found an imaginary component to eigenvalue on i = ', i
      print*, 'Re(i) + Im(i)', WR(i), WI(i)
    endif
  enddo
  print*,'n_good = ',n_good
  allocate( list_good(n_good), iorder(n_good) )

  n_good = 0
  do i = 1, n
    if( dabs(WI(i)).lt.thr ) then
      n_good += 1
      list_good(n_good) = i
      eigval(n_good) = WR(i)
    endif
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

  !
  ! -------------------------------------------------------------------------------------

  ! ---

  ! -------------------------------------------------------------------------------------
  !                               check bi-orthogonality

  allocate( S(n_real_eigv,n_real_eigv) )

  ! S = VL x VR
  call dgemm( 'T', 'N', n_real_eigv, n_real_eigv, n, 1.d0          &
            , leigvec, size(leigvec, 1), reigvec, size(reigvec, 1) &
            , 0.d0, S, size(S, 1) )

  accu_nd = 0.d0
  do i = 1, n_real_eigv
    do j = 1, n_real_eigv
      if(i==j) cycle
      accu_nd = accu_nd + S(j,i) * S(j,i)
    enddo
  enddo
  accu_nd = dsqrt(accu_nd)

  if(accu_nd .lt. thresh_biorthog_nondiag) then
    ! L x R is already bi-orthogonal

    print *, ' L & T bi-orthogonality: ok'
    deallocate( S )
    return

  else
    ! impose bi-orthogonality 

    print *, ' L & T bi-orthogonality: not imposed yet'
    print *, ' accu_nd = ', accu_nd
    call impose_biorthog_qr(n, n_real_eigv, thresh_biorthog_diag, thresh_biorthog_nondiag, leigvec, reigvec)
    deallocate( S )
  
  endif

  !
  ! -------------------------------------------------------------------------------------

  return

end subroutine non_hrmt_bieig_random_diag

! ---

subroutine non_hrmt_real_im(n, A, leigvec, reigvec, n_real_eigv, eigval)

  BEGIN_DOC
  ! 
  ! routine which returns the EIGENVALUES sorted the REAL part and corresponding LEFT/RIGHT eigenvetors 
  ! of a non hermitian matrix A(n,n)
  !
  ! n_real_eigv is the number of real eigenvalues, which might be smaller than the dimension "n" 
  !
  END_DOC

  implicit none
  integer,          intent(in)  :: n
  double precision, intent(in)  :: A(n,n)
  integer,          intent(out) :: n_real_eigv
  double precision, intent(out) :: reigvec(n,n), leigvec(n,n), eigval(n)

  integer                       :: i, j
  integer                       :: n_bad
  double precision              :: thr
  double precision              :: accu_nd

  integer,          allocatable :: iorder(:)
  double precision, allocatable :: Aw(:,:)
  double precision, allocatable :: WR(:), WI(:), VL(:,:), VR(:,:)
  double precision, allocatable :: S(:,:)
  double precision :: r

  ! -------------------------------------------------------------------------------------
  !

  print *, 'Computing the left/right eigenvectors ...'
  allocate( WR(n), WI(n), VL(n,n), VR(n,n), Aw(n,n), iorder(n))

  Aw(:,:) = A(:,:)
   do i = 1, n
     call RANDOM_NUMBER(r)
     Aw(i,i) += 10.d-10* r
   enddo
  call lapack_diag_non_sym(n, Aw, WR, WI, VL, VR)

  ! -------------------------------------------------------------------------------------
  !                  track & sort the real eigenvalues 

  i = 1
  thr    = 1.d-15
  n_real_eigv = 0
  do while (i.le.n) 
!    print*,i,dabs(WI(i))
    if( dabs(WI(i)).gt.thr ) then
      print*, 'Found an imaginary component to eigenvalue on i = ', i
      print*, 'Re(i) , Im(i)  ', WR(i), WI(i)
      iorder(i) = i
      eigval(i) = WR(i)
      i+=1
      print*, 'Re(i+1),Im(i+1)',WR(i), WI(i)
      iorder(i) = i
      eigval(i) = WR(i)
      i+=1
    else  
      n_real_eigv += 1
      iorder(i) = i
      eigval(i) = WR(i)
      i+=1
    endif
  enddo
  call dsort(eigval, iorder, n)
  reigvec(:,:) = 0.d0 
  leigvec(:,:) = 0.d0 
  do i = 1, n
    do j = 1, n
      reigvec(j,i) = VR(j,iorder(i))
      leigvec(j,i) = VL(j,iorder(i))
    enddo
  enddo

  deallocate( iorder )
  deallocate( VL, VR )

  !
  ! -------------------------------------------------------------------------------------

  ! ---

  ! -------------------------------------------------------------------------------------
  !                               check bi-orthogonality

  allocate( S(n,n) )

  ! S = VL x VR
  call dgemm( 'T', 'N', n, n, n, 1.d0          &
            , leigvec, size(leigvec, 1), reigvec, size(reigvec, 1) &
            , 0.d0, S, size(S, 1) )

  accu_nd = 0.d0
  do i = 1, n
    do j = 1, n
      if(i==j) cycle
      accu_nd = accu_nd + S(j,i) * S(j,i)
    enddo
  enddo
  accu_nd = dsqrt(accu_nd)

  deallocate( S )

end subroutine non_hrmt_real_im

! ---

subroutine non_hrmt_generalized_real_im(n, A, B, leigvec, reigvec, n_real_eigv, eigval)

  BEGIN_DOC
  ! 
  ! routine which returns the EIGENVALUES sorted the REAL part and corresponding LEFT/RIGHT eigenvetors 
  ! for A R = lambda B R and A^\dagger L = lambda B^\dagger L
  !
  ! n_real_eigv is the number of real eigenvalues, which might be smaller than the dimension "n" 
  !
  END_DOC

  implicit none
  integer,          intent(in)  :: n
  double precision, intent(in)  :: A(n,n),B(n,n)
  integer,          intent(out) :: n_real_eigv
  double precision, intent(out) :: reigvec(n,n), leigvec(n,n), eigval(n)

  integer                       :: i, j
  integer                       :: n_bad
  double precision              :: thr
  double precision              :: accu_nd

  integer,          allocatable :: iorder(:)
  double precision, allocatable :: Aw(:,:),Bw(:,:)
  double precision, allocatable :: WR(:), WI(:), VL(:,:), VR(:,:), beta(:)
  double precision, allocatable :: S(:,:)
  double precision :: r

  ! -------------------------------------------------------------------------------------
  !

  print *, 'Computing the left/right eigenvectors ...'
  allocate( WR(n), WI(n), VL(n,n), VR(n,n), Aw(n,n), Bw(n,n),iorder(n),beta(n))

  Aw(:,:) = A(:,:)
  Bw(:,:) = B(:,:)
  call lapack_diag_general_non_sym(n,Aw,Bw,WR,beta,WI,VL,VR)

  ! -------------------------------------------------------------------------------------
  !                  track & sort the real eigenvalues 

  i = 1
  thr    = 1.d-10
  n_real_eigv = 0
  do while (i.le.n) 
    if( dabs(WI(i)).gt.thr ) then
      print*, 'Found an imaginary component to eigenvalue on i = ', i
      print*, 'Re(i) , Im(i)  ', WR(i), WI(i)
      iorder(i) = i
      eigval(i) = WR(i)/(beta(i) + 1.d-10)
      i+=1
      print*, 'Re(i+1),Im(i+1)',WR(i), WI(i)
      iorder(i) = i
      eigval(i) = WR(i)/(beta(i) + 1.d-10)
      i+=1
    else  
      n_real_eigv += 1
      iorder(i) = i
      eigval(i) = WR(i)/(beta(i) + 1.d-10)
      i+=1
    endif
  enddo
  call dsort(eigval, iorder, n)
  reigvec(:,:) = 0.d0 
  leigvec(:,:) = 0.d0 
  do i = 1, n
    do j = 1, n
      reigvec(j,i) = VR(j,iorder(i))
      leigvec(j,i) = VL(j,iorder(i))
    enddo
  enddo

  deallocate( iorder )
  deallocate( VL, VR )

  !
  ! -------------------------------------------------------------------------------------

  ! ---

  ! -------------------------------------------------------------------------------------
  !                               check bi-orthogonality

  allocate( S(n,n) )

  ! S = VL x VR
  call dgemm( 'T', 'N', n, n, n, 1.d0          &
            , leigvec, size(leigvec, 1), reigvec, size(reigvec, 1) &
            , 0.d0, S, size(S, 1) )

  accu_nd = 0.d0
  do i = 1, n
    do j = 1, n
      if(i==j) cycle
      accu_nd = accu_nd + S(j,i) * S(j,i)
    enddo
  enddo
  accu_nd = dsqrt(accu_nd)

  deallocate( S )

end subroutine non_hrmt_generalized_real_im

! ---

subroutine non_hrmt_bieig_fullvect(n, A, leigvec, reigvec, n_real_eigv, eigval)

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
  integer,          intent(out) :: n_real_eigv
  double precision, intent(out) :: reigvec(n,n), leigvec(n,n), eigval(n)

  integer                       :: i, j
  integer                       :: n_good
  double precision              :: thr
  double precision              :: accu_nd

  integer,          allocatable :: iorder(:)
  double precision, allocatable :: Aw(:,:)
  double precision, allocatable :: WR(:), WI(:), VL(:,:), VR(:,:)
  double precision, allocatable :: S(:,:)
  double precision, allocatable :: eigval_sorted(:)


  ! -------------------------------------------------------------------------------------
  !

  print *, 'Computing the left/right eigenvectors ...'

  allocate( WR(n), WI(n), VL(n,n), VR(n,n), Aw(n,n) )
  Aw(:,:) = A(:,:)

  call lapack_diag_non_sym_new(n, Aw, WR, WI, VL, VR)

  deallocate( Aw )

  !
  ! -------------------------------------------------------------------------------------

  ! ---

  ! -------------------------------------------------------------------------------------
  !                  track & sort the real eigenvalues 

  allocate( eigval_sorted(n), iorder(n) )

  n_good = 0
  thr    = 1.d-10

  do i = 1, n

    iorder(i) = i
    eigval_sorted(i) = WR(i)

    if(dabs(WI(i)) .gt. thr) then
      print*, ' Found an imaginary component to eigenvalue on i = ', i
      print*, ' Re(i) + Im(i)', WR(i), WI(i)
    else
      n_good += 1
    endif

  enddo

  n_real_eigv = n_good 

  call dsort(eigval_sorted, iorder, n)
      
  reigvec(:,:) = 0.d0 
  leigvec(:,:) = 0.d0 
  do i = 1, n
    eigval(i) = WR(i)
    do j = 1, n
      reigvec(j,i) = VR(j,iorder(i))
      leigvec(j,i) = VL(j,iorder(i))
    enddo
  enddo

  deallocate( eigval_sorted, iorder )
  deallocate( WR, WI )
  deallocate( VL, VR )

  !
  ! -------------------------------------------------------------------------------------

  ! ---

  ! -------------------------------------------------------------------------------------
  !                               check bi-orthogonality

  allocate( S(n,n) )

  ! S = VL x VR
  call dgemm( 'T', 'N', n, n, n, 1.d0                              &
            , leigvec, size(leigvec, 1), reigvec, size(reigvec, 1) &
            , 0.d0, S, size(S, 1) )

  accu_nd = 0.d0
  do i = 1, n
    do j = 1, n
      if(i==j) cycle
      accu_nd = accu_nd + S(j,i) * S(j,i)
    enddo
  enddo
  accu_nd = dsqrt(accu_nd)

  if(accu_nd .lt. thresh_biorthog_nondiag) then
    ! L x R is already bi-orthogonal

    !print *, ' L & T bi-orthogonality: ok'
    deallocate( S )
    return

  else
    ! impose bi-orthogonality 

    !print *, ' L & T bi-orthogonality: not imposed yet'
    !print *, ' accu_nd = ', accu_nd
    call impose_biorthog_qr(n, n, thresh_biorthog_diag, thresh_biorthog_nondiag, leigvec, reigvec)
    deallocate( S )
  
  endif

  !
  ! -------------------------------------------------------------------------------------

  return

end subroutine non_hrmt_bieig_fullvect

! ---


subroutine split_matrix_degen(aw,n,shift)
 implicit none
 BEGIN_DOC
 ! subroutines that splits the degeneracies of a matrix by adding a splitting of magnitude thr * n_degen/2
 !
 ! WARNING !! THE MATRIX IS ASSUMED TO BE PASSED WITH INCREASING DIAGONAL ELEMENTS
 END_DOC
 double precision,intent(inout) :: Aw(n,n)
 double precision,intent(in)    :: shift
 integer, intent(in) :: n
 integer :: i,j,n_degen
 logical :: keep_on
 i=1
 do while(i.lt.n)
  if(dabs(Aw(i,i)-Aw(i+1,i+1)).lt.shift)then
   j=1
   keep_on = .True.
   do while(keep_on)
    if(i+j.gt.n)then
     keep_on = .False.
     exit
    endif
    if(dabs(Aw(i,i)-Aw(i+j,i+j)).lt.shift)then
     j+=1
    else
     keep_on=.False.
     exit
    endif
   enddo
   n_degen = j
   j=0
   keep_on = .True.
   do while(keep_on)
    if(i+j+1.gt.n)then
     keep_on = .False.
     exit
    endif
    if(dabs(Aw(i+j,i+j)-Aw(i+j+1,i+j+1)).lt.shift)then
     Aw(i+j,i+j) += (j-n_degen/2) * shift
     j+=1
    else 
     keep_on = .False.
     exit
    endif
   enddo
   Aw(i+n_degen-1,i+n_degen-1) += (n_degen-1-n_degen/2) * shift
   i+=n_degen
  else 
   i+=1
  endif
 enddo

end

subroutine give_degen(a,n,shift,list_degen,n_degen_list)
 implicit none
 BEGIN_DOC
 ! returns n_degen_list :: the number of degenerated SET of elements (i.e. with |A(i)-A(i+1)| below shift)
 !
 ! for each of these sets, list_degen(1,i) = first degenerate element of the set i, 
 !
 !                         list_degen(2,i) = last degenerate element of the set i.
 END_DOC
 double precision,intent(in) :: A(n)
 double precision,intent(in)    :: shift
 integer, intent(in) :: n
 integer, intent(out) :: list_degen(2,n),n_degen_list
 integer :: i,j,n_degen,k
 logical :: keep_on
 double precision,allocatable :: Aw(:)
 list_degen = -1
 allocate(Aw(n))
 Aw = A
 i=1
 k = 0
 do while(i.lt.n)
  if(dabs(Aw(i)-Aw(i+1)).lt.shift)then
   k+=1
   j=1
   list_degen(1,k) = i
   keep_on = .True.
   do while(keep_on)
    if(i+j.gt.n)then
     keep_on = .False.
     exit
    endif
    if(dabs(Aw(i)-Aw(i+j)).lt.shift)then
     j+=1
    else
     keep_on=.False.
     exit
    endif
   enddo
   n_degen = j
   list_degen(2,k) = list_degen(1,k)-1 + n_degen
   j=0
   keep_on = .True.
   do while(keep_on)
    if(i+j+1.gt.n)then
     keep_on = .False.
     exit
    endif
    if(dabs(Aw(i+j)-Aw(i+j+1)).lt.shift)then
     Aw(i+j) += (j-n_degen/2) * shift
     j+=1
    else 
     keep_on = .False.
     exit
    endif
   enddo
   Aw(i+n_degen-1) += (n_degen-1-n_degen/2) * shift
   i+=n_degen
  else 
   i+=1
  endif
 enddo
 n_degen_list = k

end

subroutine cancel_small_elmts(aw,n,shift)
 implicit none
 BEGIN_DOC
 ! subroutines that splits the degeneracies of a matrix by adding a splitting of magnitude thr * n_degen/2
 !
 ! WARNING !! THE MATRIX IS ASSUMED TO BE PASSED WITH INCREASING DIAGONAL ELEMENTS
 END_DOC
 double precision,intent(inout) :: Aw(n,n)
 double precision,intent(in)    :: shift
 integer, intent(in) :: n
 integer :: i,j
 do i = 1, n
  do j = 1, n
   if(dabs(Aw(j,i)).lt.shift)then
    Aw(j,i) = 0.d0
   endif
  enddo
 enddo
end

subroutine check_bi_ortho(reigvec,leigvec,n,S,accu_nd)
 implicit none
 integer, intent(in) :: n
 double precision,intent(in) :: reigvec(n,n),leigvec(n,n)
 double precision, intent(out) :: S(n,n),accu_nd
 BEGIN_DOC
! retunrs the overlap matrix S = Leigvec^T Reigvec 
!
! and the square root of the sum of the squared off-diagonal elements of S
 END_DOC
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
