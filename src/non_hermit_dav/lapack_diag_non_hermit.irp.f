subroutine lapack_diag_non_sym(n, A, WR, WI, VL, VR)

  BEGIN_DOC
  ! You enter with a general non hermitian matrix A(n,n) 
  !
  ! You get out with the real WR and imaginary part WI of the eigenvalues 
  !
  ! Eigvalue(n) = WR(n) + i * WI(n)
  !
  ! And the left VL and right VR eigenvectors 
  !
  ! VL(i,j) = <i|Psi_left(j)>  :: projection on the basis element |i> on the jth left  eigenvector 
  !
  ! VR(i,j) = <i|Psi_right(j)> :: projection on the basis element |i> on the jth right eigenvector 
  !
  ! The real part of the matrix A can be written as A = VR D VL^T
  !
  END_DOC

  implicit none

  integer,          intent(in)  :: n
  double precision, intent(in)  :: A(n,n)
  double precision, intent(out) :: WR(n), WI(n), VL(n,n), VR(n,n)

  integer                       :: lda, ldvl, ldvr, LWORK, INFO
  double precision, allocatable :: Atmp(:,:), WORK(:)

  lda  = n
  ldvl = n
  ldvr = n

  allocate( Atmp(n,n) )
  Atmp(1:n,1:n) = A(1:n,1:n)

  allocate(WORK(1))
  LWORK = -1 ! to ask for the optimal size of WORK
  call dgeev('V', 'V', n, Atmp, lda, WR, WI, VL, ldvl, VR, ldvr, WORK, LWORK, INFO)
  if(INFO.gt.0)then
    print*,'dgeev failed !!',INFO
    stop
  endif
  LWORK = max(int(WORK(1)), 1) ! this is the optimal size of WORK 
  deallocate(WORK)

  allocate(WORK(LWORK))

  ! Actual diagonalization 
  call dgeev('V', 'V', n, Atmp, lda, WR, WI, VL, ldvl, VR, ldvr, WORK, LWORK, INFO)
  if(INFO.ne.0) then
    print*,'dgeev failed !!', INFO
    stop
  endif

  deallocate(Atmp, WORK)

end subroutine lapack_diag_non_sym


subroutine non_sym_diag_inv_right(n,A,leigvec,reigvec,n_real_eigv,eigval)
 implicit none
 BEGIN_DOC
! routine which returns the sorted REAL EIGENVALUES ONLY and corresponding LEFT/RIGHT eigenvetors 
!
! of a non hermitian matrix A(n,n)
!
! n_real_eigv is the number of real eigenvalues, which might be smaller than the dimension "n" 
 END_DOC
 integer, intent(in) :: n
 double precision, intent(in) :: A(n,n)
 double precision, intent(out) :: reigvec(n,n),leigvec(n,n),eigval(n)
 double precision, allocatable :: Aw(:,:)
 integer, intent(out) :: n_real_eigv
 print*,'Computing the left/right eigenvectors ...'
 character*1 :: JOBVL,JOBVR
 JOBVL = "V" ! computes the left  eigenvectors 
 JOBVR = "V" ! computes the right eigenvectors 
 double precision, allocatable :: WR(:),WI(:),Vl(:,:),VR(:,:),S(:,:),inv_reigvec(:,:)
 integer :: i,j
 integer :: n_good
 integer, allocatable :: list_good(:), iorder(:)
 double precision :: thr
 thr = 1.d-10
 ! Eigvalue(n) = WR(n) + i * WI(n)
 allocate(WR(n),WI(n),VL(n,n),VR(n,n),Aw(n,n))
 Aw = A
 do i = 1, n
  do j = i+1, n
   if(dabs(Aw(j,j)-Aw(i,i)).lt.thr)then
     Aw(j,j)+= thr
     Aw(i,i)-= thr
!    if(Aw(j,i) * A(i,j) .lt.0d0  )then
!     if(dabs(Aw(j,i) * A(i,j)).lt.thr**(1.5d0))then
!      print*,Aw(j,j),Aw(i,i)
!      print*,Aw(j,i) , A(i,j)
      Aw(j,i) = 0.d0
      Aw(i,j) = Aw(j,i)
!     endif
!    endif
   endif
  enddo
 enddo
 call lapack_diag_non_sym(n,Aw,WR,WI,VL,VR)
 ! You track the real eigenvalues 
 n_good = 0
! do i = 1, n
!  write(*,'(100(F16.12,X))')A(:,i)
! enddo
 do i = 1, n
  print*,'Im part of lambda = ',dabs(WI(i))
  if(dabs(WI(i)).lt.thr)then
   n_good += 1
  else
   print*,'Found an imaginary component to eigenvalue'
   print*,'Re(i) + Im(i)',WR(i),WI(i)
   write(*,'(100(F10.5,X))')VR(:,i)
   write(*,'(100(F10.5,X))')VR(:,i+1)
   write(*,'(100(F10.5,X))')VL(:,i)
   write(*,'(100(F10.5,X))')VL(:,i+1)
  endif
 enddo
 allocate(list_good(n_good),iorder(n_good))
 n_good = 0
 do i = 1, n
  if(dabs(WI(i)).lt.thr)then
   n_good += 1
   list_good(n_good) = i
   eigval(n_good) = WR(i)
  endif
 enddo
 n_real_eigv = n_good 
 do i = 1, n_good
  iorder(i) = i
 enddo
 ! You sort the real eigenvalues 
 call dsort(eigval,iorder,n_good)
 do i = 1, n_real_eigv
  do j = 1, n
   reigvec(j,i) = VR(j,list_good(iorder(i)))
   leigvec(j,i) = VL(j,list_good(iorder(i)))
  enddo
 enddo
 allocate(inv_reigvec(n_real_eigv,n_real_eigv))
! call get_pseudo_inverse(reigvec,n_real_eigv,n_real_eigv,n_real_eigv,inv_reigvec,n_real_eigv,thr)
! do i = 1, n_real_eigv
!  do j = 1, n
!   leigvec(j,i) = inv_reigvec(i,j)
!  enddo
! enddo
 allocate( S(n_real_eigv,n_real_eigv) )

  ! S = VL x VR
  call dgemm( 'T', 'N', n_real_eigv, n_real_eigv, n_real_eigv, 1.d0                              &
            , leigvec, size(leigvec, 1), reigvec, size(reigvec, 1) &
            , 0.d0, S, size(S, 1) )
   do i = 1,n_real_eigv
    write(*,'(100(F10.5,X))')S(:,i)
   enddo
! call lapack_diag_non_sym(n,S,WR,WI,VL,VR)
! print*,'Eigenvalues of S'
! do i = 1, n
!  print*,WR(i),dabs(WI(i))
! enddo
  call dgemm( 'T', 'N', n_real_eigv, n_real_eigv, n_real_eigv, 1.d0                              &
            , leigvec, size(leigvec, 1), reigvec, size(reigvec, 1) &
            , 0.d0, S, size(S, 1) )
! call get_inv_half_svd(S, n_real_eigv, inv_reigvec)

  double precision :: accu_d,accu_nd
  accu_nd = 0.d0
  accu_d = 0.d0
  do i = 1, n_real_eigv
    do j = 1, n_real_eigv
      if(i==j) then
       accu_d += S(j,i) * S(j,i)
      else
       accu_nd = accu_nd + S(j,i) * S(j,i)
      endif
    enddo
  enddo
  accu_nd = dsqrt(accu_nd)

  print*,'accu_nd = ',accu_nd
  if( accu_nd .lt. 1d-10 ) then
    ! L x R is already bi-orthogonal
    !print *, ' L & T bi-orthogonality: ok'
    return
  else
   print*,'PB with bi-orthonormality!!'
   stop
  endif
end

subroutine lapack_diag_non_sym_new(n, A, WR, WI, VL, VR)

  BEGIN_DOC
  !
  ! You enter with a general non hermitian matrix A(n,n) 
  !
  ! You get out with the real WR and imaginary part WI of the eigenvalues 
  !
  ! Eigvalue(n) = WR(n) + i * WI(n)
  !
  ! And the left VL and right VR eigenvectors 
  !
  ! VL(i,j) = <i|Psi_left(j)>  :: projection on the basis element |i> on the jth left  eigenvector 
  !
  ! VR(i,j) = <i|Psi_right(j)> :: projection on the basis element |i> on the jth right eigenvector 
  !
  END_DOC

  implicit none
  integer,          intent(in)  :: n
  double precision, intent(in)  :: A(n,n)
  double precision, intent(out) :: WR(n), WI(n), VL(n,n), VR(n,n)

  character*1                   :: JOBVL,JOBVR,BALANC,SENSE
  integer                       :: ILO, IHI
  integer                       :: lda, ldvl, ldvr, LWORK, INFO
  double precision              :: ABNRM
  integer,          allocatable :: IWORK(:)
  double precision, allocatable :: WORK(:), SCALE_array(:), RCONDE(:), RCONDV(:)
  double precision, allocatable :: Atmp(:,:)

  allocate( Atmp(n,n) )
  Atmp(1:n,1:n) = A(1:n,1:n)

  JOBVL  = "V" ! computes the left  eigenvectors 
  JOBVR  = "V" ! computes the right eigenvectors 
  BALANC = "B" ! Diagonal scaling and Permutation for optimization
  SENSE  = "B"
  lda  = n
  ldvl = n
  ldvr = n
  allocate(WORK(1),SCALE_array(n),RCONDE(n),RCONDV(n),IWORK(2*n-2))
  LWORK = -1 ! to ask for the optimal size of WORK
  call dgeevx(BALANC,JOBVL,JOBVR,SENSE,&  ! CHARACTERS 
              n,Atmp,lda,              &  ! MATRIX TO DIAGONALIZE
              WR,WI,                   &  ! REAL AND IMAGINARY PART OF EIGENVALUES 
              VL,ldvl,VR,ldvr,         &  ! LEFT AND RIGHT EIGENVECTORS 
              ILO,IHI,SCALE_array,ABNRM,RCONDE,RCONDV, & ! OUTPUTS OF OPTIMIZATION
              WORK,LWORK,IWORK,INFO)

  !if(INFO.gt.0)then
  ! print*,'dgeev failed !!',INFO
  if( INFO.ne.0 ) then
    print *, 'dgeevx failed !!', INFO
    stop
  endif

  LWORK = max(int(work(1)), 1) ! this is the optimal size of WORK 
  deallocate(WORK)
  allocate(WORK(LWORK))
  ! Actual dnon_hrmt_real_diag_newiagonalization 
  call dgeevx(BALANC,JOBVL,JOBVR,SENSE,&  ! CHARACTERS 
              n,Atmp,lda,              &  ! MATRIX TO DIAGONALIZE
              WR,WI,                   &  ! REAL AND IMAGINARY PART OF EIGENVALUES 
              VL,ldvl,VR,ldvr,         &  ! LEFT AND RIGHT EIGENVECTORS 
              ILO,IHI,SCALE_array,ABNRM,RCONDE,RCONDV, & ! OUTPUTS OF OPTIMIZATION
              WORK,LWORK,IWORK,INFO)

  !if(INFO.ne.0)then
  ! print*,'dgeev failed !!',INFO
  if( INFO.ne.0 ) then
    print *, 'dgeevx failed !!', INFO
    stop
  endif

  deallocate( Atmp )
  deallocate( WORK, SCALE_array, RCONDE, RCONDV, IWORK )

end subroutine lapack_diag_non_sym_new

! ---

subroutine lapack_diag_non_sym_right(n, A, WR, WI, VR)

  implicit none

  integer,          intent(in)  :: n
  double precision, intent(in)  :: A(n,n)
  double precision, intent(out) :: WR(n), WI(n), VR(n,n)

  integer                       :: i, lda, ldvl, ldvr, LWORK, INFO
  double precision, allocatable :: Atmp(:,:), WORK(:), VL(:,:)

  lda  = n
  ldvl = 1
  ldvr = n

  allocate( Atmp(n,n), VL(1,1) )
  Atmp(1:n,1:n) = A(1:n,1:n)

  allocate(WORK(1))
  LWORK = -1
  call dgeev('N', 'V', n, Atmp, lda, WR, WI, VL, ldvl, VR, ldvr, WORK, LWORK, INFO)
  if(INFO.gt.0)then
    print*,'dgeev failed !!',INFO
    stop
  endif

  LWORK = max(int(WORK(1)), 1) ! this is the optimal size of WORK 
  deallocate(WORK)

  allocate(WORK(LWORK))

  ! Actual diagonalization 
  call dgeev('N', 'V', n, Atmp, lda, WR, WI, VL, ldvl, VR, ldvr, WORK, LWORK, INFO)
  if(INFO.ne.0) then
    print*,'dgeev failed !!', INFO
    stop
  endif

  deallocate(Atmp, WORK, VL)

! print *, ' JOBL = F'
! print *, ' eigenvalues'
! do i = 1, n
!   write(*, '(1000(F16.10,X))') WR(i), WI(i)
! enddo
! print *, ' right eigenvect' 
! do i = 1, n
!   write(*, '(1000(F16.10,X))') VR(:,i)
! enddo

end subroutine lapack_diag_non_sym_right

! ---

subroutine non_hrmt_real_diag(n, A, leigvec, reigvec, n_real_eigv, eigval)

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

  integer                       :: i, j, n_good
  double precision              :: thr, threshold, accu_d, accu_nd
  integer,          allocatable :: list_good(:), iorder(:)
  double precision, allocatable :: Aw(:,:)
  double precision, allocatable :: WR(:), WI(:), Vl(:,:), VR(:,:), S(:,:), S_inv_half_tmp(:,:)

  print*, ' Computing the left/right eigenvectors with lapack ...'

  ! Eigvalue(n) = WR(n) + i * WI(n)
  allocate(WR(n), WI(n), VL(n,n), VR(n,n), Aw(n,n))
  Aw = A
  !print *, ' matrix to diagonalize', Aw
  call lapack_diag_non_sym(n, Aw, WR, WI, VL, VR)

  ! ---
  ! You track the real eigenvalues 

  thr = 1d-15

  n_good = 0
  do i = 1, n
    if(dabs(WI(i)).lt.thr) then
      n_good += 1
    else
      print*, ' Found an imaginary component to eigenvalue'
      print*, ' Re(i) + Im(i)', WR(i), WI(i)
    endif
  enddo

  allocate(list_good(n_good), iorder(n_good))
  n_good = 0
  do i = 1, n
    if(dabs(WI(i)).lt.thr) then
      n_good += 1
      list_good(n_good) = i
      eigval(n_good) = WR(i)
    endif
  enddo
  n_real_eigv = n_good
  do i = 1, n_good
   iorder(i) = i
  enddo

  ! You sort the real eigenvalues 
  call dsort(eigval, iorder, n_good)
  do i = 1, n_real_eigv
    do j = 1, n
      reigvec(j,i) = VR(j,list_good(iorder(i)))
      leigvec(j,i) = Vl(j,list_good(iorder(i)))
    enddo
  enddo

! print *, ' ordered eigenvalues'
! print *, ' right eigenvect' 
! do i = 1, n
!   print *, i, eigval(i)
!   write(*, '(1000(F16.10,X))') reigvec(:,i)
! enddo

  ! ---

  allocate( S(n_real_eigv,n_real_eigv), S_inv_half_tmp(n_real_eigv,n_real_eigv) )

  ! S = VL x VR
  call dgemm( 'T', 'N', n_real_eigv, n_real_eigv, n_real_eigv, 1.d0 &
            , leigvec, size(leigvec, 1), reigvec, size(reigvec, 1)  &
            , 0.d0, S, size(S, 1) )

  accu_nd = 0.d0
  accu_d  = 0.d0
  do i = 1, n_real_eigv
    do j = 1, n_real_eigv
      if(i==j) then
        accu_d += S(j,i)
      else
        accu_nd = accu_nd + S(j,i) * S(j,i)
      endif
    enddo
  enddo
  accu_nd = dsqrt(accu_nd)

  threshold = 1.d-15
  if( (accu_nd .gt. threshold) .or. (dabs(accu_d-dble(n_real_eigv)) .gt. threshold) ) then

    print*, ' sum of off-diag S elements = ', accu_nd
    print*, ' Should be zero '
    print*, ' sum of     diag S elements = ', accu_d
    print*, ' Should be ',n
    print*, ' Not bi-orthonormal !!'
    print*, ' Notice that if you are interested in ground state it is not a problem :)'
  endif

end subroutine non_hrmt_real_diag

! ---

subroutine lapack_diag_general_non_sym(n, A, B, WR, WI, VL, VR)

  BEGIN_DOC
  ! You enter with a general non hermitian matrix A(n,n) and another B(n,n)
  !
  ! You get out with the real WR and imaginary part WI of the eigenvalues 
  !
  ! Eigvalue(n) = (WR(n) + i * WI(n))
  !
  ! And the left VL and right VR eigenvectors 
  !
  ! VL(i,j) = <i|Psi_left(j)>  :: projection on the basis element |i> on the jth left  eigenvector 
  !
  ! VR(i,j) = <i|Psi_right(j)> :: projection on the basis element |i> on the jth right eigenvector 
  END_DOC

  implicit none
  integer,          intent(in)  :: n
  double precision, intent(in)  :: A(n,n), B(n,n)
  double precision, intent(out) :: WR(n), WI(n), VL(n,n), VR(n,n)

  integer                       :: lda, ldvl, ldvr, LWORK, INFO
  integer                       :: n_good
  double precision, allocatable :: WORK(:)
  double precision, allocatable :: Atmp(:,:)

  lda  = n
  ldvl = n
  ldvr = n

  allocate( Atmp(n,n) )
  Atmp(1:n,1:n) = A(1:n,1:n)

  allocate(WORK(1))
  LWORK = -1 
  call dgeev('V', 'V', n, Atmp, lda, WR, WI, VL, ldvl, VR, ldvr, WORK, LWORK, INFO)
  if(INFO.gt.0) then
    print*,'dgeev failed !!',INFO
    stop
  endif

  LWORK = max(int(WORK(1)), 1) 
  deallocate(WORK)

  allocate(WORK(LWORK))

  call dgeev('V', 'V', n, Atmp, lda, WR, WI, VL, ldvl, VR, ldvr, WORK, LWORK, INFO)
  if(INFO.ne.0) then
    print*,'dgeev failed !!', INFO
    stop
  endif

  deallocate( WORK, Atmp )

end subroutine lapack_diag_general_non_sym

! ---

subroutine non_hrmt_general_real_diag(n, A, B, reigvec, leigvec, n_real_eigv, eigval)

  BEGIN_DOC
  ! routine which returns the sorted REAL EIGENVALUES ONLY and corresponding LEFT/RIGHT eigenvetors 
  !
  ! of a non hermitian matrix A(n,n) and B(n,n) 
  !
  ! A reigvec = eigval * B * reigvec
  !
  ! (A)^\dagger leigvec = eigval * B * leigvec
  !
  ! n_real_eigv is the number of real eigenvalues, which might be smaller than the dimension "n" 
  END_DOC

  implicit none
  integer,          intent(in)  :: n
  double precision, intent(in)  :: A(n,n), B(n,n)
  integer,          intent(out) :: n_real_eigv
  double precision, intent(out) :: reigvec(n,n), leigvec(n,n), eigval(n)

  integer                       :: i, j
  integer                       :: n_good
  integer, allocatable          :: list_good(:), iorder(:)
  double precision, allocatable :: WR(:), WI(:), Vl(:,:), VR(:,:)
  double precision, allocatable :: Aw(:,:), Bw(:,:)

  print*,'Computing the left/right eigenvectors ...'

  allocate(WR(n), WI(n), VL(n,n), VR(n,n), Aw(n,n), Bw(n,n))
  Aw = A
  Bw = B

  call lapack_diag_general_non_sym(n, A, B, WR, WI, VL, VR)

  ! You track the real eigenvalues 
  n_good = 0
  do i = 1, n
    if(dabs(WI(i)) .lt. 1.d-12) then
      n_good += 1
    else
      print*,'Found an imaginary component to eigenvalue'
      print*,'Re(i) + Im(i)',WR(i),WI(i)
    endif
  enddo

  allocate(list_good(n_good), iorder(n_good))
  n_good = 0
  do i = 1, n
    if(dabs(WI(i)).lt.1.d-12)then
      n_good += 1
      list_good(n_good) = i
      eigval(n_good) = WR(i)
    endif
  enddo
  n_real_eigv = n_good 
  do i = 1, n_good
   iorder(i) = i
  enddo

  ! You sort the real eigenvalues 
  call dsort(eigval, iorder, n_good)
  print*,'n_real_eigv = ', n_real_eigv
  print*,'n           = ', n
  do i = 1, n_real_eigv
    print*,i,'eigval(i) = ', eigval(i) 
    do j = 1, n
      reigvec(j,i) = VR(j,list_good(iorder(i)))
      leigvec(j,i) = Vl(j,list_good(iorder(i)))
    enddo
  enddo

end subroutine non_hrmt_general_real_diag

! ---

subroutine impose_biorthog_qr(m, n, thr_d, thr_nd, Vl, Vr)

  implicit none 
  integer,          intent(in)    :: m, n
  double precision, intent(in)    :: thr_d, thr_nd
  double precision, intent(inout) :: Vl(m,n), Vr(m,n)

  integer                         :: i, j
  integer                         :: LWORK, INFO
  double precision                :: accu_nd, accu_d
  double precision, allocatable   :: TAU(:), WORK(:)
  double precision, allocatable   :: S(:,:), R(:,:), tmp(:,:)

  ! ---

  call check_biorthog_binormalize(m, n, Vl, Vr, thr_d, thr_nd, .false.)
  
  ! ---
  
  allocate(S(n,n))
  call dgemm( 'T', 'N', n, n, m, 1.d0          &
            , Vl, size(Vl, 1), Vr, size(Vr, 1) &
            , 0.d0, S, size(S, 1) )

  accu_nd = 0.d0
  accu_d  = 0.d0
  do i = 1, n
    do j = 1, n
      if(i==j) then
        accu_d += S(j,i)
      else
        accu_nd = accu_nd + S(j,i) * S(j,i)
      endif
    enddo
  enddo
  accu_nd = dsqrt(accu_nd)

  if((accu_nd .lt. thr_nd) .and. (dabs(accu_d-dble(n))/dble(n) .lt. thr_d)) then
    print *, ' bi-orthogonal vectors without QR !'
    deallocate(S)
    return
  endif

  ! -------------------------------------------------------------------------------------
  !                           QR factorization of S: S = Q x R


  print *, ' apply QR decomposition ...'

  allocate( TAU(n), WORK(1) )

  LWORK = -1
  call dgeqrf(n, n, S, n, TAU, WORK, LWORK, INFO)
  if(INFO .ne. 0) then
    print*,'dgeqrf failed !!', INFO
    stop
  endif

  LWORK = max(n, int(WORK(1)))
  deallocate(WORK)

  allocate( WORK(LWORK) )
  call dgeqrf(n, n, S, n, TAU, WORK, LWORK, INFO)
  if(INFO .ne. 0) then
    print*,'dgeqrf failed !!', INFO
    stop
  endif

  ! save the upper triangular R
  allocate( R(n,n) )
  R(:,:) = S(:,:)

  ! get Q
  LWORK = -1
  call dorgqr(n, n, n, S, n, TAU, WORK, LWORK, INFO)
  if(INFO .ne. 0) then
    print*,'dorgqr failed !!', INFO
    stop
  endif

  LWORK = max(n, int(WORK(1)))
  deallocate(WORK)

  allocate( WORK(LWORK) )
  call dorgqr(n, n, n, S, n, TAU, WORK, LWORK, INFO)
  if(INFO .ne. 0) then
    print*,'dorgqr failed !!', INFO
    stop
  endif

  deallocate( WORK, TAU )

  !
  ! -------------------------------------------------------------------------------------

  ! ---

  ! -------------------------------------------------------------------------------------
  !                               get bi-orhtog left & right vectors:
  !                                           Vr' = Vr x inv(R) 
  !                                           Vl' = inv(Q) x Vl =  Q.T   x Vl 

  ! Q.T x Vl, where Q = S

  allocate( tmp(n,m) )
  call dgemm( 'T', 'T', n, m, n, 1.d0        &
            , S, size(S, 1), Vl, size(Vl, 1) &
            , 0.d0, tmp, size(tmp, 1) )

  do i = 1, n
    do j = 1, m
      Vl(j,i) = tmp(i,j)
    enddo
  enddo
  deallocate(tmp)

  ! ---

  ! inv(R) 
  !print *, ' inversing upper triangular matrix ...'
  call dtrtri("U", "N", n, R, n, INFO)
  if(INFO .ne. 0) then
    print*,'dtrtri failed !!', INFO
    stop
  endif
  !print *, ' inversing upper triangular matrix OK' 

  do i = 1, n-1
    do j = i+1, n
      R(j,i) = 0.d0
    enddo
  enddo

  !print *, ' inv(R):'
  !do i = 1, n
  !  write(*, '(1000(F16.10,X))') R(i,:)
  !enddo

  ! Vr x inv(R) 
  allocate( tmp(m,n) )
  call dgemm( 'N', 'N', m, n, n, 1.d0        &
            , Vr, size(Vr, 1), R, size(R, 1) &
            , 0.d0, tmp, size(tmp, 1) )
  deallocate( R )

  do i = 1, n
    do j = 1, m
      Vr(j,i) = tmp(j,i)
    enddo
  enddo
  deallocate(tmp)

  return
end subroutine impose_biorthog_qr

! ---

subroutine impose_biorthog_lu(m, n, Vl, Vr, S)

  implicit none 
  integer, intent(in)             :: m, n
  double precision, intent(inout) :: Vl(m,n), Vr(m,n), S(n,n)

  integer                         :: i, j
  integer                         :: INFO
  double precision                :: nrm
  integer,          allocatable   :: IPIV(:)
  double precision, allocatable   :: L(:,:), tmp(:,:), vectmp(:)
  !double precision, allocatable   :: T(:,:), ll(:,:), rr(:,:), tt(:,:)

  !allocate( T(n,n) )
  !T(:,:) = S(:,:)

  print *, ' apply LU decomposition ...'

  ! -------------------------------------------------------------------------------------
  !                           LU factorization of S: S = P x L x U

  allocate( IPIV(n) )

  call dgetrf(n, n, S, n, IPIV, INFO)
  if(INFO .ne. 0) then
    print*, 'dgetrf failed !!', INFO
    stop
  endif

  ! check | S - P x L x U |
  !allocate( ll(n,n), rr(n,n), tmp(n,n) )
  !ll = S
  !rr = S
  !do i = 1, n-1
  !  ll(i,i) = 1.d0
  !  do j = i+1, n
  !    ll(i,j) = 0.d0
  !    rr(j,i) = 0.d0
  !  enddo
  !enddo
  !ll(n,n) = 1.d0
  !call dgemm( 'N', 'N', n, n, n, 1.d0          &
  !          , ll, size(ll, 1), rr, size(rr, 1) &
  !          , 0.d0, tmp, size(tmp, 1) )
  ! deallocate(ll, rr)
  !allocate( vectmp(n) )
  !do j = n-1, 1, -1
  !  i = IPIV(j)
  !  if(i.ne.j) then
  !    print *, j, i
  !    vectmp(:) = tmp(i,:)
  !    tmp(i,:)  = tmp(j,:)
  !    tmp(j,:)  = vectmp(:)
  !  endif
  !enddo
  !deallocate( vectmp )
  !nrm = 0.d0
  !do i = 1, n
  !  do j = 1, n
  !    nrm += dabs(tmp(j,i) - T(j,i))
  !  enddo
  !enddo
  !deallocate( tmp )
  !print*, '|L.T x R - S| =', nrm
  !stop

  ! ------
  ! inv(L) 
  ! ------

  allocate( L(n,n) )
  L(:,:) = S(:,:)

  call dtrtri("L", "U", n, L, n, INFO)
  if(INFO .ne. 0) then
    print*,  'dtrtri failed !!', INFO
    stop
  endif
  do i = 1, n-1
    L(i,i) = 1.d0
    do j = i+1, n
      L(i,j) = 0.d0
    enddo
  enddo
  L(n,n) = 1.d0

  ! ------
  ! inv(U) 
  ! ------
  
  call dtrtri("U", "N", n, S, n, INFO)
  if(INFO .ne. 0) then
    print*,  'dtrtri failed !!', INFO
    stop
  endif

  do i = 1, n-1
    do j = i+1, n
      S(j,i) = 0.d0
    enddo
  enddo

  !
  ! -------------------------------------------------------------------------------------

  ! ---

  ! -------------------------------------------------------------------------------------
  !                               get bi-orhtog left & right vectors:
  !                                           Vr' = Vr x inv(U) 
  !                                           Vl' = inv(L) x inv(P) x Vl

  ! inv(P) x Vl
  allocate( vectmp(n) )
  do j = n-1, 1, -1
    i = IPIV(j)
    if(i.ne.j) then
      vectmp(:) = L(:,j)
      L(:,j)    = L(:,i)
      L(:,i)    = vectmp(:)
    endif
  enddo
  deallocate( vectmp )

  ! Vl'
  allocate( tmp(m,n) )
  call dgemm( 'N', 'T', m, n, n, 1.d0        &
            , Vl, size(Vl, 1), L, size(L, 1) &
            , 0.d0, tmp, size(tmp, 1) )
  deallocate(L)

  Vl = tmp
  deallocate(tmp)

  ! ---

  ! Vr x inv(U) 
  allocate( tmp(m,n) )
  call dgemm( 'N', 'N', m, n, n, 1.d0        &
            , Vr, size(Vr, 1), S, size(S, 1) &
            , 0.d0, tmp, size(tmp, 1) )
  Vr = tmp
  deallocate(tmp)

  !allocate( tmp(n,n) )
  !call dgemm( 'T', 'N', n, n, m, 1.d0          &
  !          , Vl, size(Vl, 1), Vr, size(Vr, 1) &
  !          , 0.d0, tmp, size(tmp, 1) )
  !nrm = 0.d0
  !do i = 1, n
  !  do j = 1, n
  !    nrm += dabs(tmp(j,i))
  !  enddo
  !enddo
  !deallocate( tmp )
  !print*, '|L.T x R| =', nrm
  !stop

  return
end subroutine impose_biorthog_lu

! ---

subroutine check_EIGVEC(n, m, A, eigval, leigvec, reigvec, thr_diag, thr_norm, stop_ifnot)

  implicit none
  integer,          intent(in)  :: n, m
  logical,          intent(in)  :: stop_ifnot
  double precision, intent(in)  :: A(n,n), eigval(m), leigvec(n,m), reigvec(n,m), thr_diag, thr_norm
 
  integer                       :: i, j
  double precision              :: tmp, tmp_abs, tmp_nrm, tmp_rel, tmp_dif
  double precision              :: V_nrm, U_nrm
  double precision, allocatable :: Mtmp(:,:)

  allocate( Mtmp(n,m) )
  
  ! ---

  Mtmp = 0.d0
  call dgemm( 'N', 'N', n, m, n, 1.d0                  &
            , A, size(A, 1), reigvec, size(reigvec, 1) &
            , 0.d0, Mtmp, size(Mtmp, 1) )

  V_nrm   = 0.d0
  tmp_nrm = 0.d0
  tmp_abs = 0.d0
  do j = 1, m
    
    tmp   = 0.d0
    U_nrm = 0.d0
    do i = 1, n
      tmp     = tmp     + dabs(Mtmp(i,j) - eigval(j) * reigvec(i,j))
      tmp_nrm = tmp_nrm + dabs(Mtmp(i,j))
      U_nrm   = U_nrm   + reigvec(i,j) * reigvec(i,j)
    enddo

    tmp_abs = tmp_abs + tmp
    V_nrm   = V_nrm   + U_nrm 
    !write(*,'(I4,X,(100(F25.16,X)))') j,eigval(j), tmp, U_nrm

  enddo

  if(tmp_abs.lt.10.d-10)then
   tmp_rel = thr_diag/10.d0
  else
   tmp_rel = tmp_abs / tmp_nrm
  endif
  tmp_dif = dabs(V_nrm - dble(m))

  if( stop_ifnot .and. ((tmp_rel .gt. thr_diag) .or. (tmp_dif .gt. thr_norm)) ) then
    print *, ' error in right-eigenvectors'
    print *, ' err tol   = ',thr_diag, thr_norm
    print *, '(tmp_rel .gt. thr_diag) = ',(tmp_rel .gt. thr_diag)
    print *, '(tmp_dif .gt. thr_norm) = ',(tmp_dif .gt. thr_norm)
    print *, ' err estim = ', tmp_abs, tmp_rel
    print *, ' CR norm   = ', V_nrm 
    stop
  endif

  ! ---

  Mtmp = 0.d0
  call dgemm( 'T', 'N', n, m, n, 1.d0                  &
            , A, size(A, 1), leigvec, size(leigvec, 1) &
            , 0.d0, Mtmp, size(Mtmp, 1) )

  V_nrm   = 0.d0
  tmp_nrm = 0.d0
  tmp_abs = 0.d0
  do j = 1, m

    tmp   = 0.d0
    U_nrm = 0.d0
    do i = 1, n
      tmp     = tmp     + dabs(Mtmp(i,j) - eigval(j) * leigvec(i,j))
      tmp_nrm = tmp_nrm + dabs(Mtmp(i,j))
      U_nrm   = U_nrm   + leigvec(i,j) * leigvec(i,j)
    enddo

    tmp_abs = tmp_abs + tmp
    V_nrm   = V_nrm   + U_nrm 
    !write(*,'(I4,X,(100(F25.16,X)))') j,eigval(j), tmp, U_nrm

  enddo

  if(tmp_abs.lt.10.d-10)then
   tmp_rel = thr_diag/10.d0
  else
   tmp_rel = tmp_abs / tmp_nrm
  endif
  if( stop_ifnot .and. ((tmp_rel .gt. thr_diag) .or. (tmp_dif .gt. thr_norm)) ) then
    print *, ' error in left-eigenvectors'
    print *, ' err tol   = ',thr_diag, thr_norm
    print *, '(tmp_rel .gt. thr_diag) = ',(tmp_rel .gt. thr_diag)
    print *, '(tmp_dif .gt. thr_norm) = ',(tmp_dif .gt. thr_norm)
    print *, ' err estim = ', tmp_abs, tmp_rel
    print *, ' CR norm   = ', V_nrm 
    stop
  endif

  ! ---

  deallocate( Mtmp )

end subroutine check_EIGVEC

! ---

subroutine check_degen(n, m, eigval, leigvec, reigvec)

  implicit none
  integer,          intent(in)    :: n, m
  double precision, intent(in)    :: eigval(m)
  double precision, intent(inout) :: leigvec(n,m), reigvec(n,m)
 
  integer                         :: i, j
  double precision                :: ei, ej, de, de_thr, accu_nd
  double precision, allocatable   :: S(:,:)

  de_thr = 1d-6

  do i = 1, m-1
    ei = eigval(i)

    do j = i+1, m
      ej = eigval(j)
      de = dabs(ei - ej)

      if(de .lt. de_thr) then

        leigvec(:,i) = 0.d0
        leigvec(:,j) = 0.d0
        leigvec(i,i) = 1.d0
        leigvec(j,j) = 1.d0

        reigvec(:,i) = 0.d0
        reigvec(:,j) = 0.d0
        reigvec(i,i) = 1.d0
        reigvec(j,j) = 1.d0

      endif

    enddo
  enddo

  ! ---

  allocate( S(m,m) )

  ! S = VL x VR
  call dgemm( 'T', 'N', m, m, n, 1.d0                              &
            , leigvec, size(leigvec, 1), reigvec, size(reigvec, 1) &
            , 0.d0, S, size(S, 1) )

  accu_nd = 0.d0
  do i = 1, m
    do j = 1, m
      if(i==j) cycle
      accu_nd = accu_nd + S(j,i) * S(j,i)
    enddo
  enddo
  accu_nd = dsqrt(accu_nd)

  deallocate( S )

  print *, ' check_degen: L & T bi-orthogonality: ok'
  print *, ' accu_nd = ', accu_nd

  if( accu_nd .lt. 1d-8 ) then
    return
  else
    stop
  endif

end subroutine check_degen

! ---

subroutine impose_weighted_orthog_svd(n, m, W, C)

  implicit none

  integer,          intent(in)    :: n, m
  double precision, intent(inout) :: C(n,m), W(n,n)

  integer                         :: i, j, num_linear_dependencies
  double precision                :: threshold
  double precision, allocatable   :: S(:,:), tmp(:,:)
  double precision, allocatable   :: U(:,:), Vt(:,:), D(:)

  !print *, ' apply SVD to orthogonalize & normalize weighted vectors'

  ! ---

  ! C.T x W x C
  allocate(S(m,m))
  allocate(tmp(m,n))
  call dgemm( 'T', 'N', m, n, n, 1.d0      &
            , C, size(C, 1), W, size(W, 1) &
            , 0.d0, tmp, size(tmp, 1) )
  call dgemm( 'N', 'N', m, m, n, 1.d0          &
            , tmp, size(tmp, 1), C, size(C, 1) &
            , 0.d0, S, size(S, 1) )
  deallocate(tmp)

  !print *, ' overlap bef SVD: '
  !do i = 1, m
  !  write(*, '(1000(F16.10,X))') S(i,:)
  !enddo

  ! ---
 
  allocate(U(m,m), Vt(m,m), D(m))

  call svd(S, m, U, m, D, Vt, m, m, m)

  deallocate(S)

  threshold               = 1.d-6
  num_linear_dependencies = 0
  do i = 1, m
    if(abs(D(i)) <= threshold) then
      D(i) = 0.d0
      num_linear_dependencies += 1
    else
      ASSERT (D(i) > 0.d0)
      D(i) = 1.d0 / dsqrt(D(i))
    endif
  enddo
  if(num_linear_dependencies > 0) then
    write(*,*) ' linear dependencies = ', num_linear_dependencies
    write(*,*) ' m                   = ', m
    stop
  endif

  ! ---

  allocate(tmp(n,m))

  ! tmp <-- C x U
  call dgemm( 'N', 'N', n, m, m, 1.d0      &
            , C, size(C, 1), U, size(U, 1) &
            , 0.d0, tmp, size(tmp, 1) )

  deallocate(U, Vt)

  ! C <-- tmp x sigma^-0.5
  do j = 1, m
    do i = 1, n
      C(i,j) = tmp(i,j) * D(j)
    enddo
  enddo

  deallocate(D, tmp)

  ! ---

  ! C.T x W x C
  allocate(S(m,m))
  allocate(tmp(m,n))
  call dgemm( 'T', 'N', m, n, n, 1.d0      &
            , C, size(C, 1), W, size(W, 1) &
            , 0.d0, tmp, size(tmp, 1) )
  call dgemm( 'N', 'N', m, m, n, 1.d0          &
            , tmp, size(tmp, 1), C, size(C, 1) &
            , 0.d0, S, size(S, 1) )
  deallocate(tmp)

  !print *, ' overlap aft SVD: '
  !do i = 1, m
  !  write(*, '(1000(F16.10,X))') S(i,:)
  !enddo

  deallocate(S)

  ! ---

end subroutine impose_weighted_orthog_svd

! ---

subroutine impose_orthog_svd(n, m, C)

  implicit none

  integer,          intent(in)    :: n, m
  double precision, intent(inout) :: C(n,m)

  integer                         :: i, j, num_linear_dependencies
  double precision                :: threshold
  double precision, allocatable   :: S(:,:), tmp(:,:)
  double precision, allocatable   :: U(:,:), Vt(:,:), D(:)

  !print *, ' apply SVD to orthogonalize & normalize vectors'

  ! ---

  allocate(S(m,m))

  ! S = C.T x C
  call dgemm( 'T', 'N', m, m, n, 1.d0      &
            , C, size(C, 1), C, size(C, 1) &
            , 0.d0, S, size(S, 1) )

  !print *, ' eigenvec overlap bef SVD: '
  !do i = 1, m
  !  write(*, '(1000(F16.10,X))') S(i,:)
  !enddo

  ! ---
 
  allocate(U(m,m), Vt(m,m), D(m))

  call svd(S, m, U, m, D, Vt, m, m, m)

  deallocate(S)

  threshold               = 1.d-6
  num_linear_dependencies = 0
  do i = 1, m
    if(abs(D(i)) <= threshold) then
      write(*,*) ' D(i) = ', D(i)
      D(i) = 0.d0
      num_linear_dependencies += 1
    else
      ASSERT (D(i) > 0.d0)
      D(i) = 1.d0 / dsqrt(D(i))
    endif
  enddo
  if(num_linear_dependencies > 0) then
    write(*,*) ' linear dependencies = ', num_linear_dependencies
    write(*,*) ' m                   = ', m
    write(*,*) ' try with Graham-Schmidt'
    stop
  endif

  ! ---

  allocate(tmp(n,m))

  ! tmp <-- C x U
  call dgemm( 'N', 'N', n, m, m, 1.d0      &
            , C, size(C, 1), U, size(U, 1) &
            , 0.d0, tmp, size(tmp, 1) )

  deallocate(U, Vt)

  ! C <-- tmp x sigma^-0.5
  do j = 1, m
    do i = 1, n
      C(i,j) = tmp(i,j) * D(j)
    enddo
  enddo

  deallocate(D, tmp)

  ! ---

  allocate(S(m,m))

  ! S = C.T x C
  call dgemm( 'T', 'N', m, m, n, 1.d0      &
            , C, size(C, 1), C, size(C, 1) &
            , 0.d0, S, size(S, 1) )

  !print *, ' eigenvec overlap aft SVD: '
  !do i = 1, m
  !  write(*, '(1000(F16.10,X))') S(i,:)
  !enddo

  deallocate(S)

  ! ---

end subroutine impose_orthog_svd

! ---

subroutine impose_orthog_svd_overlap(n, m, C, overlap)

  implicit none

  integer,          intent(in)    :: n, m
  double precision, intent(in   ) :: overlap(n,n)
  double precision, intent(inout) :: C(n,m)

  integer                         :: i, j, num_linear_dependencies
  double precision                :: threshold
  double precision, allocatable   :: S(:,:), tmp(:,:), Stmp(:,:)
  double precision, allocatable   :: U(:,:), Vt(:,:), D(:)

  print *, ' apply SVD to orthogonalize vectors'

  ! ---

  ! S = C.T x overlap x C
  allocate(S(m,m), Stmp(n,m))
  call dgemm( 'N', 'N', n, m, n, 1.d0                  &
            , overlap, size(overlap, 1), C, size(C, 1) &
            , 0.d0, Stmp, size(Stmp, 1) )
  call dgemm( 'T', 'N', m, m, n, 1.d0            &
            , C, size(C, 1), Stmp, size(Stmp, 1) &
            , 0.d0, S, size(S, 1) )
  deallocate(Stmp)

  !print *, ' eigenvec overlap bef SVD: '
  !do i = 1, m
  !  write(*, '(1000(F16.10,X))') S(i,:)
  !enddo

  ! ---
 
  allocate(U(m,m), Vt(m,m), D(m))

  call svd(S, m, U, m, D, Vt, m, m, m)

  deallocate(S)

  threshold               = 1.d-6
  num_linear_dependencies = 0
  do i = 1, m
    if(abs(D(i)) <= threshold) then
      D(i) = 0.d0
      num_linear_dependencies += 1
    else
      ASSERT (D(i) > 0.d0)
      D(i) = 1.d0 / dsqrt(D(i))
    endif
  enddo
  if(num_linear_dependencies > 0) then
    write(*,*) ' linear dependencies = ', num_linear_dependencies
    write(*,*) ' m                   = ', m
    stop
  endif

  ! ---

  allocate(tmp(n,m))

  ! tmp <-- C x U
  call dgemm( 'N', 'N', n, m, m, 1.d0      &
            , C, size(C, 1), U, size(U, 1) &
            , 0.d0, tmp, size(tmp, 1) )

  deallocate(U, Vt)

  ! C <-- tmp x sigma^-0.5
  do j = 1, m
    do i = 1, n
      C(i,j) = tmp(i,j) * D(j)
    enddo
  enddo

  deallocate(D, tmp)

  ! ---

  ! S = C.T x overlap x C
  allocate(S(m,m), Stmp(n,m))
  call dgemm( 'N', 'N', n, m, n, 1.d0                  &
            , overlap, size(overlap, 1), C, size(C, 1) &
            , 0.d0, Stmp, size(Stmp, 1) )
  call dgemm( 'T', 'N', m, m, n, 1.d0            &
            , C, size(C, 1), Stmp, size(Stmp, 1) &
            , 0.d0, S, size(S, 1) )
  deallocate(Stmp)

  !print *, ' eigenvec overlap aft SVD: '
  !do i = 1, m
  !  write(*, '(1000(F16.10,X))') S(i,:)
  !enddo
  deallocate(S)

end subroutine impose_orthog_svd_overlap

! ---

subroutine impose_orthog_GramSchmidt(n, m, C)

  implicit none

  integer,          intent(in)    :: n, m
  double precision, intent(inout) :: C(n,m)

  integer                         :: i, j, k
  double precision                :: Ojk, Ojj, fact_ct
  double precision, allocatable   :: S(:,:)

  print *, ''
  print *, ' apply Gram-Schmidt to orthogonalize & normalize vectors'
  print *, ''

  ! ---

  allocate(S(m,m))
  call dgemm( 'T', 'N', m, m, n, 1.d0      &
            , C, size(C, 1), C, size(C, 1) &
            , 0.d0, S, size(S, 1) )

  print *, ' eigenvec overlap bef Gram-Schmidt: '
  do i = 1, m
    write(*, '(1000(F16.10,X))') S(i,:)
  enddo

  ! ---

  do k = 2, m
    do j = 1, k-1

      Ojk = 0.d0    
      Ojj = 0.d0    
      do i = 1, n
        Ojk = Ojk + C(i,j) * C(i,k)
        Ojj = Ojj + C(i,j) * C(i,j)
      enddo
      fact_ct = Ojk / Ojj

      do i = 1, n
        C(i,k) = C(i,k) - fact_ct * C(i,j)
      enddo

    enddo
  enddo

  do k = 1, m
    fact_ct = 0.d0    
    do i = 1, n
      fact_ct = fact_ct + C(i,k) * C(i,k)
    enddo
    fact_ct = dsqrt(fact_ct)
    do i = 1, n
      C(i,k) = C(i,k) / fact_ct
    enddo
  enddo

  ! ---

  call dgemm( 'T', 'N', m, m, n, 1.d0      &
            , C, size(C, 1), C, size(C, 1) &
            , 0.d0, S, size(S, 1) )

  print *, ' eigenvec overlap aft Gram-Schmidt: '
  do i = 1, m
    write(*, '(1000(F16.10,X))') S(i,:)
  enddo

  deallocate(S)

  ! ---

end subroutine impose_orthog_GramSchmidt

! ---

subroutine impose_orthog_ones(n, deg_num, C)


  implicit none

  integer,          intent(in)    :: n
  integer,          intent(in)    :: deg_num(n)
  double precision, intent(inout) :: C(n,n)

  integer                         :: i, j, ii, di, dj

  print *, ''
  print *, ' orthogonalize vectors by hand'
  print *, ''

  do i = 1, n-1
    di = deg_num(i)

    if(di .gt. 1) then

      do ii = 1, di
        C(:     ,i+ii-1) = 0.d0
        C(i+ii-1,i+ii-1) = 1.d0
      enddo

      do j = i+di+1, n
        dj = deg_num(j)
        if(dj .eq. di) then
          do ii = 1, dj
            C(:,     j+ii-1) = 0.d0
            C(j+ii-1,j+ii-1) = 1.d0
          enddo
        endif
      enddo

    endif
  enddo 

end subroutine impose_orthog_ones

! ---

subroutine impose_orthog_degen_eigvec(n, e0, C0)

  implicit none

  integer,          intent(in)    :: n
  double precision, intent(in)    :: e0(n)
  double precision, intent(inout) :: C0(n,n)

  integer                         :: i, j, k, m
  double precision                :: ei, ej, de, de_thr
  integer,          allocatable   :: deg_num(:)
  double precision, allocatable   :: C(:,:)

  ! ---

  allocate( deg_num(n) )
  do i = 1, n
    deg_num(i) = 1
  enddo

  de_thr = thr_degen_tc

  do i = 1, n-1
    ei = e0(i)

    ! already considered in degen vectors
    if(deg_num(i).eq.0) cycle

    do j = i+1, n
      ej = e0(j)
      de = dabs(ei - ej)

      if(de .lt. de_thr) then
        deg_num(i) = deg_num(i) + 1 
        deg_num(j) = 0
      endif

    enddo
  enddo

  
  !do i = 1, n
  !  if(deg_num(i) .gt. 1) then
  !    print *, ' degen on', i, deg_num(i)
  !  endif
  !enddo

  ! ---

!  call impose_orthog_ones(n, deg_num, C0)

  do i = 1, n
    m = deg_num(i)

    if(m .gt. 1) then
    !if(m.eq.3) then
  
      allocate(C(n,m))
      do j = 1, m
        C(1:n,j) = C0(1:n,i+j-1)
      enddo

      ! ---

      ! C <= C U sigma^-0.5
      call impose_orthog_svd(n, m, C)

      ! ---

      ! C = I
      !C = 0.d0
      !do j = 1, m
      !  C(i+j-1,j) = 1.d0
      !enddo

      ! ---

!      call impose_orthog_GramSchmidt(n, m, C)

      ! ---

      do j = 1, m
        C0(1:n,i+j-1) = C(1:n,j)
      enddo
      deallocate(C)

    endif
  enddo

end subroutine impose_orthog_degen_eigvec 

! ---

subroutine get_halfinv_svd(n, S)

  implicit none

  integer,          intent(in)    :: n
  double precision, intent(inout) :: S(n,n)

  integer                         :: num_linear_dependencies
  integer                         :: i, j, k
  double precision                :: accu_d, accu_nd, thresh
  double precision, parameter     :: threshold = 1.d-6
  double precision, allocatable   :: U(:,:), Vt(:,:), D(:)
  double precision, allocatable   :: S0(:,:), Stmp(:,:), Stmp2(:,:)

  allocate( S0(n,n) )
  S0(1:n,1:n) = S(1:n,1:n)

  allocate(U(n,n), Vt(n,n), D(n))
  call svd(S, n, U, n, D, Vt, n, n, n)

  num_linear_dependencies = 0
  do i = 1, n
    if(abs(D(i)) <= threshold) then
      D(i) = 0.d0
      num_linear_dependencies += 1
    else
      ASSERT (D(i) > 0.d0)
      D(i) = 1.d0 / dsqrt(D(i))
    endif
  enddo
  write(*,*) ' linear dependencies', num_linear_dependencies

  S(:,:) = 0.d0
  do k = 1, n
    if(D(k) /= 0.d0) then
      do j = 1, n
        do i = 1, n
          S(i,j) = S(i,j) + U(i,k) * D(k) * Vt(k,j)
        enddo
      enddo
    endif
  enddo
  deallocate(U, D, Vt)

  allocate( Stmp(n,n), Stmp2(n,n) )
  Stmp  = 0.d0
  Stmp2 = 0.d0
  ! S^-1/2 x S
  call dgemm( 'N', 'N', n, n, n, 1.d0        &
            , S, size(S, 1), S0, size(S0, 1) &
            , 0.d0, Stmp, size(Stmp, 1) )
  ! ( S^-1/2 x S ) x S^-1/2
  call dgemm( 'N', 'N', n, n, n, 1.d0            &
            , Stmp, size(Stmp, 1), S, size(S, 1) &
            , 0.d0, Stmp2, size(Stmp2, 1) )

  accu_nd = 0.d0
  accu_d  = 0.d0
  thresh  = 1.d-10
  do i = 1, n
    do j = 1, n
      if(i==j) then
       accu_d += Stmp2(j,i)
      else 
       accu_nd = accu_nd + Stmp2(j,i) * Stmp2(j,i)
      endif
    enddo
  enddo
  accu_nd = dsqrt(accu_nd)
  if( accu_nd.gt.thresh .or. dabs(accu_d-dble(n)).gt.thresh) then
    print*, ' after S^-1/2: sum of off-diag S elements = ', accu_nd
    print*, ' after S^-1/2: sum of     diag S elements = ', accu_d
    do i = 1, n
      write(*,'(1000(F16.10,X))') Stmp2(i,:)
    enddo
    stop
  endif

  deallocate(S0, Stmp, Stmp2)

end subroutine get_halfinv_svd

! ---

subroutine check_biorthog_binormalize(n, m, Vl, Vr, thr_d, thr_nd, stop_ifnot)

  implicit none
  
  integer,          intent(in)    :: n, m
  logical,          intent(in)    :: stop_ifnot
  double precision, intent(in)    :: thr_d, thr_nd
  double precision, intent(inout) :: Vl(n,m), Vr(n,m)

  integer                         :: i, j
  double precision                :: accu_d, accu_nd, s_tmp
  double precision, allocatable   :: S(:,:)

  !print *, ' check bi-orthonormality'

  ! ---

  allocate(S(m,m))
  call dgemm( 'T', 'N', m, m, n, 1.d0          &
            , Vl, size(Vl, 1), Vr, size(Vr, 1) &
            , 0.d0, S, size(S, 1) )
  !print *, ' overlap matrix before:'
  !do i = 1, m
  !  write(*,'(1000(F16.10,X))') S(i,:)
  !enddo

  ! S(i,i) = -1
  do i = 1, m
    if(S(i,i) .lt. 0.d0) then
    !if( (S(i,i) + 1.d0) .lt. thr_d ) then
      do j = 1, n
        Vl(j,i) = -1.d0 * Vl(j,i)
      enddo
      !S(i,i) = 1.d0
      S(i,i) = -S(i,i)
    endif
  enddo

  accu_d  = 0.d0
  accu_nd = 0.d0
  do i = 1, m
    do j = 1, m
      if(i==j) then
        accu_d = accu_d + S(i,i)
      else
        accu_nd = accu_nd + S(j,i) * S(j,i)
      endif
    enddo
  enddo
  accu_nd = dsqrt(accu_nd) / dble(m)
  !print*, '    diag acc bef = ', accu_d
  !print*, ' nondiag acc bef = ', accu_nd

  ! ---

  if( (accu_nd .lt. thr_nd) .and. (dabs(accu_d-dble(m))/dble(m) .gt. thr_d) ) then

    do i = 1, m
      if(S(i,i) <= 0.d0) then
        print *, ' overap negative'
        print *, i, S(i,i)
        exit
      endif
      if(dabs(S(i,i) - 1.d0) .gt. thr_d) then
        s_tmp = 1.d0 / dsqrt(S(i,i))
        do j = 1, n
          Vl(j,i) = Vl(j,i) * s_tmp 
          Vr(j,i) = Vr(j,i) * s_tmp 
        enddo
      endif

    enddo

  endif

  ! ---

  call dgemm( 'T', 'N', m, m, n, 1.d0          &
            , Vl, size(Vl, 1), Vr, size(Vr, 1) &
            , 0.d0, S, size(S, 1) )
  !print *, ' overlap matrix after:'
  !do i = 1, m
  !  write(*,'(1000(F16.10,X))') S(i,:)
  !enddo

  accu_d  = 0.d0
  accu_nd = 0.d0
  do i = 1, m
    do j = 1, m
      if(i==j) then
        accu_d = accu_d + S(i,i)
      else
        accu_nd = accu_nd + S(j,i) * S(j,i)
      endif
    enddo
  enddo
  accu_nd = dsqrt(accu_nd) / dble(m)
  !print *, '    diag acc aft = ', accu_d
  !print *, ' nondiag acc aft = ', accu_nd

  deallocate(S)

  ! ---

  if( stop_ifnot .and. ((accu_nd .gt. thr_nd) .or. (dabs(accu_d-dble(m))/dble(m) .gt. thr_d)) ) then
    print *, accu_nd, thr_nd 
    print *, dabs(accu_d-dble(m))/dble(m), thr_d
    print *, ' biorthog_binormalize failed !'
    stop
  endif

end subroutine check_biorthog_binormalize

! ---

subroutine check_weighted_biorthog(n, m, W, Vl, Vr, thr_d, thr_nd, accu_d, accu_nd, S, stop_ifnot)

  implicit none
  
  integer,          intent(in)  :: n, m
  double precision, intent(in)  :: Vl(n,m), Vr(n,m), W(n,n)
  double precision, intent(in)  :: thr_d, thr_nd
  logical,          intent(in)  :: stop_ifnot
  double precision, intent(out) :: accu_d, accu_nd, S(m,m)

  integer                       :: i, j
  double precision, allocatable :: SS(:,:), tmp(:,:)

  print *, ' check weighted bi-orthogonality'

  ! ---

  allocate(tmp(m,n))
  call dgemm( 'T', 'N', m, n, n, 1.d0        &
            , Vl, size(Vl, 1), W, size(W, 1) &
            , 0.d0, tmp, size(tmp, 1) )
  call dgemm( 'N', 'N', m, m, n, 1.d0            &
            , tmp, size(tmp, 1), Vr, size(Vr, 1) &
            , 0.d0, S, size(S, 1) )
  deallocate(tmp)

  !print *, ' overlap matrix:'
  !do i = 1, m
  !  write(*,'(1000(F16.10,X))') S(i,:)
  !enddo

  accu_d  = 0.d0
  accu_nd = 0.d0
  do i = 1, m
    do j = 1, m
      if(i==j) then
        accu_d = accu_d + dabs(S(i,i))
      else
        accu_nd = accu_nd + S(j,i) * S(j,i)
      endif
    enddo
  enddo
  accu_nd = dsqrt(accu_nd)

  print *, ' accu_nd = ', accu_nd
  print *, ' accu_d  = ', dabs(accu_d-dble(m))/dble(m)

  ! ---

  if( stop_ifnot .and. ((accu_nd .gt. thr_nd) .or. dabs(accu_d-dble(m))/dble(m) .gt. thr_d) ) then
    print *, ' non bi-orthogonal vectors !'
    print *, ' accu_nd = ', accu_nd
    print *, ' accu_d  = ', dabs(accu_d-dble(m))/dble(m)
    !print *, ' overlap matrix:'
    !do i = 1, m
    !  write(*,'(1000(F16.10,X))') S(i,:)
    !enddo
    stop
  endif

end subroutine check_weighted_biorthog

! ---

subroutine check_biorthog(n, m, Vl, Vr, accu_d, accu_nd, S, thr_d, thr_nd, stop_ifnot)

  implicit none
  
  integer,          intent(in)  :: n, m
  double precision, intent(in)  :: Vl(n,m), Vr(n,m)
  logical,          intent(in)  :: stop_ifnot
  double precision, intent(in)  :: thr_d, thr_nd
  double precision, intent(out) :: accu_d, accu_nd, S(m,m)

  integer                       :: i, j
  double precision, allocatable :: SS(:,:)

  !print *, ' check bi-orthogonality'

  ! ---

  call dgemm( 'T', 'N', m, m, n, 1.d0          &
            , Vl, size(Vl, 1), Vr, size(Vr, 1) &
            , 0.d0, S, size(S, 1) )

  !print *, ' overlap matrix:'
  !do i = 1, m
  !  write(*,'(1000(F16.10,X))') S(i,:)
  !enddo

  accu_d  = 0.d0
  accu_nd = 0.d0
  do i = 1, m
    do j = 1, m
      if(i==j) then
        accu_d = accu_d + dabs(S(i,i))
      else
        accu_nd = accu_nd + S(j,i) * S(j,i)
      endif
    enddo
  enddo
  accu_nd = dsqrt(accu_nd) / dble(m)

  !print *, ' accu_nd = ', accu_nd
  !print *, ' accu_d  = ', dabs(accu_d-dble(m))/dble(m)

  ! ---

  if(stop_ifnot .and. ((accu_nd .gt. thr_nd) .or. dabs(accu_d-dble(m))/dble(m) .gt. thr_d)) then
    print *, ' non bi-orthogonal vectors !'
    print *, ' accu_nd = ', accu_nd
    print *, ' accu_d  = ', dabs(accu_d-dble(m))/dble(m)
    !print *, ' overlap matrix:'
    !do i = 1, m
    !  write(*,'(1000(F16.10,X))') S(i,:)
    !enddo
    stop
  endif

end subroutine check_biorthog

! ---

subroutine check_orthog(n, m, V, accu_d, accu_nd, S)

  implicit none
  
  integer,          intent(in)  :: n, m
  double precision, intent(in)  :: V(n,m)
  double precision, intent(out) :: accu_d, accu_nd, S(m,m)

  integer                       :: i, j

  S = 0.d0
  call dgemm( 'T', 'N', m, m, n, 1.d0      &
            , V, size(V, 1), V, size(V, 1) &
            , 0.d0, S, size(S, 1) )

  !print *, ''
  !print *, ' overlap matrix:'
  !do i = 1, m
  !  write(*,'(1000(F16.10,X))') S(i,:)
  !enddo
  !print *, ''

  accu_d  = 0.d0
  accu_nd = 0.d0
  do i = 1, m
    do j = 1, m
      if(i==j) then
        accu_d = accu_d + dabs(S(i,i))
      else
        accu_nd = accu_nd + S(j,i) * S(j,i)
      endif
    enddo
  enddo
  accu_nd = dsqrt(accu_nd)

  !print*, '    diag acc: ', accu_d
  !print*, ' nondiag acc: ', accu_nd

end subroutine check_orthog

! ---

subroutine impose_biorthog_degen_eigvec(n, e0, L0, R0)

  implicit none

  integer,          intent(in)    :: n
  double precision, intent(in)    :: e0(n)
  double precision, intent(inout) :: L0(n,n), R0(n,n)

  logical                         :: complex_root
  integer                         :: i, j, k, m
  double precision                :: ei, ej, de, de_thr
  double precision                :: accu_d, accu_nd
  integer,          allocatable   :: deg_num(:)
  double precision, allocatable   :: L(:,:), R(:,:), S(:,:), S_inv_half(:,:)

  ! ---

  allocate( deg_num(n) )
  do i = 1, n
    deg_num(i) = 1
  enddo

  de_thr = thr_degen_tc

  do i = 1, n-1
    ei = e0(i)

    ! already considered in degen vectors
    if(deg_num(i).eq.0) cycle

    do j = i+1, n
      ej = e0(j)
      de = dabs(ei - ej)

      if(de .lt. de_thr) then
        deg_num(i) = deg_num(i) + 1 
        deg_num(j) = 0
      endif

    enddo
  enddo
  
  !do i = 1, n
  !  if(deg_num(i) .gt. 1) then
  !    print *, ' degen on', i, deg_num(i), e0(i)
  !  endif
  !enddo

  ! ---

  do i = 1, n
    m = deg_num(i)

    if(m .gt. 1) then
  
      allocate(L(n,m))
      allocate(R(n,m))

      do j = 1, m
        L(1:n,j) = L0(1:n,i+j-1)
        R(1:n,j) = R0(1:n,i+j-1)
      enddo

      ! ---

      call impose_orthog_svd(n, m, L)
      call impose_orthog_svd(n, m, R)
      !call impose_orthog_GramSchmidt(n, m, L)
      !call impose_orthog_GramSchmidt(n, m, R)

      ! ---

      !allocate(S(m,m))
      !call dgemm( 'T', 'N', m, m, n, 1.d0      &
      !          , L, size(L, 1), R, size(R, 1) &
      !          , 0.d0, S, size(S, 1) )
      !allocate(S_inv_half(m,m))
      !call get_inv_half_nonsymmat_diago(S, m, S_inv_half, complex_root)
      !if(complex_root) then
      !  print*, ' complex roots in inv_half !!! '
      !  stop
      !endif
      !call bi_ortho_s_inv_half(m, L, R, S_inv_half)
      !deallocate(S, S_inv_half)

      call impose_biorthog_svd(n, m, L, R)

      !call impose_biorthog_qr(n, m, thr_d, thr_nd, L, R)

      ! ---

      do j = 1, m
        L0(1:n,i+j-1) = L(1:n,j)
        R0(1:n,i+j-1) = R(1:n,j)
      enddo

      deallocate(L, R)

    endif
  enddo

end subroutine impose_biorthog_degen_eigvec 

! ---

subroutine impose_orthog_biorthog_degen_eigvec(n, thr_d, thr_nd, e0, L0, R0)

  implicit none

  integer,          intent(in)    :: n
  double precision, intent(in)    :: thr_d, thr_nd
  double precision, intent(in)    :: e0(n)
  double precision, intent(inout) :: L0(n,n), R0(n,n)

  integer                         :: i, j, k, m
  double precision                :: ei, ej, de, de_thr
  double precision                :: accu_d, accu_nd
  integer,          allocatable   :: deg_num(:)
  double precision, allocatable   :: L(:,:), R(:,:), S(:,:)

  ! ---

  allocate( deg_num(n) )
  do i = 1, n
    deg_num(i) = 1
  enddo

  de_thr = thr_degen_tc

  do i = 1, n-1
    ei = e0(i)

    ! already considered in degen vectors
    if(deg_num(i).eq.0) cycle

    do j = i+1, n
      ej = e0(j)
      de = dabs(ei - ej)

      if(de .lt. de_thr) then
        deg_num(i) = deg_num(i) + 1 
        deg_num(j) = 0
      endif

    enddo
  enddo
  
  do i = 1, n
    if(deg_num(i).gt.1) then
      print *, ' degen on', i, deg_num(i)
    endif
  enddo

  ! ---

  do i = 1, n
    m = deg_num(i)

    if(m .gt. 1) then
  
      allocate(L(n,m))
      allocate(R(n,m))

      do j = 1, m
        L(1:n,j) = L0(1:n,i+j-1)
        R(1:n,j) = R0(1:n,i+j-1)
      enddo

      ! ---

      call impose_orthog_svd(n, m, L)
      call impose_orthog_svd(n, m, R)

      ! ---
  
      call impose_biorthog_qr(n, m, thr_d, thr_nd, L, R)

      allocate(S(m,m))
      call check_biorthog(n, m, L, R, accu_d, accu_nd, S, thr_d, thr_nd, .true.)
      !call check_biorthog(n, m, L, L, accu_d, accu_nd, S, thr_d, thr_nd, .true.)
      !call check_biorthog(n, m, R, R, accu_d, accu_nd, S, thr_d, thr_nd, .false.)
      deallocate(S)

      ! ---

      do j = 1, m
        L0(1:n,i+j-1) = L(1:n,j)
        R0(1:n,i+j-1) = R(1:n,j)
      enddo

      deallocate(L, R)

    endif
  enddo

end subroutine impose_orthog_biorthog_degen_eigvec 

! ---

subroutine impose_unique_biorthog_degen_eigvec(n, thr_d, thr_nd, e0, C0, W0, L0, R0)

  implicit none

  integer,          intent(in)    :: n
  double precision, intent(in)    :: thr_d, thr_nd
  double precision, intent(in)    :: e0(n), W0(n,n), C0(n,n)
  double precision, intent(inout) :: L0(n,n), R0(n,n)

  logical                         :: complex_root
  integer                         :: i, j, k, m
  double precision                :: ei, ej, de, de_thr
  integer,          allocatable   :: deg_num(:)
  double precision, allocatable   :: L(:,:), R(:,:), C(:,:)
  double precision, allocatable   :: S(:,:), S_inv_half(:,:), tmp(:,:)

  ! ---

  allocate( deg_num(n) )
  do i = 1, n
    deg_num(i) = 1
  enddo

  de_thr = thr_degen_tc

  do i = 1, n-1
    ei = e0(i)

    ! already considered in degen vectors
    if(deg_num(i).eq.0) cycle

    do j = i+1, n
      ej = e0(j)
      de = dabs(ei - ej)

      if(de .lt. de_thr) then
        deg_num(i) = deg_num(i) + 1 
        deg_num(j) = 0
      endif

    enddo
  enddo
  
  !do i = 1, n
  !  if(deg_num(i) .gt. 1) then
  !    print *, ' degen on', i, deg_num(i)
  !  endif
  !enddo

  ! ---

  do i = 1, n
    m = deg_num(i)

    if(m .gt. 1) then
  
      allocate(L(n,m))
      allocate(R(n,m))
      allocate(C(n,m))

      do j = 1, m
        L(1:n,j) = L0(1:n,i+j-1)
        R(1:n,j) = R0(1:n,i+j-1)
        C(1:n,j) = C0(1:n,i+j-1)
      enddo

      ! ---

      call impose_orthog_svd(n, m, L)
      call impose_orthog_svd(n, m, R)

      ! ---


      ! TODO:
      ! select C correctly via overlap
      ! or via selecting degen in HF

      !call max_overlap_qr(n, m, C, L)
      !call max_overlap_qr(n, m, C, R)


      allocate(tmp(m,n))
      allocate(S(m,m))

      call dgemm( 'T', 'N', m, n, n, 1.d0        &
                , L, size(L, 1), W0, size(W0, 1) &
                , 0.d0, tmp, size(tmp, 1) )
      call dgemm( 'N', 'N', m, m, n, 1.d0          &
                , tmp, size(tmp, 1), C, size(C, 1) &
                , 0.d0, S, size(S, 1) )

      call max_overlap_qr(n, m, S, L)
      !call max_overlap_invprod(n, m, S, L)

      call dgemm( 'T', 'N', m, n, n, 1.d0        &
                , C, size(C, 1), W0, size(W0, 1) &
                , 0.d0, tmp, size(tmp, 1) )
      call dgemm( 'N', 'N', m, m, n, 1.d0          &
                , tmp, size(tmp, 1), R, size(R, 1) &
                , 0.d0, S, size(S, 1) )

      call max_overlap_qr(n, m, S, R)
      !call max_overlap_invprod(n, m, S, R)

      deallocate(S, tmp)

      ! ---
  
      allocate(S(m,m), S_inv_half(m,m))
      call dgemm( 'T', 'N', m, m, n, 1.d0      &
                , L, size(L, 1), R, size(R, 1) &
                , 0.d0, S, size(S, 1) )
      call get_inv_half_nonsymmat_diago(S, m, S_inv_half, complex_root)
      if(complex_root)then
        call impose_biorthog_svd(n, m, L, R)
        !call impose_biorthog_qr(n, m, thr_d, thr_nd, L, R)
      else
        call bi_ortho_s_inv_half(m, L, R, S_inv_half)
      endif
      deallocate(S, S_inv_half)

      ! ---

      do j = 1, m
        L0(1:n,i+j-1) = L(1:n,j)
        R0(1:n,i+j-1) = R(1:n,j)
      enddo

      deallocate(L, R, C)

    endif
  enddo

end subroutine impose_unique_biorthog_degen_eigvec 

! ---

subroutine max_overlap_qr(m, n, S0, V)

  implicit none 
  integer,          intent(in)    :: m, n
  double precision, intent(in)    :: S0(n,n)
  double precision, intent(inout) :: V(m,n)

  integer                         :: i, j
  integer                         :: LWORK, INFO
  double precision, allocatable   :: TAU(:), WORK(:)
  double precision, allocatable   :: S(:,:), tmp(:,:)

  allocate(S(n,n))
  S = S0

  ! ---

  allocate( TAU(n), WORK(1) )

  LWORK = -1
  call dgeqrf(n, n, S, n, TAU, WORK, LWORK, INFO)
  if(INFO .ne. 0) then
    print*,'dgeqrf failed !!', INFO
    stop
  endif

  LWORK = max(n, int(WORK(1)))
  deallocate(WORK)

  allocate( WORK(LWORK) )
  call dgeqrf(n, n, S, n, TAU, WORK, LWORK, INFO)
  if(INFO .ne. 0) then
    print*,'dgeqrf failed !!', INFO
    stop
  endif

  ! get Q in S matrix
  LWORK = -1
  call dorgqr(n, n, n, S, n, TAU, WORK, LWORK, INFO)
  if(INFO .ne. 0) then
    print*,'dorgqr failed !!', INFO
    stop
  endif

  LWORK = max(n, int(WORK(1)))
  deallocate(WORK)

  allocate( WORK(LWORK) )
  call dorgqr(n, n, n, S, n, TAU, WORK, LWORK, INFO)
  if(INFO .ne. 0) then
    print*,'dorgqr failed !!', INFO
    stop
  endif

  deallocate( WORK, TAU )

  ! ---

  ! V0.T <-- Q.T x V0.T, where Q = S

  allocate( tmp(n,m) )

  call dgemm( 'T', 'T', n, m, n, 1.d0      &
            , S, size(S, 1), V, size(V, 1) &
            , 0.d0, tmp, size(tmp, 1) )

  deallocate(S)

  do i = 1, n
    do j = 1, m
      V(j,i) = tmp(i,j)
    enddo
  enddo

  deallocate(tmp)

  ! ---

  return
end subroutine max_overlap_qr

! ---

subroutine max_overlap_invprod(n, m, S, V)

  implicit none 
  integer,          intent(in)    :: m, n
  double precision, intent(in)    :: S(m,m)
  double precision, intent(inout) :: V(n,m)

  integer                         :: i
  double precision, allocatable   :: invS(:,:), tmp(:,:)

  allocate(invS(m,m))
  call get_inverse(S, size(S, 1), m, invS, size(invS, 1))
  print *, ' overlap '
  do i = 1, m
    write(*, '(1000(F16.10,X))') S(i,:)
  enddo
  print *, ' inv overlap '
  do i = 1, m
    write(*, '(1000(F16.10,X))') invS(i,:)
  enddo

  allocate(tmp(n,m))
  tmp = V

  call dgemm( 'N', 'N', n, m, m, 1.d0        &
            , tmp, size(tmp, 1), invS, size(invS, 1) &
            , 0.d0, V, size(V, 1) )

  deallocate(tmp, invS)

  return
end subroutine max_overlap_invprod

! ---

subroutine impose_biorthog_svd(n, m, L, R)

  implicit none

  integer,          intent(in)    :: n, m
  double precision, intent(inout) :: L(n,m), R(n,m)

  integer                         :: i, j, num_linear_dependencies
  double precision                :: threshold
  double precision, allocatable   :: S(:,:), tmp(:,:)
  double precision, allocatable   :: U(:,:), V(:,:), Vt(:,:), D(:)

  ! ---

  allocate(S(m,m))

  call dgemm( 'T', 'N', m, m, n, 1.d0      &
            , L, size(L, 1), R, size(R, 1) &
            , 0.d0, S, size(S, 1) )

  !print *, ' overlap bef SVD: '
  !do i = 1, m
  !  write(*, '(1000(F16.10,X))') S(i,:)
  !enddo

  ! ---
 
  allocate(U(m,m), Vt(m,m), D(m))

  call svd(S, m, U, m, D, Vt, m, m, m)

  deallocate(S)

  threshold               = 1.d-6
  num_linear_dependencies = 0
  do i = 1, m
    if(abs(D(i)) <= threshold) then
      D(i) = 0.d0
      num_linear_dependencies += 1
    else
      ASSERT (D(i) > 0.d0)
      D(i) = 1.d0 / dsqrt(D(i))
    endif
  enddo
  if(num_linear_dependencies > 0) then
    write(*,*) ' linear dependencies = ', num_linear_dependencies
    write(*,*) ' m                   = ', m
    stop
  endif

  allocate(V(m,m))
  do i = 1, m
    do j = 1, m
      V(j,i) = Vt(i,j)
    enddo
  enddo
  deallocate(Vt)

  ! ---

  allocate(tmp(n,m))

  ! tmp <-- R x V
  call dgemm( 'N', 'N', n, m, m, 1.d0      &
            , R, size(R, 1), V, size(V, 1) &
            , 0.d0, tmp, size(tmp, 1) )
  deallocate(V)
  ! R <-- tmp x sigma^-0.5
  do j = 1, m
    do i = 1, n
      R(i,j) = tmp(i,j) * D(j)
    enddo
  enddo

  ! tmp <-- L x U
  call dgemm( 'N', 'N', n, m, m, 1.d0      &
            , L, size(L, 1), U, size(U, 1) &
            , 0.d0, tmp, size(tmp, 1) )
  deallocate(U)
  ! L <-- tmp x sigma^-0.5
  do j = 1, m
    do i = 1, n
      L(i,j) = tmp(i,j) * D(j)
    enddo
  enddo

  deallocate(D, tmp)

  ! ---

  allocate(S(m,m))
  call dgemm( 'T', 'N', m, m, n, 1.d0      &
            , L, size(L, 1), R, size(R, 1) &
            , 0.d0, S, size(S, 1) )

  !print *, ' overlap aft SVD: '
  !do i = 1, m
  !  write(*, '(1000(F16.10,X))') S(i,:)
  !enddo

  deallocate(S)

  ! ---

end subroutine impose_biorthog_svd

! ---

subroutine impose_weighted_biorthog_qr(m, n, thr_d, thr_nd, Vl, W, Vr)

  implicit none 
  integer,          intent(in)    :: m, n
  double precision, intent(in)    :: thr_d, thr_nd
  double precision, intent(inout) :: Vl(m,n), W(m,m), Vr(m,n)

  integer                         :: i, j
  integer                         :: LWORK, INFO
  double precision                :: accu_nd, accu_d
  double precision, allocatable   :: TAU(:), WORK(:)
  double precision, allocatable   :: S(:,:), R(:,:), tmp(:,:), Stmp(:,:)


  call check_weighted_biorthog_binormalize(m, n, Vl, W, Vr, thr_d, thr_nd, .false.)
  
  ! ---
  
  allocate(Stmp(n,m), S(n,n))
  call dgemm( 'T', 'N', n, m, m, 1.d0        &
            , Vl, size(Vl, 1), W, size(W, 1) &
            , 0.d0, Stmp, size(Stmp, 1) )
  call dgemm( 'N', 'N', n, n, m, 1.d0              &
            , Stmp, size(Stmp, 1), Vr, size(Vr, 1) &
            , 0.d0, S, size(S, 1) )
  deallocate(Stmp)

  accu_nd = 0.d0
  accu_d  = 0.d0
  do i = 1, n
    do j = 1, n
      if(i==j) then
        accu_d += S(j,i)
      else
        accu_nd = accu_nd + S(j,i) * S(j,i)
      endif
    enddo
  enddo
  accu_nd = dsqrt(accu_nd)

  if((accu_nd .lt. thr_nd) .and. (dabs(accu_d-dble(n))/dble(n) .lt. thr_d)) then
    print *, ' bi-orthogonal vectors without QR !'
    deallocate(S)
    return
  endif

  ! -------------------------------------------------------------------------------------
  !                           QR factorization of S: S = Q x R


  print *, ' apply QR decomposition ...'

  allocate( TAU(n), WORK(1) )

  LWORK = -1
  call dgeqrf(n, n, S, n, TAU, WORK, LWORK, INFO)
  if(INFO .ne. 0) then
    print*,'dgeqrf failed !!', INFO
    stop
  endif

  LWORK = max(n, int(WORK(1)))
  deallocate(WORK)

  allocate( WORK(LWORK) )
  call dgeqrf(n, n, S, n, TAU, WORK, LWORK, INFO)
  if(INFO .ne. 0) then
    print*,'dgeqrf failed !!', INFO
    stop
  endif

  ! save the upper triangular R
  allocate( R(n,n) )
  R(:,:) = S(:,:)

  ! get Q
  LWORK = -1
  call dorgqr(n, n, n, S, n, TAU, WORK, LWORK, INFO)
  if(INFO .ne. 0) then
    print*,'dorgqr failed !!', INFO
    stop
  endif

  LWORK = max(n, int(WORK(1)))
  deallocate(WORK)

  allocate( WORK(LWORK) )
  call dorgqr(n, n, n, S, n, TAU, WORK, LWORK, INFO)
  if(INFO .ne. 0) then
    print*,'dorgqr failed !!', INFO
    stop
  endif

  deallocate( WORK, TAU )

  !
  ! -------------------------------------------------------------------------------------

  ! ---

  ! -------------------------------------------------------------------------------------
  !                               get bi-orhtog left & right vectors:
  !                                           Vr' = Vr x inv(R) 
  !                                           Vl' = inv(Q) x Vl =  Q.T   x Vl 

  ! Q.T x Vl, where Q = S

  allocate( tmp(n,m) )
  call dgemm( 'T', 'T', n, m, n, 1.d0        &
            , S, size(S, 1), Vl, size(Vl, 1) &
            , 0.d0, tmp, size(tmp, 1) )

  do i = 1, n
    do j = 1, m
      Vl(j,i) = tmp(i,j)
    enddo
  enddo
  deallocate(tmp)

  ! ---

  ! inv(R) 
  !print *, ' inversing upper triangular matrix ...'
  call dtrtri("U", "N", n, R, n, INFO)
  if(INFO .ne. 0) then
    print*,'dtrtri failed !!', INFO
    stop
  endif
  !print *, ' inversing upper triangular matrix OK' 

  do i = 1, n-1
    do j = i+1, n
      R(j,i) = 0.d0
    enddo
  enddo

  !print *, ' inv(R):'
  !do i = 1, n
  !  write(*, '(1000(F16.10,X))') R(i,:)
  !enddo

  ! Vr x inv(R) 
  allocate( tmp(m,n) )
  call dgemm( 'N', 'N', m, n, n, 1.d0        &
            , Vr, size(Vr, 1), R, size(R, 1) &
            , 0.d0, tmp, size(tmp, 1) )
  deallocate( R )

  do i = 1, n
    do j = 1, m
      Vr(j,i) = tmp(j,i)
    enddo
  enddo
  deallocate(tmp)

  call check_weighted_biorthog_binormalize(m, n, Vl, W, Vr, thr_d, thr_nd, .false.)

  return
end subroutine impose_weighted_biorthog_qr

! ---

subroutine check_weighted_biorthog_binormalize(n, m, Vl, W, Vr, thr_d, thr_nd, stop_ifnot)

  implicit none
  
  integer,          intent(in)    :: n, m
  logical,          intent(in)    :: stop_ifnot
  double precision, intent(in)    :: thr_d, thr_nd
  double precision, intent(inout) :: Vl(n,m), W(n,n), Vr(n,m)

  integer                         :: i, j
  double precision                :: accu_d, accu_nd, s_tmp
  double precision, allocatable   :: S(:,:), Stmp(:,:)

  print *, ' check weighted bi-orthonormality'

  ! ---

  allocate(Stmp(m,n), S(m,m))
  call dgemm( 'T', 'N', m, n, n, 1.d0        &
            , Vl, size(Vl, 1), W, size(W, 1) &
            , 0.d0, Stmp, size(Stmp, 1) )
  call dgemm( 'N', 'N', m, m, n, 1.d0              &
            , Stmp, size(Stmp, 1), Vr, size(Vr, 1) &
            , 0.d0, S, size(S, 1) )
  deallocate(Stmp)
  !print *, ' overlap matrix before:'
  !do i = 1, m
  !  write(*,'(1000(F16.10,X))') S(i,:)
  !enddo

  ! S(i,i) = -1
  do i = 1, m
    if( (S(i,i) + 1.d0) .lt. thr_d ) then
      do j = 1, n
        Vl(j,i) = -1.d0 * Vl(j,i)
      enddo
      S(i,i) = 1.d0
    endif
  enddo

  accu_d  = 0.d0
  accu_nd = 0.d0
  do i = 1, m
    do j = 1, m
      if(i==j) then
        accu_d = accu_d + S(i,i)
      else
        accu_nd = accu_nd + S(j,i) * S(j,i)
      endif
    enddo
  enddo
  accu_nd = dsqrt(accu_nd) / dble(m)
  print*, '    diag acc: ', accu_d
  print*, ' nondiag acc: ', accu_nd

  ! ---

  if( (accu_nd .lt. thr_nd) .and. (dabs(accu_d-dble(m))/dble(m) .gt. thr_d) ) then

    do i = 1, m
      print *, i, S(i,i)
      if(dabs(S(i,i) - 1.d0) .gt. thr_d) then
        s_tmp = 1.d0 / dsqrt(S(i,i))
        do j = 1, n
          Vl(j,i) = Vl(j,i) * s_tmp 
          Vr(j,i) = Vr(j,i) * s_tmp 
        enddo
      endif
    enddo

  endif

  ! ---

  allocate(Stmp(m,n))
  call dgemm( 'T', 'N', m, n, n, 1.d0        &
            , Vl, size(Vl, 1), W, size(W, 1) &
            , 0.d0, Stmp, size(Stmp, 1) )
  call dgemm( 'N', 'N', m, m, n, 1.d0              &
            , Stmp, size(Stmp, 1), Vr, size(Vr, 1) &
            , 0.d0, S, size(S, 1) )
  deallocate(Stmp)
  !print *, ' overlap matrix after:'
  !do i = 1, m
  !  write(*,'(1000(F16.10,X))') S(i,:)
  !enddo

  accu_d  = 0.d0
  accu_nd = 0.d0
  do i = 1, m
    do j = 1, m
      if(i==j) then
        accu_d = accu_d + S(i,i)
      else
        accu_nd = accu_nd + S(j,i) * S(j,i)
      endif
    enddo
  enddo
  accu_nd = dsqrt(accu_nd) / dble(m)
  print *, '    diag acc: ', accu_d
  print *, ' nondiag acc: ', accu_nd

  deallocate(S)

  ! ---

  if( stop_ifnot .and. ((accu_nd .gt. thr_nd) .or. (dabs(accu_d-dble(m))/dble(m) .gt. thr_d)) ) then
    print *, accu_nd, thr_nd 
    print *, dabs(accu_d-dble(m))/dble(m), thr_d
    print *, ' weighted biorthog_binormalize failed !'
    stop
  endif

end subroutine check_weighted_biorthog_binormalize

! ---

subroutine impose_weighted_biorthog_svd(n, m, overlap, L, R)

  implicit none

  integer,          intent(in)    :: n, m
  double precision, intent(in)    :: overlap(n,n)
  double precision, intent(inout) :: L(n,m), R(n,m)

  integer                         :: i, j, num_linear_dependencies
  double precision                :: threshold
  double precision, allocatable   :: S(:,:), tmp(:,:),Stmp(:,:)
  double precision, allocatable   :: U(:,:), V(:,:), Vt(:,:), D(:)

  ! ---

  allocate(S(m,m),Stmp(n,m))

  ! S = C.T x overlap x C
  call dgemm( 'N', 'N', n, m, n, 1.d0      &
            , overlap, size(overlap, 1), R, size(R, 1) &
            , 0.d0, Stmp, size(Stmp, 1) )
  call dgemm( 'T', 'N', m, m, n, 1.d0      &
            , L, size(L, 1), Stmp, size(Stmp, 1) &
            , 0.d0, S, size(S, 1) )
  deallocate(Stmp)

  !print *, ' overlap bef SVD: '
  !do i = 1, m
  !  write(*, '(1000(F25.16,X))') S(i,:)
  !enddo

  ! ---
 
  allocate(U(m,m), Vt(m,m), D(m))

  call svd(S, m, U, m, D, Vt, m, m, m)

  deallocate(S)

  threshold               = 1.d-6
  num_linear_dependencies = 0
  do i = 1, m
    if(abs(D(i)) <= threshold) then
      D(i) = 0.d0
      num_linear_dependencies += 1
    else
      ASSERT (D(i) > 0.d0)
      D(i) = 1.d0 / dsqrt(D(i))
    endif
  enddo
  if(num_linear_dependencies > 0) then
    write(*,*) ' linear dependencies = ', num_linear_dependencies
    write(*,*) ' m                   = ', m
    stop
  endif

  allocate(V(m,m))
  do i = 1, m
    do j = 1, m
      V(j,i) = Vt(i,j)
    enddo
  enddo
  deallocate(Vt)

  ! ---

  allocate(tmp(n,m))

  ! tmp <-- R x V
  call dgemm( 'N', 'N', n, m, m, 1.d0      &
            , R, size(R, 1), V, size(V, 1) &
            , 0.d0, tmp, size(tmp, 1) )
  deallocate(V)
  ! R <-- tmp x sigma^-0.5
  do j = 1, m
    do i = 1, n
      R(i,j) = tmp(i,j) * D(j)
    enddo
  enddo

  ! tmp <-- L x U
  call dgemm( 'N', 'N', n, m, m, 1.d0      &
            , L, size(L, 1), U, size(U, 1) &
            , 0.d0, tmp, size(tmp, 1) )
  deallocate(U)
  ! L <-- tmp x sigma^-0.5
  do j = 1, m
    do i = 1, n
      L(i,j) = tmp(i,j) * D(j)
    enddo
  enddo

  deallocate(D, tmp)

  ! ---

  allocate(S(m,m),Stmp(n,m))
  ! S = C.T x overlap x C
  call dgemm( 'N', 'N', n, m, n, 1.d0      &
            , overlap, size(overlap, 1), R, size(R, 1) &
            , 0.d0, Stmp, size(Stmp, 1) )
  call dgemm( 'T', 'N', m, m, n, 1.d0      &
            , L, size(L, 1), Stmp, size(Stmp, 1) &
            , 0.d0, S, size(S, 1) )
  deallocate(Stmp)

  !print *, ' overlap aft SVD with overlap: '
  !do i = 1, m
  !  write(*, '(1000(F16.10,X))') S(i,:)
  !enddo

  deallocate(S)

  return
end subroutine impose_weighted_biorthog_svd

! ---

