
subroutine get_inv_half_svd(matrix, n, matrix_inv_half)

  BEGIN_DOC
  !   :math:`X = S^{-1/2}` obtained by SVD
  END_DOC

  implicit none

  integer,          intent(in)  :: n
  double precision, intent(in)  :: matrix(n,n)
  double precision, intent(out) :: matrix_inv_half(n,n)

  integer                       :: num_linear_dependencies
  integer                       :: LDA, LDC
  integer                       :: info, i, j, k
  double precision, parameter   :: threshold = 1.d-6
  double precision, allocatable :: U(:,:),Vt(:,:), D(:),matrix_half(:,:),D_half(:)

  double precision :: accu_d,accu_nd

  LDA = size(matrix, 1)
  LDC = size(matrix_inv_half, 1)
  if(LDA .ne. LDC) then
    print*, ' LDA != LDC'
    stop
  endif

  print*, ' n   = ', n
  print*, ' LDA = ', LDA
  print*, ' LDC = ', LDC

  double precision,allocatable :: WR(:),WI(:),VL(:,:),VR(:,:)
  allocate(WR(n),WI(n),VL(n,n),VR(n,n))
  call lapack_diag_non_sym(n,matrix,WR,WI,VL,VR)
  do i = 1, n
    print*,'WR,WI',WR(i),WI(i)
  enddo


  allocate(U(LDC,n), Vt(LDA,n), D(n))

  call svd(matrix, LDA, U, LDC, D, Vt, LDA, n, n)
  double precision, allocatable :: tmp1(:,:),tmp2(:,:),D_mat(:,:)
  allocate(tmp1(n,n),tmp2(n,n),D_mat(n,n),matrix_half(n,n),D_half(n))
  D_mat = 0.d0
  do i = 1,n
   D_mat(i,i) = D(i)
  enddo
  ! matrix = U D Vt
  ! tmp1 = U D
  tmp1 = 0.d0
  call dgemm( 'N', 'N', n, n, n, 1.d0                                            &
            , U, size(U, 1), D_mat, size(D_mat, 1) &
            , 0.d0, tmp1, size(tmp1, 1) )
  ! tmp2 = tmp1 X Vt = matrix
  tmp2 = 0.d0
  call dgemm( 'N', 'N', n, n, n, 1.d0                                            &
            , tmp1, size(tmp1, 1), Vt, size(Vt, 1) &
            , 0.d0, tmp2, size(tmp2, 1) )
  print*,'Checking the recomposition of the matrix'
  accu_nd = 0.d0
  accu_d  = 0.d0
  do i = 1, n
   accu_d += dabs(tmp2(i,i) - matrix(i,i))
   do j = 1, n
    if(i==j)cycle 
    accu_nd += dabs(tmp2(j,i) - matrix(j,i))
   enddo
  enddo
  print*,'accu_d  =',accu_d
  print*,'accu_nd =',accu_nd
  print*,'passed the recomposition'

  num_linear_dependencies = 0
  do i = 1, n
    if(abs(D(i)) <= threshold) then
      D(i) = 0.d0
      num_linear_dependencies += 1
    else
      ASSERT (D(i) > 0.d0)
      D_half(i) = dsqrt(D(i))
      D(i) = 1.d0 / dsqrt(D(i))
    endif
  enddo
  write(*,*) ' linear dependencies', num_linear_dependencies

  matrix_inv_half = 0.d0
  matrix_half = 0.d0
  do k = 1, n
    if(D(k) /= 0.d0) then
      do j = 1, n
        do i = 1, n
!          matrix_inv_half(i,j) = matrix_inv_half(i,j) + U(i,k) * D(k) * Vt(k,j)
           matrix_inv_half(i,j) = matrix_inv_half(i,j) + U(i,k) * D(k) * Vt(j,k)
           matrix_half(i,j) = matrix_half(i,j) + U(i,k) * D_half(k) * Vt(j,k)
        enddo
      enddo
    endif
  enddo
  print*,'testing S^1/2 * S^1/2= S'
  ! tmp1 = S^1/2 X S^1/2
  tmp1 = 0.d0
  call dgemm( 'N', 'N', n, n, n, 1.d0                                            &
            , matrix_half, size(matrix_half, 1), matrix_half, size(matrix_half, 1) &
            , 0.d0, tmp1, size(tmp1, 1) )
  accu_nd = 0.d0
  accu_d  = 0.d0
  do i = 1, n
   accu_d += dabs(tmp1(i,i) - matrix(i,i))
   do j = 1, n
    if(i==j)cycle 
    accu_nd += dabs(tmp1(j,i) - matrix(j,i))
   enddo
  enddo
  print*,'accu_d  =',accu_d
  print*,'accu_nd =',accu_nd

!  print*,'S inv half'
!  do i = 1, n
!   write(*, '(1000(F16.10,X))') matrix_inv_half(i,:)
!  enddo

  double precision, allocatable :: pseudo_inverse(:,:),identity(:,:)
  allocate( pseudo_inverse(n,n),identity(n,n))
  call get_pseudo_inverse(matrix,n,n,n,pseudo_inverse,n,threshold)
 
  ! S^-1 X S = 1
!  identity = 0.d0
!  call dgemm( 'N', 'N', n, n, n, 1.d0                                            &
!            , matrix, size(matrix, 1), pseudo_inverse, size(pseudo_inverse, 1) &
!            , 0.d0, identity, size(identity, 1) )
  print*,'Checking  S^-1/2 X S^-1/2 = S^-1 ?'
  ! S^-1/2 X S^-1/2 = S^-1 ?
  tmp1 = 0.d0
  call dgemm( 'N', 'N', n, n, n, 1.d0                                            &
            ,matrix_inv_half, size(matrix_inv_half, 1), matrix_inv_half, size(matrix_inv_half, 1) &
            , 0.d0, tmp1, size(tmp1, 1) )
  accu_nd = 0.d0
  accu_d  = 0.d0
  do i = 1, n
   accu_d += dabs(1.d0 - pseudo_inverse(i,i))
   do j = 1, n
    if(i==j)cycle 
    accu_nd += dabs(tmp1(j,i) - pseudo_inverse(j,i))
   enddo
  enddo
  print*,'accu_d  =',accu_d
  print*,'accu_nd =',accu_nd

  stop
!
!  ! ( S^-1/2 x S ) x S^-1/2
!  Stmp2 = 0.d0
!  call dgemm( 'N', 'N', n, n, n, 1.d0                                        &
!            , Stmp, size(Stmp, 1), matrix_inv_half, size(matrix_inv_half, 1) &
!            , 0.d0, Stmp2, size(Stmp2, 1) )

  ! S^-1/2 x ( S^-1/2 x S )
!  Stmp2 = 0.d0
!  call dgemm( 'N', 'N', n, n, n, 1.d0                                        &
!            , matrix_inv_half, size(matrix_inv_half, 1), Stmp, size(Stmp, 1) &
!            , 0.d0, Stmp2, size(Stmp2, 1) )
 
!  do i = 1, n
!    do j = 1, n
!      if(i==j) then
!       accu_d += Stmp2(j,i)
!      else 
!       accu_nd = accu_nd + Stmp2(j,i) * Stmp2(j,i)
!      endif
!    enddo
!  enddo
!  accu_nd = dsqrt(accu_nd)
!  print*, ' after S^-1/2: sum of off-diag S elements = ', accu_nd
!  print*, ' after S^-1/2: sum of     diag S elements = ', accu_d
!  do i = 1, n
!   write(*,'(1000(F16.10,X))') Stmp2(i,:)
!  enddo

  !double precision  :: thresh
  !thresh = 1.d-10
  !if( accu_nd.gt.thresh .or. dabs(accu_d-dble(n)).gt.thresh) then
  !  stop
  !endif

end subroutine get_inv_half_svd

! ---

subroutine get_inv_half_nonsymmat_diago(matrix, n, matrix_inv_half, complex_root)

  BEGIN_DOC
  ! input:  S = matrix
  ! output: S^{-1/2} = matrix_inv_half obtained by diagonalization
  !
  ! S = VR D VL^T
  !   = VR D^{1/2} D^{1/2} VL^T
  !   = VR D^{1/2} VL^T VR D^{1/2} VL^T
  !   = S^{1/2} S^{1/2} with S = VR D^{1/2} VL^T 
  !
  ! == > S^{-1/2} = VR D^{-1/2} VL^T
  !
  END_DOC

  implicit none

  integer,          intent(in)  :: n
  double precision, intent(in)  :: matrix(n,n)
  logical,          intent(out) :: complex_root
  double precision, intent(out) :: matrix_inv_half(n,n)

  integer                       :: i, j
  double precision              :: accu_d, accu_nd
  double precision, allocatable :: WR(:), WI(:), VL(:,:), VR(:,:), S(:,:), S_diag(:)
  double precision, allocatable :: tmp1(:,:), D_mat(:,:)

  complex_root = .False.

  matrix_inv_half = 0.D0
  print*,'Computing S^{-1/2}'

  allocate(WR(n), WI(n), VL(n,n), VR(n,n))
  call lapack_diag_non_sym(n, matrix, WR, WI, VL, VR)

  allocate(S(n,n))
  call check_biorthog(n, n, VL, VR, accu_d, accu_nd, S)
  print*,'accu_nd S^{-1/2}',accu_nd
  if(accu_nd.gt.1.d-10) then
    complex_root = .True. ! if vectors are not bi-orthogonal return 
    print*,'Eigenvectors of S are not bi-orthonormal, skipping S^{-1/2}'
    return
  endif

  allocate(S_diag(n))
  do i = 1, n
    S_diag(i) = 1.d0/dsqrt(S(i,i))
    if(dabs(WI(i)).gt.1.d-20.or.WR(i).lt.0.d0)then ! check that eigenvalues are real and positive 
     complex_root = .True.
     print*,'Eigenvalues of S have imaginary part '
     print*,'WR(i),WI(i)',WR(i), WR(i)
     print*,'Skipping S^{-1/2}'
     return
    endif
  enddo
  deallocate(S)

  if(complex_root) return

  ! normalization of vectors 
  do i = 1, n
    if(S_diag(i).eq.1.d0) cycle
    do j = 1,n
      VL(j,i) *= S_diag(i)
      VR(j,i) *= S_diag(i)
    enddo
  enddo
  deallocate(S_diag)

  allocate(tmp1(n,n), D_mat(n,n))

  D_mat = 0.d0
  do i = 1, n
    D_mat(i,i) = 1.d0/dsqrt(WR(i))
  enddo
  deallocate(WR, WI)

  ! tmp1 = VR D^{-1/2} 
  tmp1 = 0.d0
  call dgemm( 'N', 'N', n, n, n, 1.d0                &
            , VR, size(VR, 1), D_mat, size(D_mat, 1) &
            , 0.d0, tmp1, size(tmp1, 1) )
  deallocate(VR, D_mat)

  ! S^{-1/2} = tmp1 X VL^T 
  matrix_inv_half = 0.d0
  call dgemm( 'N', 'T', n, n, n, 1.d0              &
            , tmp1, size(tmp1, 1), VL, size(VL, 1) &
            , 0.d0, matrix_inv_half, size(matrix_inv_half, 1) )
  deallocate(tmp1, VL)

end

! ---

subroutine bi_ortho_s_inv_half(n,leigvec,reigvec,S_nh_inv_half)
 implicit  none
 integer, intent(in) :: n
 double precision, intent(in) :: S_nh_inv_half(n,n)
 double precision, intent(inout) :: leigvec(n,n),reigvec(n,n)
 BEGIN_DOC 
 ! bi-orthonormalization of left and right vectors 
 ! 
 ! S = VL^T VR 
 !
 ! S^{-1/2} S S^{-1/2} = 1 = S^{-1/2} VL^T VR S^{-1/2} = VL_new^T VR_new
 !
 ! VL_new = VL (S^{-1/2})^T
 ! 
 ! VR_new = VR S^{^{-1/2}}
 END_DOC
 double precision,allocatable :: vl_tmp(:,:),vr_tmp(:,:)
 print*,'Bi-orthonormalization using S^{-1/2}'
 allocate(vl_tmp(n,n),vr_tmp(n,n))
 vl_tmp = leigvec
 vr_tmp = reigvec
 ! VL_new = VL (S^{-1/2})^T
 call dgemm( 'N', 'T', n, n, n, 1.d0                                            &
             , vl_tmp, size(vl_tmp, 1), S_nh_inv_half, size(S_nh_inv_half, 1) &
             , 0.d0, leigvec, size(leigvec, 1) )
 ! VR_new = VR S^{^{-1/2}}
 call dgemm( 'N', 'N', n, n, n, 1.d0                                            &
             , vr_tmp, size(vr_tmp, 1), S_nh_inv_half, size(S_nh_inv_half, 1) &
             , 0.d0, reigvec, size(reigvec, 1) )
 double precision :: accu_d, accu_nd
 double precision,allocatable :: S(:,:) 
 allocate(S(n,n)) 
 call check_biorthog(n, n, leigvec, reigvec, accu_d, accu_nd, S)
 if(dabs(accu_d - n).gt.1.d-10 .or. accu_nd .gt.1.d-8 )then
  print*,'Pb in bi_ortho_s_inv_half !!'
  print*,'accu_d =',accu_d
  print*,'accu_nd =',accu_nd
  stop
 endif
end
