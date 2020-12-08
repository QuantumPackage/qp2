subroutine svd(A,LDA,U,LDU,D,Vt,LDVt,m,n)
  implicit none
  BEGIN_DOC
  ! Compute A = U.D.Vt
  !
  ! LDx : leftmost dimension of x
  !
  ! Dimension of A is m x n
  !
  END_DOC

  integer, intent(in)             :: LDA, LDU, LDVt, m, n
  double precision, intent(in)    :: A(LDA,n)
  double precision, intent(out)   :: U(LDU,min(m,n))
  double precision,intent(out)    :: Vt(LDVt,n)
  double precision,intent(out)    :: D(min(m,n))
  double precision,allocatable    :: work(:)
  integer                         :: info, lwork, i, j, k

  double precision,allocatable    :: A_tmp(:,:)
  allocate (A_tmp(LDA,n))
  A_tmp(:,:) = A(:,:)

  ! Find optimal size for temp arrays
  allocate(work(1))
  lwork = -1
  call dgesvd('A','A', m, n, A_tmp, LDA,                             &
      D, U, LDU, Vt, LDVt, work, lwork, info)
  ! /!\ int(WORK(1)) becomes negative when WORK(1) > 2147483648
  lwork = max(int(work(1)), 5*MIN(M,N))
  deallocate(work)

  allocate(work(lwork))
  call dgesvd('A','A', m, n, A_tmp, LDA,                             &
      D, U, LDU, Vt, LDVt, work, lwork, info)
  deallocate(A_tmp,work)

  if (info /= 0) then
    print *,  info, ': SVD failed'
    stop
  endif

  do j=1,m
    do i=1,m
      if (dabs(U(i,j)) < 1.d-14)  U(i,j) = 0.d0
    enddo
  enddo

  do j=1,n
    do i=1,n
      if (dabs(Vt(i,j)) < 1.d-14) Vt(i,j) = 0.d0
    enddo
  enddo

end

subroutine svd_symm(A,LDA,U,LDU,D,Vt,LDVt,m,n)
  implicit none
  BEGIN_DOC
  ! Compute A = U.D.Vt
  !
  ! LDx : leftmost dimension of x
  !
  ! Dimension of A is m x n
  !
  END_DOC

  integer, intent(in)             :: LDA, LDU, LDVt, m, n
  double precision, intent(in)    :: A(LDA,n)
  double precision, intent(out)   :: U(LDU,min(m,n))
  double precision,intent(out)    :: Vt(LDVt,n)
  double precision,intent(out)    :: D(min(m,n))
  double precision,allocatable    :: work(:)
  integer                         :: info, lwork, i, j, k

  double precision,allocatable    :: A_tmp(:,:)
  allocate (A_tmp(LDA,n))
  A_tmp(:,:) = A(:,:)

  ! Find optimal size for temp arrays
  allocate(work(1))
  lwork = -1
  call dgesvd('A','A', m, n, A_tmp, LDA,                             &
      D, U, LDU, Vt, LDVt, work, lwork, info)
  ! /!\ int(WORK(1)) becomes negative when WORK(1) > 2147483648
  lwork = max(int(work(1)), 5*MIN(M,N))
  deallocate(work)

  allocate(work(lwork))
  call dgesvd('A','A', m, n, A_tmp, LDA,                             &
      D, U, LDU, Vt, LDVt, work, lwork, info)
  deallocate(A_tmp,work)

  if (info /= 0) then
    print *,  info, ': SVD failed'
    stop
  endif

  ! Iterative refinement
  ! --------------------
  ! https://doi.org/10.1016/j.cam.2019.112512

  integer :: iter
  double precision,allocatable    :: R(:,:), S(:,:), T(:,:)
  double precision,allocatable    :: sigma(:), F(:,:), G(:,:)
  double precision :: alpha, beta, x, thresh

  allocate (R(m,m), S(n,n), T(m,n), sigma(m), F(m,m), G(n,n), A_tmp(m,n))
  sigma = 0.d0
  R = 0.d0
  S = 0.d0
  T = 0.d0
  F = 0.d0
  G = 0.d0

  thresh = 1.d-8
  call restore_symmetry(m,m,U,size(U,1),thresh)
  call restore_symmetry(n,n,Vt,size(Vt,1),thresh)

  do iter=1,4
    do k=1,n
      A_tmp(1:m,k) = D(k) * U(1:m,k)
    enddo

    call dgemm('N','N',m,n,n,1.d0,A_tmp,size(A_tmp,1), &
       Vt,size(Vt,1),0.d0,R,size(R,1))

    print *, maxval(dabs(R(1:m,1:n) - A(1:m,1:n)))


    call dgemm('T','N',m,m,m,-1.d0,U,size(U,1), &
      U,size(U,1),0.d0,R,size(R,1))
    do i=1,m
      R(i,i) = R(i,i) + 1.d0
    enddo

    call dgemm('N','T',n,n,n,-1.d0,Vt,size(Vt,1), &
      Vt,size(Vt,1),0.d0,S,size(S,1))
    do i=1,n
      S(i,i) = S(i,i) + 1.d0
    enddo


    call dgemm('T','N',m,n,m,1.d0,U,size(U,1), &
       A,size(A,1),0.d0,A_tmp,size(A_tmp,1))

    call dgemm('N','T',m,n,n,1.d0,A_tmp,size(A_tmp,1), &
       Vt,size(Vt,1),0.d0,T,size(T,1))


    do i=1,n
      sigma(i) = T(i,i)/(1.d0 - (R(i,i)+S(i,i))*0.5d0)
      F(i,i) = 0.5d0*R(i,i)
      G(i,i) = 0.5d0*S(i,i)
    enddo

    do j=1,n
      do i=1,n
        if (i == j) cycle
        alpha = T(i,j) + sigma(j) * R(i,j)
        beta  = T(i,j) + sigma(j) * S(i,j)
        x = 1.d0 / (sigma(j)*sigma(j) - sigma(i)*sigma(i))
        F(i,j) = (alpha * sigma(j) + beta * sigma(i)) * x
        G(i,j) = (alpha * sigma(i) + beta * sigma(j)) * x
      enddo
    enddo

    do i=1,n
      x = 1.d0/sigma(i)
      do j=n+1,m
        F(i,j) = -T(j,i) * x
      enddo
    enddo

    do i=n+1,m
      do j=1,n
        F(i,j) = R(i,j) - F(j,i)
      enddo
    enddo

    do i=n+1,m
      do j=n+1,m
        F(i,j) = R(i,j)*0.5d0
      enddo
    enddo

    D(1:min(n,m)) = sigma(1:min(n,m))
    call dgemm('N','N',m,m,m,1.d0,U,size(U,1),F,size(F,1), &
       0.d0, A_tmp, size(A_tmp,1))
    do j=1,m
      do i=1,m
        U(i,j) = U(i,j) + A_tmp(i,j)
      enddo
    enddo

    call dgemm('T','N',n,n,n,1.d0,G,size(G,1),Vt,size(Vt,1), &
       0.d0, A_tmp, size(A_tmp,1))
    do j=1,n
      do i=1,n
        Vt(i,j) = Vt(i,j) + A_tmp(i,j)
      enddo
    enddo

    thresh = 0.01d0 * thresh
    call restore_symmetry(m,m,U,size(U,1),thresh)
    call restore_symmetry(n,n,Vt,size(Vt,1),thresh)

  enddo

  deallocate(A_tmp,R,S,F,G,sigma)
end

subroutine eigSVD(A,LDA,U,LDU,D,Vt,LDVt,m,n)
  implicit none
  BEGIN_DOC
! Algorithm 3 of https://arxiv.org/pdf/1810.06860.pdf
!
! A(m,n) = U(m,n) D(n) Vt(n,n) with m>n
  END_DOC
  integer, intent(in)            :: LDA, LDU, LDVt, m, n
  double precision, intent(in)   :: A(LDA,n)
  double precision, intent(out)  :: U(LDU,n)
  double precision,intent(out)   :: Vt(LDVt,n)
  double precision,intent(out)   :: D(n)

  integer                        :: i,j,k

  if (m<n) then
    stop -1
    call svd(A,LDA,U,LDU,D,Vt,LDVt,m,n)
    return
  endif

  double precision, allocatable :: B(:,:), V(:,:)
  allocate(B(n,n))
  ! B = - At . A
  call dgemm('T','N',n,n,m,-1.d0,A,size(A,1),A,size(A,1),0.d0,B,size(B,1))

  ! V, D = eig(B)
  allocate(V(n,n))
  call lapack_diagd(D,V,B,n,n)
  deallocate(B)
  do j=1,n
   do i=1,n
     Vt(i,j) = V(j,i)
   enddo
  enddo

  ! S = sqrt(-D)
  ! U = A.V.S^-1
  ! U = A.(S^-1.vt)t

  do k=1,n
    if (D(k) >= 0.d0) then
      exit
    endif
    D(k) = dsqrt(-D(k))
    call dscal(n, 1.d0/D(k), V(1,k), 1)
  enddo
  D(k:n) = 0.d0
  k=k-1
  call dgemm('N','N',m,n,k,1.d0,A,size(A,1),V,size(V,1),0.d0,U,size(U,1))

end


subroutine randomized_svd(A,LDA,U,LDU,D,Vt,LDVt,m,n,q,r)
  implicit none
  include 'constants.include.F'
  BEGIN_DOC
! Randomized SVD: rank r, q power iterations
!
! 1. Sample column space of A with P: Z = A.P where P is random with r+p columns.
!
! 2. Power iterations : Z <- X . (Xt.Z)
!
! 3. Z = Q.R
!
! 4. Compute SVD on projected Qt.X = U' . S. Vt
!
! 5. U = Q U'
  END_DOC

  integer, intent(in)             :: LDA, LDU, LDVt, m, n, q, r
  double precision, intent(in)    :: A(LDA,n)
  double precision, intent(out)   :: U(LDU,r)
  double precision,intent(out)    :: Vt(LDVt,r)
  double precision,intent(out)    :: D(r)
  integer                         :: i, j, k

  double precision,allocatable    :: Z(:,:), P(:,:), Y(:,:), UY(:,:)
  double precision :: r1,r2
  allocate(P(n,r), Z(m,r))

  ! P is a normal random matrix (n,r)
  do k=1,r
    do i=1,n
      call random_number(r1)
      call random_number(r2)
      r1 = dsqrt(-2.d0*dlog(r1))
      r2 = dtwo_pi*r2
      P(i,k) = r1*dcos(r2)
    enddo
  enddo

  ! Z(m,r) = A(m,n).P(n,r)
  call dgemm('N','N',m,r,n,1.d0,A,size(A,1),P,size(P,1),0.d0,Z,size(Z,1))

  ! Power iterations
  do k=1,q
    ! P(n,r) = At(n,m).Z(m,r)
    call dgemm('T','N',n,r,m,1.d0,A,size(A,1),Z,size(Z,1),0.d0,P,size(P,1))
    ! Z(m,r) = A(m,n).P(n,r)
    call dgemm('N','N',m,r,n,1.d0,A,size(A,1),P,size(P,1),0.d0,Z,size(Z,1))
  enddo

  deallocate(P)

  ! QR factorization of Z
  call ortho_svd(Z,size(Z,1),m,r)

  allocate(Y(r,n), UY(r,r))
  ! Y(r,n) = Zt(r,m).A(m,n)
  call dgemm('T','N',r,n,m,1.d0,Z,size(Z,1),A,size(A,1),0.d0,Y,size(Y,1))

  ! SVD of Y
  call svd(Y,size(Y,1),UY,size(UY,1),D,Vt,size(Vt,1),r,n)
  deallocate(Y)

  ! U(m,r) = Z(m,r).UY(r,r)
  call dgemm('N','N',m,r,r,1.d0,Z,size(Z,1),UY,size(UY,1),0.d0,U,size(U,1))
  deallocate(UY,Z)

end

subroutine svd_complex(A,LDA,U,LDU,D,Vt,LDVt,m,n)
  implicit none
  BEGIN_DOC
  ! Compute A = U.D.Vt
  !
  ! LDx : leftmost dimension of x
  !
  ! Dimension of A is m x n
  ! A,U,Vt are complex*16
  ! D is double precision
  END_DOC

  integer, intent(in)             :: LDA, LDU, LDVt, m, n
  complex*16, intent(in)          :: A(LDA,n)
  complex*16, intent(out)         :: U(LDU,m)
  complex*16, intent(out)         :: Vt(LDVt,n)
  double precision,intent(out)    :: D(min(m,n))
  complex*16,allocatable          :: work(:)
  double precision,allocatable    :: rwork(:)
  integer                         :: info, lwork, i, j, k, lrwork

  complex*16,allocatable          :: A_tmp(:,:)
  allocate (A_tmp(LDA,n))
  A_tmp = A
  lrwork = 5*min(m,n)

  ! Find optimal size for temp arrays
  allocate(work(1),rwork(lrwork))
  lwork = -1
  call zgesvd('A','A', m, n, A_tmp, LDA,                             &
      D, U, LDU, Vt, LDVt, work, lwork, rwork, info)
  lwork = int(work(1))
  deallocate(work)

  allocate(work(lwork))
  call zgesvd('A','A', m, n, A_tmp, LDA,                             &
      D, U, LDU, Vt, LDVt, work, lwork, rwork, info)
  deallocate(work,rwork,A_tmp)

  if (info /= 0) then
    print *,  info, ': SVD failed'
    stop
  endif

end

subroutine ortho_canonical_complex(overlap,LDA,N,C,LDC,m,cutoff)
  implicit none
  BEGIN_DOC
  ! Compute C_new=C_old.U.s^-1/2 canonical orthogonalization.
  !
  ! overlap : overlap matrix
  !
  ! LDA : leftmost dimension of overlap array
  !
  ! N : Overlap matrix is NxN (array is (LDA,N) )
  !
  ! C : Coefficients of the vectors to orthogonalize. On exit,
  !     orthogonal vectors
  !
  ! LDC : leftmost dimension of C
  !
  ! m : Coefficients matrix is MxN, ( array is (LDC,N) )
  !
  END_DOC

  integer, intent(in)            :: lda, ldc, n
  integer, intent(out)           :: m
  complex*16, intent(in)   :: overlap(lda,n)
  double precision, intent(in)   :: cutoff
  complex*16, intent(inout) :: C(ldc,n)
  complex*16, allocatable  :: U(:,:)
  complex*16, allocatable  :: Vt(:,:)
  double precision, allocatable  :: D(:)
  complex*16, allocatable  :: S(:,:)
  !DIR$ ATTRIBUTES ALIGN : 64    :: U, Vt, D
  integer                        :: info, i, j
  double precision :: local_cutoff

  if (n < 2) then
    return
  endif

  allocate (U(ldc,n), Vt(lda,n), D(n), S(lda,n))

  call svd_complex(overlap,lda,U,ldc,D,Vt,lda,n,n)

  D(:) = dsqrt(D(:))
  local_cutoff = dsqrt(cutoff)*D(1) ! such that D(i)/D(1) > dsqrt(cutoff) is kept
  m=n
  do i=1,n
    if ( D(i) >= local_cutoff ) then
      D(i) = 1.d0/D(i)
    else
      m = i-1
      print *,  'Removed Linear dependencies below:', local_cutoff
      exit
    endif
  enddo
  do i=m+1,n
    D(i) = 0.d0
  enddo

  do j=1,n
    do i=1,n
      S(i,j) = U(i,j)*D(j)
    enddo
  enddo

  do j=1,n
    do i=1,n
      U(i,j) = C(i,j)
    enddo
  enddo

  call zgemm('N','N',n,n,n,(1.d0,0.d0),U,size(U,1),S,size(S,1),(0.d0,0.d0),C,size(C,1))
  deallocate (U, Vt, D, S)

end


subroutine ortho_qr_complex(A,LDA,m,n)
  implicit none
  BEGIN_DOC
  ! Orthogonalization using Q.R factorization
  !
  ! A : matrix to orthogonalize
  !
  ! LDA : leftmost dimension of A
  !
  ! n : Number of rows of A
  !
  ! m : Number of columns of A
  !
  END_DOC
  integer, intent(in)            :: m,n, LDA
  complex*16, intent(inout) :: A(LDA,n)

  integer                        :: lwork, info
  integer, allocatable           :: jpvt(:)
  complex*16, allocatable  :: tau(:), work(:)

  allocate (jpvt(n), tau(n), work(1))
  LWORK=-1
  call  zgeqrf( m, n, A, LDA, TAU, WORK, LWORK, INFO )
  LWORK=2*int(WORK(1))
  deallocate(WORK)
  allocate(WORK(LWORK))
  call  zgeqrf(m, n, A, LDA, TAU, WORK, LWORK, INFO )
  call  zungqr(m, n, n, A, LDA, tau, WORK, LWORK, INFO)
  deallocate(WORK,jpvt,tau)
end

subroutine ortho_qr_unblocked_complex(A,LDA,m,n)
  implicit none
  BEGIN_DOC
  ! Orthogonalization using Q.R factorization
  !
  ! A : matrix to orthogonalize
  !
  ! LDA : leftmost dimension of A
  !
  ! n : Number of rows of A
  !
  ! m : Number of columns of A
  !
  END_DOC
  integer, intent(in)            :: m,n, LDA
  double precision, intent(inout) :: A(LDA,n)

  integer                        :: info
  integer, allocatable           :: jpvt(:)
  double precision, allocatable  :: tau(:), work(:)

  print *, irp_here, ': TO DO'
  stop -1

!  allocate (jpvt(n), tau(n), work(n))
!  call  dgeqr2( m, n, A, LDA, TAU, WORK, INFO )
!  call dorg2r(m, n, n, A, LDA, tau, WORK, INFO)
!  deallocate(WORK,jpvt,tau)
end

subroutine ortho_lowdin_complex(overlap,LDA,N,C,LDC,m,cutoff)
  implicit none
  BEGIN_DOC
  ! Compute C_new=C_old.S^-1/2 orthogonalization.
  !
  ! overlap : overlap matrix
  !
  ! LDA : leftmost dimension of overlap array
  !
  ! N : Overlap matrix is NxN (array is (LDA,N) )
  !
  ! C : Coefficients of the vectors to orthogonalize. On exit,
  !     orthogonal vectors
  !
  ! LDC : leftmost dimension of C
  !
  ! M : Coefficients matrix is MxN, ( array is (LDC,N) )
  !
  END_DOC

  integer, intent(in)            :: LDA, ldc, n, m
  complex*16, intent(in)   :: overlap(lda,n)
  complex*16, intent(inout) :: C(ldc,n)
  complex*16, allocatable  :: U(:,:)
  complex*16, allocatable  :: Vt(:,:)
  double precision, allocatable  :: D(:)
  complex*16, allocatable  :: S(:,:)
  double precision, intent(in) :: cutoff
  double precision :: local_cutoff
  integer                        :: info, i, j, k, mm

  if (n < 2) then
    return
  endif

  allocate(U(ldc,n),Vt(lda,n),S(lda,n),D(n))

  call svd_complex(overlap,lda,U,ldc,D,Vt,lda,n,n)
  D(:) = dsqrt(D(:))
  local_cutoff = dsqrt(cutoff)*D(1) ! such that D(i)/D(1) > dsqrt(cutoff) is kept
  mm=n
  do i=1,n
    if ( D(i) >= local_cutoff) then
      D(i) = 1.d0/D(i)
    else
      mm = mm-1
      D(i) = 0.d0
    endif
    do j=1,n
      S(j,i) = (0.d0,0.d0)
    enddo
  enddo

  if (mm < n) then
    print *,  'Removed Linear dependencies below ', local_cutoff
  endif

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP SHARED(S,U,D,Vt,n,C,m,local_cutoff) &
  !$OMP PRIVATE(i,j,k)

  do k=1,n
    if (D(k) /= 0.d0) then
      !$OMP DO SCHEDULE(STATIC)
      do j=1,n
        do i=1,n
          S(i,j) = S(i,j) + U(i,k)*D(k)*Vt(k,j)
        enddo
      enddo
      !$OMP END DO NOWAIT
    endif
  enddo

  !$OMP BARRIER
  !$OMP DO
  do j=1,n
    do i=1,m
      U(i,j) = C(i,j)
    enddo
  enddo
  !$OMP END DO

  !$OMP END PARALLEL

  call zgemm('N','N',m,n,n,(1.d0,0.d0),U,size(U,1),S,size(S,1),(0.d0,0.d0),C,size(C,1))

  deallocate(U,Vt,S,D)
end

subroutine get_inverse_complex(A,LDA,m,C,LDC)
  implicit none
  BEGIN_DOC
  ! Returns the inverse of the square matrix A
  END_DOC
  integer, intent(in)            :: m, LDA, LDC
  complex*16, intent(in)   :: A(LDA,m)
  complex*16, intent(out)  :: C(LDC,m)

  integer                        :: info,lwork
  integer, allocatable           :: ipiv(:)
  complex*16,allocatable   :: work(:)
  allocate (ipiv(m), work(m*m))
  lwork = size(work)
  C(1:m,1:m) = A(1:m,1:m)
  call zgetrf(m,m,C,size(C,1),ipiv,info)
  if (info /= 0) then
    print *,  info
    stop 'error in inverse (zgetrf)'
  endif
  call zgetri(m,C,size(C,1),ipiv,work,lwork,info)
  if (info /= 0) then
    print *,  info
    stop 'error in inverse (zgetri)'
  endif
  deallocate(ipiv,work)
end


subroutine get_pseudo_inverse_complex(A,LDA,m,n,C,LDC,cutoff)
  implicit none
  BEGIN_DOC
  ! Find C = A^-1
  END_DOC
  integer, intent(in)            :: m,n, LDA, LDC
  complex*16, intent(in)   :: A(LDA,n)
  double precision, intent(in)   :: cutoff
  complex*16, intent(out)  :: C(LDC,m)

  double precision, allocatable  :: D(:), rwork(:)
  complex*16, allocatable  :: U(:,:), Vt(:,:), work(:), A_tmp(:,:)
  integer                        :: info, lwork
  integer                        :: i,j,k
  allocate (D(n),U(m,n),Vt(n,n),work(1),A_tmp(m,n),rwork(5*n))
  do j=1,n
    do i=1,m
      A_tmp(i,j) = A(i,j)
    enddo
  enddo
  lwork = -1
  call zgesvd('S','A', m, n, A_tmp, m,D,U,m,Vt,n,work,lwork,rwork,info)
  if (info /= 0) then
    print *,  info, ': SVD failed'
    stop
  endif
  lwork = int(real(work(1)))
  deallocate(work)
  allocate(work(lwork))
  call zgesvd('S','A', m, n, A_tmp, m,D,U,m,Vt,n,work,lwork,rwork,info)
  if (info /= 0) then
    print *,  info, ':: SVD failed'
    stop 1
  endif

  do i=1,n
    if (D(i) > cutoff*D(1)) then
      D(i) = 1.d0/D(i)
    else
      D(i) = 0.d0
    endif
  enddo

  C = (0.d0,0.d0)
  do i=1,m
    do j=1,n
      do k=1,n
        C(j,i) = C(j,i) + U(i,k) * D(k) * Vt(k,j)
      enddo
    enddo
  enddo

  deallocate(U,D,Vt,work,A_tmp,rwork)

end

subroutine lapack_diagd_diag_in_place_complex(eigvalues,eigvectors,nmax,n)
  implicit none
  BEGIN_DOC
  ! Diagonalize matrix H(complex)
  !
  ! H is untouched between input and ouptut
  !
  ! eigevalues(i) = ith lowest eigenvalue of the H matrix
  !
  ! eigvectors(i,j) = <i|psi_j> where i is the basis function and psi_j is the j th eigenvector
  !
  END_DOC
  integer, intent(in)            :: n,nmax
!  double precision, intent(out)  :: eigvectors(nmax,n)
  complex*16, intent(inout)        :: eigvectors(nmax,n)
  double precision, intent(out)  :: eigvalues(n)
!  double precision, intent(in)   :: H(nmax,n)
  complex*16,allocatable         :: work(:)
  integer         ,allocatable   :: iwork(:)
!  complex*16,allocatable         :: A(:,:)
  double precision, allocatable  :: rwork(:)
  integer                        :: lrwork, lwork, info, i,j,l,k, liwork

! print*,'Diagonalization by jacobi'
! print*,'n = ',n

  lwork = 2*n*n + 2*n
  lrwork = 2*n*n + 5*n+ 1
  liwork = 5*n + 3
  allocate (work(lwork),iwork(liwork),rwork(lrwork))

  lwork = -1
  liwork = -1
  lrwork = -1
  ! get optimal work size
  call ZHEEVD( 'V', 'U', n, eigvectors, nmax, eigvalues, work, lwork, &
    rwork, lrwork, iwork, liwork, info )
  if (info < 0) then
    print *, irp_here, ': ZHEEVD: the ',-info,'-th argument had an illegal value'
    stop 2
  endif
  lwork  = int( real(work(1)))
  liwork = iwork(1)
  lrwork = int(rwork(1))
  deallocate (work,iwork,rwork)

  allocate (work(lwork),iwork(liwork),rwork(lrwork))
  call ZHEEVD( 'V', 'U', n, eigvectors, nmax, eigvalues, work, lwork, &
    rwork, lrwork, iwork, liwork, info )
  deallocate(work,iwork,rwork)


  if (info < 0) then
    print *, irp_here, ': ZHEEVD: the ',-info,'-th argument had an illegal value'
    stop 2
  else if( info > 0  ) then
     write(*,*)'ZHEEVD Failed; calling ZHEEV'
     lwork = 2*n - 1
     lrwork = 3*n - 2
     allocate(work(lwork),rwork(lrwork))
     lwork = -1
     call ZHEEV('V','L',n,eigvectors,nmax,eigvalues,work,lwork,rwork,info)
     if (info < 0) then
       print *, irp_here, ': ZHEEV: the ',-info,'-th argument had an illegal value'
       stop 2
     endif
     lwork = int(work(1))
     deallocate(work)
     allocate(work(lwork))
     call ZHEEV('V','L',n,eigvectors,nmax,eigvalues,work,lwork,rwork,info)
     if (info /= 0 ) then
       write(*,*)'ZHEEV Failed'
       stop 1
     endif
    deallocate(work,rwork)
  end if

end

subroutine lapack_diagd_diag_complex(eigvalues,eigvectors,H,nmax,n)
  implicit none
  BEGIN_DOC
  ! Diagonalize matrix H(complex)
  !
  ! H is untouched between input and ouptut
  !
  ! eigevalues(i) = ith lowest eigenvalue of the H matrix
  !
  ! eigvectors(i,j) = <i|psi_j> where i is the basis function and psi_j is the j th eigenvector
  !
  END_DOC
  integer, intent(in)            :: n,nmax
!  double precision, intent(out)  :: eigvectors(nmax,n)
  complex*16, intent(out)        :: eigvectors(nmax,n)
  double precision, intent(out)  :: eigvalues(n)
!  double precision, intent(in)   :: H(nmax,n)
  complex*16, intent(in)         :: H(nmax,n)
  double precision, allocatable  :: eigenvalues(:)
  complex*16,allocatable         :: work(:)
  integer         ,allocatable   :: iwork(:)
  complex*16,allocatable         :: A(:,:)
  double precision, allocatable  :: rwork(:)
  integer                        :: lrwork, lwork, info, i,j,l,k, liwork

  allocate(A(nmax,n),eigenvalues(n))
! print*,'Diagonalization by jacobi'
! print*,'n = ',n

  A=H
  lwork = 2*n*n + 2*n
  lrwork = 2*n*n + 5*n+ 1
  liwork = 5*n + 3
  allocate (work(lwork),iwork(liwork),rwork(lrwork))

  lwork = -1
  liwork = -1
  lrwork = -1
  ! get optimal work size
  call ZHEEVD( 'V', 'U', n, A, nmax, eigenvalues, work, lwork, &
    rwork, lrwork, iwork, liwork, info )
  if (info < 0) then
    print *, irp_here, ': ZHEEVD: the ',-info,'-th argument had an illegal value'
    stop 2
  endif
  lwork  = int( real(work(1)))
  liwork = iwork(1)
  lrwork = int(rwork(1))
  deallocate (work,iwork,rwork)

  allocate (work(lwork),iwork(liwork),rwork(lrwork))
  call ZHEEVD( 'V', 'U', n, A, nmax, eigenvalues, work, lwork, &
    rwork, lrwork, iwork, liwork, info )
  deallocate(work,iwork,rwork)

  if (info < 0) then
    print *, irp_here, ': ZHEEVD: the ',-info,'-th argument had an illegal value'
    stop 2
  else if( info > 0  ) then
     write(*,*)'ZHEEVD Failed; calling ZHEEV'
     lwork = 2*n - 1
     lrwork = 3*n - 2
     allocate(work(lwork),rwork(lrwork))
     lwork = -1
     call ZHEEV('V','L',n,A,nmax,eigenvalues,work,lwork,rwork,info)
     if (info < 0) then
       print *, irp_here, ': ZHEEV: the ',-info,'-th argument had an illegal value'
       stop 2
     endif
     lwork = int(work(1))
     deallocate(work)
     allocate(work(lwork))
     call ZHEEV('V','L',n,A,nmax,eigenvalues,work,lwork,rwork,info)
     if (info /= 0 ) then
       write(*,*)'ZHEEV Failed'
       stop 1
     endif
    deallocate(work,rwork)
  end if

  eigvectors = (0.d0,0.d0)
  eigvalues = 0.d0
  do j = 1, n
    eigvalues(j) = eigenvalues(j)
    do i = 1, n
      eigvectors(i,j) = A(i,j)
    enddo
  enddo
  deallocate(A,eigenvalues)
end

subroutine lapack_diagd_complex(eigvalues,eigvectors,H,nmax,n)
  implicit none
  BEGIN_DOC
  ! Diagonalize matrix H(complex)
  !
  ! H is untouched between input and ouptut
  !
  ! eigevalues(i) = ith lowest eigenvalue of the H matrix
  !
  ! eigvectors(i,j) = <i|psi_j> where i is the basis function and psi_j is the j th eigenvector
  !
  END_DOC
  integer, intent(in)            :: n,nmax
!  double precision, intent(out)  :: eigvectors(nmax,n)
  complex*16, intent(out)        :: eigvectors(nmax,n)
  double precision, intent(out)  :: eigvalues(n)
!  double precision, intent(in)   :: H(nmax,n)
  complex*16, intent(in)         :: H(nmax,n)
  double precision, allocatable  :: eigenvalues(:)
  complex*16,allocatable         :: work(:)
  integer         ,allocatable   :: iwork(:)
  complex*16,allocatable         :: A(:,:)
  double precision, allocatable  :: rwork(:)
  integer                        :: lrwork, lwork, info, i,j,l,k, liwork

  allocate(A(nmax,n),eigenvalues(n))
! print*,'Diagonalization by jacobi'
! print*,'n = ',n

  A=H
  lwork = 2*n*n + 2*n
  lrwork = 2*n*n + 5*n+ 1
  liwork = 5*n + 3
  allocate (work(lwork),iwork(liwork),rwork(lrwork))

  lwork = -1
  liwork = -1
  lrwork = -1
  call ZHEEVD( 'V', 'U', n, A, nmax, eigenvalues, work, lwork, &
    rwork, lrwork, iwork, liwork, info )
  if (info < 0) then
    print *, irp_here, ': ZHEEVD: the ',-info,'-th argument had an illegal value'
    stop 2
  endif
  lwork  = max(int( work( 1 ) ),lwork)
  liwork = iwork(1)
  lrwork = max(int(rwork(1),4),lrwork)
  deallocate (work,iwork,rwork)

  allocate (work(lwork),iwork(liwork),rwork(lrwork))
  call ZHEEVD( 'V', 'U', n, A, nmax, eigenvalues, work, lwork, &
    rwork, lrwork, iwork, liwork, info )
  deallocate(work,iwork,rwork)


  if (info < 0) then
    print *, irp_here, ': ZHEEVD: the ',-info,'-th argument had an illegal value'
    stop 2
  else if( info > 0  ) then
     write(*,*)'ZHEEVD Failed'
     stop 1
  end if

  eigvectors = (0.d0,0.d0)
  eigvalues = 0.d0
  do j = 1, n
    eigvalues(j) = eigenvalues(j)
    do i = 1, n
      eigvectors(i,j) = A(i,j)
    enddo
  enddo
  deallocate(A,eigenvalues)
end

subroutine lapack_diag_complex(eigvalues,eigvectors,H,nmax,n)
  implicit none
  BEGIN_DOC
  ! Diagonalize matrix H (complex)
  !
  ! H is untouched between input and ouptut
  !
  ! eigevalues(i) = ith lowest eigenvalue of the H matrix
  !
  ! eigvectors(i,j) = <i|psi_j> where i is the basis function and psi_j is the j th eigenvector
  !
  END_DOC
  integer, intent(in)            :: n,nmax
  complex*16, intent(out)        :: eigvectors(nmax,n)
  double precision, intent(out)  :: eigvalues(n)
  complex*16, intent(in)   :: H(nmax,n)
  double precision,allocatable   :: eigenvalues(:)
  complex*16,allocatable   :: work(:)
  complex*16,allocatable   :: A(:,:)
  double precision,allocatable   :: rwork(:)
  integer                        :: lwork, info, i,j,l,k,lrwork

  allocate(A(nmax,n),eigenvalues(n))
! print*,'Diagonalization by jacobi'
! print*,'n = ',n

  A=H
  !lwork = 2*n*n + 6*n+ 1
  lwork = 2*n - 1
  lrwork = 3*n - 2
  allocate (work(lwork),rwork(lrwork))

  lwork = -1
  call ZHEEV( 'V', 'U', n, A, nmax, eigenvalues, work, lwork, &
    rwork, info )
  if (info < 0) then
    print *, irp_here, ': ZHEEV: the ',-info,'-th argument had an illegal value'
    stop 2
  endif
  lwork  = int( work( 1 ) )
  deallocate (work)

  allocate (work(lwork))
  call ZHEEV( 'V', 'U', n, A, nmax, eigenvalues, work, lwork, &
    rwork, info )
  deallocate(work,rwork)

  if (info < 0) then
    print *, irp_here, ': ZHEEV: the ',-info,'-th argument had an illegal value'
    stop 2
  else if( info > 0  ) then
     write(*,*)'ZHEEV Failed : ', info
     do i=1,n
      do j=1,n
        print *,  H(i,j)
      enddo
     enddo
     stop 1
  end if

  eigvectors = (0.d0,0.d0)
  eigvalues = 0.d0
  do j = 1, n
    eigvalues(j) = eigenvalues(j)
    do i = 1, n
      eigvectors(i,j) = A(i,j)
    enddo
  enddo
  deallocate(A,eigenvalues)
end

subroutine matrix_vector_product_complex(u0,u1,matrix,sze,lda)
 implicit none
 BEGIN_DOC
! performs u1 += u0 * matrix
 END_DOC
 integer, intent(in)             :: sze,lda
 complex*16, intent(in)    :: u0(sze)
 complex*16, intent(inout) :: u1(sze)
 complex*16, intent(in)    :: matrix(lda,sze)
 integer :: i,j
 integer                        :: incx,incy
 incx = 1
 incy = 1
 !call dsymv('U', sze, 1.d0, matrix, lda, u0, incx, 1.d0, u1, incy)
 call zhemv('U', sze, (1.d0,0.d0), matrix, lda, u0, incx, (1.d0,0.d0), u1, incy)
end

subroutine ortho_canonical(overlap,LDA,N,C,LDC,m,cutoff)
  implicit none
  BEGIN_DOC
  ! Compute C_new=C_old.U.s^-1/2 canonical orthogonalization.
  !
  ! overlap : overlap matrix
  !
  ! LDA : leftmost dimension of overlap array
  !
  ! N : Overlap matrix is NxN (array is (LDA,N) )
  !
  ! C : Coefficients of the vectors to orthogonalize. On exit,
  !     orthogonal vectors
  !
  ! LDC : leftmost dimension of C
  !
  ! m : Coefficients matrix is MxN, ( array is (LDC,N) )
  !
  END_DOC

  integer, intent(in)            :: lda, ldc, n
  integer, intent(out)           :: m
  double precision, intent(in)   :: overlap(lda,n)
  double precision, intent(in)   :: cutoff
  double precision, intent(inout) :: C(ldc,n)
  double precision, allocatable  :: U(:,:)
  double precision, allocatable  :: Vt(:,:)
  double precision, allocatable  :: D(:)
  double precision, allocatable  :: S(:,:)
  !DIR$ ATTRIBUTES ALIGN : 64    :: U, Vt, D
  integer                        :: info, i, j
  double precision :: local_cutoff

  if (n < 2) then
    return
  endif

  allocate (U(ldc,n), Vt(lda,n), D(n), S(lda,n))

  call svd(overlap,lda,U,ldc,D,Vt,lda,n,n)

  D(:) = dsqrt(D(:))
  local_cutoff = dsqrt(cutoff)*D(1) ! such that D(i)/D(1) > dsqrt(cutoff) is kept
  m=n
  do i=1,n
    if ( D(i) >= local_cutoff ) then
      D(i) = 1.d0/D(i)
    else
      m = i-1
      print *,  'Removed Linear dependencies below:', local_cutoff
      exit
    endif
  enddo
  do i=m+1,n
    D(i) = 0.d0
  enddo

  do j=1,n
    do i=1,n
      S(i,j) = U(i,j)*D(j)
    enddo
  enddo

  do j=1,n
    do i=1,n
      U(i,j) = C(i,j)
    enddo
  enddo

  call dgemm('N','N',n,n,n,1.d0,U,size(U,1),S,size(S,1),0.d0,C,size(C,1))
  deallocate (U, Vt, D, S)

end


subroutine ortho_svd(A,LDA,m,n)
  implicit none
  BEGIN_DOC
  ! Orthogonalization via fast SVD
  !
  ! A : matrix to orthogonalize
  !
  ! LDA : leftmost dimension of A
  !
  ! m : Number of rows of A
  !
  ! n : Number of columns of A
  !
  END_DOC
  integer, intent(in)            :: m,n, LDA
  double precision, intent(inout) :: A(LDA,n)
  if (m < n) then
    call ortho_qr(A,LDA,m,n)
  endif
  double precision, allocatable :: U(:,:), D(:), Vt(:,:)
  allocate(U(m,n), D(n), Vt(n,n))
  call SVD(A,LDA,U,size(U,1),D,Vt,size(Vt,1),m,n)
  A(1:m,1:n) = U(1:m,1:n)
  deallocate(U,D, Vt)

end

subroutine ortho_qr(A,LDA,m,n)
  implicit none
  BEGIN_DOC
  ! Orthogonalization using Q.R factorization
  !
  ! A : matrix to orthogonalize
  !
  ! LDA : leftmost dimension of A
  !
  ! m : Number of rows of A
  !
  ! n : Number of columns of A
  !
  END_DOC
  integer, intent(in)            :: m,n, LDA
  double precision, intent(inout) :: A(LDA,n)

  integer                        :: LWORK, INFO
  double precision, allocatable  :: TAU(:), WORK(:)

  allocate (TAU(min(m,n)), WORK(1))

  LWORK=-1
  call dgeqrf( m, n, A, LDA, TAU, WORK, LWORK, INFO )
  ! /!\ int(WORK(1)) becomes negative when WORK(1) > 2147483648
  LWORK=max(n,int(WORK(1)))

  deallocate(WORK)
  allocate(WORK(LWORK))
  call dgeqrf(m, n, A, LDA, TAU, WORK, LWORK, INFO )

  LWORK=-1
  call dorgqr(m, n, n, A, LDA, TAU, WORK, LWORK, INFO)
  ! /!\ int(WORK(1)) becomes negative when WORK(1) > 2147483648
  LWORK=max(n,int(WORK(1)))

  deallocate(WORK)
  allocate(WORK(LWORK))
  call dorgqr(m, n, n, A, LDA, TAU, WORK, LWORK, INFO)

  deallocate(WORK,TAU)
end

subroutine ortho_qr_unblocked(A,LDA,m,n)
  implicit none
  BEGIN_DOC
  ! Orthogonalization using Q.R factorization
  !
  ! A : matrix to orthogonalize
  !
  ! LDA : leftmost dimension of A
  !
  ! n : Number of rows of A
  !
  ! m : Number of columns of A
  !
  END_DOC
  integer, intent(in)            :: m,n, LDA
  double precision, intent(inout) :: A(LDA,n)

  integer                        :: info
  double precision, allocatable  :: TAU(:), WORK(:)

  allocate (TAU(n), WORK(n))
  call dgeqr2( m, n, A, LDA, TAU, WORK, INFO )
  call dorg2r(m, n, n, A, LDA, TAU, WORK, INFO)
  deallocate(WORK,TAU)
end

subroutine ortho_lowdin(overlap,LDA,N,C,LDC,m,cutoff)
  implicit none
  BEGIN_DOC
  ! Compute C_new=C_old.S^-1/2 orthogonalization.
  !
  ! overlap : overlap matrix
  !
  ! LDA : leftmost dimension of overlap array
  !
  ! N : Overlap matrix is NxN (array is (LDA,N) )
  !
  ! C : Coefficients of the vectors to orthogonalize. On exit,
  !     orthogonal vectors
  !
  ! LDC : leftmost dimension of C
  !
  ! M : Coefficients matrix is MxN, ( array is (LDC,N) )
  !
  END_DOC

  integer, intent(in)            :: LDA, ldc, n, m
  double precision, intent(in)   :: overlap(lda,n)
  double precision, intent(in)   :: cutoff
  double precision, intent(inout) :: C(ldc,n)
  double precision, allocatable  :: U(:,:)
  double precision, allocatable  :: Vt(:,:)
  double precision, allocatable  :: D(:)
  double precision, allocatable  :: S(:,:)
  integer                        :: info, i, j, k, mm
  double precision :: local_cutoff

  if (n < 2) then
    return
  endif

  allocate(U(ldc,n),Vt(lda,n),S(lda,n),D(n))

  call svd(overlap,lda,U,ldc,D,Vt,lda,n,n)
  D(:) = dsqrt(D(:))
  local_cutoff = dsqrt(cutoff)*D(1) ! such that D(i)/D(1) > dsqrt(cutoff) is kept
  mm=n
  do i=1,n
    if ( D(i) >= local_cutoff) then
      D(i) = 1.d0/D(i)
    else
      mm = mm-1
      D(i) = 0.d0
    endif
    do j=1,n
      S(j,i) = 0.d0
    enddo
  enddo

  if (mm < n) then
    print *,  'Removed Linear dependencies below ', local_cutoff
  endif

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP SHARED(S,U,D,Vt,n,C,m,cutoff) &
  !$OMP PRIVATE(i,j,k)

  do k=1,n
    if (D(k) /= 0.d0) then
      !$OMP DO
      do j=1,n
        do i=1,n
          S(i,j) = S(i,j) + U(i,k)*D(k)*Vt(k,j)
        enddo
      enddo
      !$OMP END DO NOWAIT
    endif
  enddo

  !$OMP BARRIER
  !$OMP DO
  do j=1,n
    do i=1,m
      U(i,j) = C(i,j)
    enddo
  enddo
  !$OMP END DO

  !$OMP END PARALLEL

  call dgemm('N','N',m,n,n,1.d0,U,size(U,1),S,size(S,1),0.d0,C,size(C,1))

  deallocate(U,Vt,S,D)
end



subroutine get_inverse(A,LDA,m,C,LDC)
  implicit none
  BEGIN_DOC
  ! Returns the inverse of the square matrix A
  END_DOC
  integer, intent(in)            :: m, LDA, LDC
  double precision, intent(in)   :: A(LDA,m)
  double precision, intent(out)  :: C(LDC,m)

  integer                        :: info,lwork
  integer, allocatable           :: ipiv(:)
  double precision,allocatable   :: work(:)
  allocate (ipiv(m), work(m*m))
  lwork = size(work)
  C(1:m,1:m) = A(1:m,1:m)
  call dgetrf(m,m,C,size(C,1),ipiv,info)
  if (info /= 0) then
    print *,  info
    stop 'error in inverse (dgetrf)'
  endif
  call dgetri(m,C,size(C,1),ipiv,work,lwork,info)
  if (info /= 0) then
    print *,  info
    stop 'error in inverse (dgetri)'
  endif
  deallocate(ipiv,work)
end

subroutine get_pseudo_inverse(A,LDA,m,n,C,LDC,cutoff)
  implicit none
  BEGIN_DOC
  ! Find C = A^-1
  END_DOC
  integer, intent(in)            :: m,n, LDA, LDC
  double precision, intent(in)   :: A(LDA,n)
  double precision, intent(in)   :: cutoff
  double precision, intent(out)  :: C(LDC,m)

  double precision, allocatable  :: U(:,:), D(:), Vt(:,:), work(:), A_tmp(:,:)
  integer                        :: info, lwork
  integer                        :: i,j,k
  allocate (D(n),U(m,n),Vt(n,n),work(1),A_tmp(m,n))
  do j=1,n
    do i=1,m
      A_tmp(i,j) = A(i,j)
    enddo
  enddo
  lwork = -1
  call dgesvd('S','A', m, n, A_tmp, m,D,U,m,Vt,n,work,lwork,info)
  if (info /= 0) then
    print *,  info, ': SVD failed'
    stop
  endif
  LWORK=max(5*min(m,n),int(WORK(1)))
  deallocate(work)
  allocate(work(lwork))
  call dgesvd('S','A', m, n, A_tmp, m,D,U,m,Vt,n,work,lwork,info)
  if (info /= 0) then
    print *,  info, ':: SVD failed'
    stop 1
  endif

  do i=1,n
    if (D(i)/D(1) > cutoff) then
      D(i) = 1.d0/D(i)
    else
      D(i) = 0.d0
    endif
  enddo

  C = 0.d0
  do i=1,m
    do j=1,n
      do k=1,n
        C(j,i) = C(j,i) + U(i,k) * D(k) * Vt(k,j)
      enddo
    enddo
  enddo

  deallocate(U,D,Vt,work,A_tmp)

end



subroutine find_rotation(A,LDA,B,m,C,n)
  implicit none
  BEGIN_DOC
  ! Find A.C = B
  END_DOC
  integer, intent(in)            :: m,n, LDA
  double precision, intent(in)   :: A(LDA,n), B(LDA,n)
  double precision, intent(out)  :: C(n,n)

  double precision, allocatable  :: A_inv(:,:)
  allocate(A_inv(LDA,n))
  call get_pseudo_inverse(A,LDA,m,n,A_inv,LDA,1.d-10)

  integer                        :: i,j,k
  call dgemm('N','N',n,n,m,1.d0,A_inv,n,B,LDA,0.d0,C,n)
  deallocate(A_inv)
end


subroutine apply_rotation(A,LDA,R,LDR,B,LDB,m,n)
  implicit none
  BEGIN_DOC
  ! Apply the rotation found by find_rotation
  END_DOC
  integer, intent(in)            :: m,n, LDA, LDB, LDR
  double precision, intent(in)   :: R(LDR,n)
  double precision, intent(in)   :: A(LDA,n)
  double precision, intent(out)  :: B(LDB,n)
  call dgemm('N','N',m,n,n,1.d0,A,LDA,R,LDR,0.d0,B,LDB)
end

subroutine lapack_diagd(eigvalues,eigvectors,H,nmax,n)
  implicit none
  BEGIN_DOC
  ! Diagonalize matrix H
  !
  ! H is untouched between input and ouptut
  !
  ! eigevalues(i) = ith lowest eigenvalue of the H matrix
  !
  ! eigvectors(i,j) = <i|psi_j> where i is the basis function and psi_j is the j th eigenvector
  !
  END_DOC
  integer, intent(in)            :: n,nmax
  double precision, intent(out)  :: eigvectors(nmax,n)
  double precision, intent(out)  :: eigvalues(n)
  double precision, intent(in)   :: H(nmax,n)
  double precision,allocatable   :: eigenvalues(:)
  double precision,allocatable   :: work(:)
  integer         ,allocatable   :: iwork(:)
  double precision,allocatable   :: A(:,:)
  integer                        :: lwork, info, i,j,l,k, liwork

  allocate(A(nmax,n),eigenvalues(n))
  ! print*,'Diagonalization by jacobi'
  ! print*,'n = ',n

  A=H
  lwork = 1
  liwork = 1
  allocate (work(lwork),iwork(liwork))

  lwork = -1
  liwork = -1
  call DSYEVD( 'V', 'U', n, A, nmax, eigenvalues, work, lwork,       &
      iwork, liwork, info )
  if (info < 0) then
    print *, irp_here, ': DSYEVD: the ',-info,'-th argument had an illegal value'
    stop 2
  endif
  ! /!\ int(WORK(1)) becomes negative when WORK(1) > 2147483648
  LWORK  = max(int(work(1)), 2*n*n + 6*n+ 1)
  liwork = max(iwork(1), 5*n + 3)
  deallocate (work,iwork)

  allocate (work(lwork),iwork(liwork))
  call DSYEVD( 'V', 'U', n, A, nmax, eigenvalues, work, lwork,       &
      iwork, liwork, info )
  deallocate(work,iwork)

  if (info < 0) then
    print *, irp_here, ': DSYEVD: the ',-info,'-th argument had an illegal value'
    stop 2
  else if( info > 0  ) then
    write(*,*)'DSYEVD Failed'
    stop 1
  end if

  eigvectors = 0.d0
  eigvalues = 0.d0
  do j = 1, n
    eigvalues(j) = eigenvalues(j)
    do i = 1, n
      eigvectors(i,j) = A(i,j)
    enddo
  enddo
  deallocate(A,eigenvalues)
end

subroutine lapack_diag(eigvalues,eigvectors,H,nmax,n)
  implicit none
  BEGIN_DOC
  ! Diagonalize matrix H
  !
  ! H is untouched between input and ouptut
  !
  ! eigevalues(i) = ith lowest eigenvalue of the H matrix
  !
  ! eigvectors(i,j) = <i|psi_j> where i is the basis function and psi_j is the j th eigenvector
  !
  END_DOC
  integer, intent(in)            :: n,nmax
  double precision, intent(out)  :: eigvectors(nmax,n)
  double precision, intent(out)  :: eigvalues(n)
  double precision, intent(in)   :: H(nmax,n)
  double precision,allocatable   :: eigenvalues(:)
  double precision,allocatable   :: work(:)
  double precision,allocatable   :: A(:,:)
  integer                        :: lwork, info, i,j,l,k, liwork

  allocate(A(nmax,n),eigenvalues(n))

  A=H
  lwork = 1
  allocate (work(lwork))

  lwork = -1
  call DSYEV( 'V', 'U', n, A, nmax, eigenvalues, work, lwork, &
    info )
  if (info < 0) then
    print *, irp_here, ': DSYEV: the ',-info,'-th argument had an illegal value'
    stop 2
  endif
  ! /!\ int(WORK(1)) becomes negative when WORK(1) > 2147483648
  LWORK  = max(int(work(1)), 2*n*n + 6*n+ 1)
  deallocate (work)

  allocate (work(lwork))
  call DSYEV( 'V', 'U', n, A, nmax, eigenvalues, work, lwork, &
    info )
  deallocate(work)

  if (info < 0) then
    print *, irp_here, ': DSYEV: the ',-info,'-th argument had an illegal value'
    stop 2
  else if( info > 0  ) then
     write(*,*)'DSYEV Failed : ', info
     do i=1,n
      do j=1,n
        print *,  H(i,j)
      enddo
     enddo
     stop 1
  end if

  eigvectors = 0.d0
  eigvalues = 0.d0
  do j = 1, n
    eigvalues(j) = eigenvalues(j)
    do i = 1, n
      eigvectors(i,j) = A(i,j)
    enddo
  enddo
  deallocate(A,eigenvalues)
end

subroutine nullify_small_elements(m,n,A,LDA,thresh)
  implicit none
  integer, intent(in) :: m,n,LDA
  double precision, intent(inout) :: A(LDA,n)
  double precision, intent(in) :: thresh
  integer :: i,j
  double precision :: amax

  ! Find max value
  amax = 0.d0
  do j=1,n
    do i=1,m
      amax = max(dabs(A(i,j)), amax)
    enddo
  enddo
  if (amax == 0.d0) return
  amax = 1.d0/amax

  ! Remove tiny elements
  do j=1,n
    do i=1,m
      if ( dabs(A(i,j) * amax) < thresh ) then
         A(i,j) = 0.d0
      endif
    enddo
  enddo

end

subroutine restore_symmetry(m,n,A,LDA,thresh)
  implicit none
  BEGIN_DOC
! Tries to find the matrix elements that are the same, and sets them
! to the average value.
! If restore_symm is False, only nullify small elements
  END_DOC
  integer, intent(in) :: m,n,LDA
  double precision, intent(inout) :: A(LDA,n)
  double precision, intent(in) :: thresh
  integer :: i,j,k,l
  logical, allocatable :: done(:,:)
  double precision :: f, g, count, thresh2
  thresh2 = dsqrt(thresh)
  call nullify_small_elements(m,n,A,LDA,thresh)

  if (.not.restore_symm) then
    return
  endif

  ! TODO:  Costs O(n^4), but can be improved to (2 n^2 * log(n)):
  ! - copy all values in a 1D array
  ! - sort 1D array
  ! - average nearby elements 
  ! - for all elements, find matching value in the sorted 1D array

  allocate(done(m,n))

  do j=1,n
    do i=1,m
      done(i,j) = A(i,j) == 0.d0
    enddo
  enddo

  do j=1,n
    do i=1,m
      if ( done(i,j) ) cycle
      done(i,j) = .True.
      count = 1.d0
      f = 1.d0/A(i,j)
      do l=1,n
        do k=1,m
          if ( done(k,l) ) cycle
          g = f * A(k,l)
          if ( dabs(dabs(g) - 1.d0) < thresh2 ) then
            count = count + 1.d0
            if (g>0.d0) then
              A(i,j) = A(i,j) + A(k,l)
            else
              A(i,j) = A(i,j) - A(k,l)
            end if
          endif
        enddo
      enddo
      if (count > 1.d0) then
        A(i,j) = A(i,j) / count
        do l=1,n
          do k=1,m
            if ( done(k,l) ) cycle
            g = f * A(k,l)
            if ( dabs(dabs(g) - 1.d0) < thresh2 ) then
              done(k,l) = .True.
              if (g>0.d0) then
                A(k,l) = A(i,j)
              else
                A(k,l) = -A(i,j)
              end if
            endif
          enddo
        enddo
      endif

    enddo
  enddo

end
