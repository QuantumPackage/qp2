subroutine svd(A,LDA,U,LDU,D,Vt,LDVt,m,n)
  implicit none
  BEGIN_DOC
  ! Compute A = U.D.Vt
  !
  ! LDx : leftmost dimension of x
  !
  ! Dimsneion of A is m x n
  !
  END_DOC

  integer, intent(in)             :: LDA, LDU, LDVt, m, n
  double precision, intent(in)    :: A(LDA,n)
  double precision, intent(out)   :: U(LDU,m)
  double precision,intent(out)    :: Vt(LDVt,n)
  double precision,intent(out)    :: D(min(m,n))
  double precision,allocatable    :: work(:)
  integer                         :: info, lwork, i, j, k

  double precision,allocatable    :: A_tmp(:,:)
  allocate (A_tmp(LDA,n))
  A_tmp = A

  ! Find optimal size for temp arrays
  allocate(work(1))
  lwork = -1
  call dgesvd('A','A', m, n, A_tmp, LDA,                             &
      D, U, LDU, Vt, LDVt, work, lwork, info)
  lwork = int(work(1))
  deallocate(work)

  allocate(work(lwork))
  call dgesvd('A','A', m, n, A_tmp, LDA,                             &
      D, U, LDU, Vt, LDVt, work, lwork, info)
  deallocate(work,A_tmp)

  if (info /= 0) then
    print *,  info, ': SVD failed'
    stop
  endif

end


subroutine ortho_canonical(overlap,LDA,N,C,LDC,m)
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
  double precision, intent(inout) :: C(ldc,n)
  double precision, allocatable  :: U(:,:)
  double precision, allocatable  :: Vt(:,:)
  double precision, allocatable  :: D(:)
  double precision, allocatable  :: S(:,:)
  !DIR$ ATTRIBUTES ALIGN : 64    :: U, Vt, D
  integer                        :: info, i, j

  if (n < 2) then
    return
  endif

  allocate (U(ldc,n), Vt(lda,n), D(n), S(lda,n))

  call svd(overlap,lda,U,ldc,D,Vt,lda,n,n)

  D(:) = dsqrt(D(:))
  m=n
  do i=1,n
    if ( D(i) >= 1.d-6 ) then
      D(i) = 1.d0/D(i)
    else
      m = i-1
      print *,  'Removed Linear dependencies below:', 1.d0/D(m)
      exit
    endif
  enddo
  do i=m+1,n
    D(i) = 0.d0
  enddo

  do i=1,m
    if ( D(i) >= 1.d5 ) then
      print *,  'Warning: Basis set may have linear dependence problems'
    endif
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


subroutine ortho_qr(A,LDA,m,n)
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

  integer                        :: lwork, info
  double precision, allocatable  :: TAU(:), WORK(:)

  allocate (TAU(n), WORK(1))

  LWORK=-1
  call dgeqrf( m, n, A, LDA, TAU, WORK, LWORK, INFO )
  LWORK=int(WORK(1))

  deallocate(WORK)
  allocate(WORK(LWORK))
  call dgeqrf(m, n, A, LDA, TAU, WORK, LWORK, INFO )

  LWORK=-1
  call dorgqr(m, n, n, A, LDA, TAU, WORK, LWORK, INFO)
  LWORK=int(WORK(1))

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

subroutine ortho_lowdin(overlap,LDA,N,C,LDC,m)
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
  double precision, intent(inout) :: C(ldc,n)
  double precision, allocatable  :: U(:,:)
  double precision, allocatable  :: Vt(:,:)
  double precision, allocatable  :: D(:)
  double precision, allocatable  :: S(:,:)
  integer                        :: info, i, j, k

  if (n < 2) then
    return
  endif

  allocate(U(ldc,n),Vt(lda,n),S(lda,n),D(n))

  call svd(overlap,lda,U,ldc,D,Vt,lda,n,n)

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP SHARED(S,U,D,Vt,n,C,m) &
  !$OMP PRIVATE(i,j,k)

  !$OMP DO
  do i=1,n
    if ( D(i) < 1.d-6 ) then
      D(i) = 0.d0
    else
      D(i) = 1.d0/dsqrt(D(i))
    endif
    do j=1,n
      S(j,i) = 0.d0
    enddo
  enddo
  !$OMP END DO

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

subroutine get_pseudo_inverse(A,LDA,m,n,C,LDC)
  implicit none
  BEGIN_DOC
  ! Find C = A^-1
  END_DOC
  integer, intent(in)            :: m,n, LDA, LDC
  double precision, intent(in)   :: A(LDA,n)
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
  lwork = int(work(1))
  deallocate(work)
  allocate(work(lwork))
  call dgesvd('S','A', m, n, A_tmp, m,D,U,m,Vt,n,work,lwork,info)
  if (info /= 0) then
    print *,  info, ':: SVD failed'
    stop 1
  endif

  do i=1,n
    if (D(i)/D(1) > 1.d-10) then
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
  call get_pseudo_inverse(A,LDA,m,n,A_inv,LDA)

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
  lwork = max(1000,2*n*n + 6*n+ 1)
  liwork = max(5*n + 3,1000)
  allocate (work(lwork),iwork(liwork))

  lwork = -1
  liwork = -1
  call DSYEVD( 'V', 'U', n, A, nmax, eigenvalues, work, lwork,       &
      iwork, liwork, info )
  if (info < 0) then
    print *, irp_here, ': DSYEVD: the ',-info,'-th argument had an illegal value'
    stop 2
  endif
  lwork  = int( work( 1 ) )
  liwork = iwork(1)
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
  lwork = 2*n*n + 6*n+ 1
  allocate (work(lwork))

  lwork = -1
  call DSYEV( 'V', 'U', n, A, nmax, eigenvalues, work, lwork, &
    info )
  if (info < 0) then
    print *, irp_here, ': DSYEV: the ',-info,'-th argument had an illegal value'
    stop 2
  endif
  lwork  = int( work( 1 ) )
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

