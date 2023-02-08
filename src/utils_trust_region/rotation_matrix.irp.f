! Rotation matrix

! *Build a rotation matrix from an antisymmetric matrix*

! Compute a rotation matrix $\textbf{R}$ from an antisymmetric matrix $$\textbf{A}$$ such as :
! $$
! \textbf{R}=\exp(\textbf{A})
! $$

! So :
! \begin{align*}
! \textbf{R}=& \exp(\textbf{A}) \\
! =& \sum_k^{\infty} \frac{1}{k!}\textbf{A}^k \\
! =& \textbf{W} \cdot \cos(\tau) \cdot \textbf{W}^{\dagger} + \textbf{W} \cdot \tau^{-1} \cdot \sin(\tau) \cdot \textbf{W}^{\dagger} \cdot \textbf{A}
! \end{align*}

! With :
! $\textbf{W}$ : eigenvectors of $\textbf{A}^2$
! $\tau$ : $\sqrt{-x}$
! $x$ : eigenvalues of $\textbf{A}^2$

! Input:
! | A(n,n) | double precision | antisymmetric matrix                                            |
! | n      | integer          | number of columns of the A matrix                               |
! | LDA    | integer          | specifies the leading dimension of A, must be at least max(1,n) |
! | LDR    | integer          | specifies the leading dimension of R, must be at least max(1,n) |

! Output:
! | R(n,n) | double precision | Rotation matrix                                      |
! | info   | integer          | if info = 0, the execution is successful             |
! |        |                  | if info = k, the k-th parameter has an illegal value |
! |        |                  | if info = -k, the algorithm failed                   |

! Internal:
! | B(n,n)        | double precision | B = A.A                                                       |
! | work(lwork,n) | double precision | work matrix for dysev, dimension max(1,lwork)                 |
! | lwork         | integer          | dimension of the syev work array >= max(1, 3n-1)              |
! | W(n,n)        | double precision | eigenvectors of B                                             |
! | e_val(n)      | double precision | eigenvalues of B                                              |
! | m_diag(n,n)   | double precision | diagonal matrix with the eigenvalues of B                     |
! | cos_tau(n,n)  | double precision | diagonal matrix with cos(tau) values                          |
! | sin_tau(n,n)  | double precision | diagonal matrix with sin cos(tau) values                      |
! | tau_m1(n,n)   | double precision | diagonal matrix with (tau)^-1 values                          |
! | part_1(n,n)   | double precision | matrix W.cos_tau.W^t                                          |
! | part_1a(n,n)  | double precision | matrix cos_tau.W^t                                            |
! | part_2(n,n)   | double precision | matrix W.tau_m1.sin_tau.W^t.A                                 |
! | part_2a(n,n)  | double precision | matrix W^t.A                                                  |
! | part_2b(n,n)  | double precision | matrix sin_tau.W^t.A                                          |
! | part_2c(n,n)  | double precision | matrix tau_m1.sin_tau.W^t.A                                   |
! | RR_t(n,n)     | double precision | R.R^t must be equal to the identity<=> R.R^t-1=0 <=> norm = 0 |
! | norm          | integer          | norm of R.R^t-1, must be equal to 0                           |
! | i,j           | integer          | indexes                                                       |

! Functions:
! | dnrm2  | double precision | Lapack function, compute the norm of a matrix |
! | disnan | logical          | Lapack function, check if an element is NaN   |



subroutine rotation_matrix(A,LDA,R,LDR,n,info,enforce_step_cancellation)

  implicit none

  BEGIN_DOC
  ! Rotation matrix to rotate the molecular orbitals.
  ! If the rotation is too large the transformation is not unitary and must be cancelled.
  END_DOC

  include 'pi.h'

  ! Variables

  ! in
  integer, intent(in)             :: n,LDA,LDR
  double precision, intent(inout) :: A(LDA,n)

  ! out
  double precision, intent(out)   :: R(LDR,n)
  integer, intent(out)            :: info
  logical, intent(out)            :: enforce_step_cancellation

  ! internal
  double precision, allocatable   :: B(:,:) 
  double precision, allocatable   :: work(:,:) 
  double precision, allocatable   :: W(:,:), e_val(:)
  double precision, allocatable   :: m_diag(:,:),cos_tau(:,:),sin_tau(:,:),tau_m1(:,:)
  double precision, allocatable   :: part_1(:,:),part_1a(:,:)
  double precision, allocatable   :: part_2(:,:),part_2a(:,:),part_2b(:,:),part_2c(:,:)
  double precision, allocatable   :: RR_t(:,:)
  integer                         :: i,j
  integer                         :: info2, lwork ! for dsyev
  double precision                :: norm, max_elem, max_elem_A, t1,t2,t3

  ! function
  double precision              :: dnrm2
  logical                       :: disnan

  print*,''
  print*,'---rotation_matrix---'

  call wall_time(t1)

  ! Allocation
  allocate(B(n,n))
  allocate(m_diag(n,n),cos_tau(n,n),sin_tau(n,n),tau_m1(n,n))
  allocate(W(n,n),e_val(n))
  allocate(part_1(n,n),part_1a(n,n))
  allocate(part_2(n,n),part_2a(n,n),part_2b(n,n),part_2c(n,n))
  allocate(RR_t(n,n))

! Pre-conditions

! Initialization
info=0
enforce_step_cancellation = .False.

! Size of matrix A must be at least 1 by 1
if (n<1) then
  info = 3
  print*, 'WARNING: invalid parameter 5'
  print*, 'n<1'
  return
endif

! Leading dimension of A must be >= n
if (LDA < n) then
  info = 25
  print*, 'WARNING: invalid parameter 2 or 5'
  print*, 'LDA < n'
  return
endif

! Leading dimension of A must be >= n
if (LDR < n) then
  info = 4
  print*, 'WARNING: invalid parameter 4'
  print*, 'LDR < n'
  return
endif

! Matrix elements of A must by non-NaN
do j = 1, n
  do i = 1, n
    if (disnan(A(i,j))) then
      info=1
      print*, 'WARNING: invalid parameter 1'
      print*, 'NaN element in A matrix'
      return
    endif
  enddo
enddo

do i = 1, n
  if (A(i,i) /= 0d0) then
    print*, 'WARNING: matrix A is not antisymmetric'
    print*, 'Non 0 element on the diagonal', i, A(i,i)
    call ABORT
  endif
enddo

do j = 1, n
  do i = 1, n
    if (A(i,j)+A(j,i)>1d-16) then
      print*, 'WANRING: matrix A is not antisymmetric'
      print*, 'A(i,j) /= - A(j,i):', i,j,A(i,j), A(j,i)
      print*, 'diff:', A(i,j)+A(j,i)
      call ABORT 
    endif
  enddo
enddo

! Fix for too big elements ! bad idea better to cancel if the error is too big
!do j = 1, n
!  do i = 1, n
!    A(i,j) = mod(A(i,j),2d0*pi)
!    if (dabs(A(i,j)) > pi) then
!      A(i,j) = 0d0
!    endif
!  enddo
!enddo

max_elem_A = 0d0
do j = 1, n
  do i = 1, n
    if (ABS(A(i,j)) > ABS(max_elem_A)) then
      max_elem_A = A(i,j)
    endif
  enddo
enddo
print*,'max element in A', max_elem_A

if (ABS(max_elem_A) > 2 * pi) then
   print*,''
   print*,'WARNING: ABS(max_elem_A) > 2 pi '
   print*,''
endif

! B=A.A
!     - Calculation of the matrix $\textbf{B} = \textbf{A}^2$
!     - Diagonalization of $\textbf{B}$ 
!     W, the eigenvectors
!     e_val, the eigenvalues


! Compute B=A.A

call dgemm('N','N',n,n,n,1d0,A,size(A,1),A,size(A,1),0d0,B,size(B,1))

! Copy B in W, diagonalization will put the eigenvectors in W
W=B

! Diagonalization of B
! Eigenvalues -> e_val
! Eigenvectors -> W
lwork = 3*n-1
allocate(work(lwork,n))

print*,'Starting diagonalization ...'

call dsyev('V','U',n,W,size(W,1),e_val,work,lwork,info2)

deallocate(work)

if (info2 == 0) then
   print*, 'Diagonalization : Done'
elseif (info2 < 0) then
   print*, 'WARNING: error in the diagonalization'
   print*, 'Illegal value of the ', info2,'-th parameter'
else
   print*, "WARNING: Diagonalization failed to converge"
endif

! Tau^-1, cos(tau), sin(tau)
!     $$\tau = \sqrt{-x}$$
!     - Calculation of $\cos(\tau)$  $\Leftrightarrow$ $\cos(\sqrt{-x})$
!     - Calculation of $\sin(\tau)$  $\Leftrightarrow$ $\sin(\sqrt{-x})$
!     - Calculation of $\tau^{-1}$ $\Leftrightarrow$ $(\sqrt{-x})^{-1}$
!     These matrices are diagonals

! Diagonal matrix m_diag
do j = 1, n
   if (e_val(j) >= -1d-12) then !0.d0) then !!! e_avl(i) must be < -1d-12 to avoid numerical problems
      e_val(j) = 0.d0
   else
      e_val(j) = - e_val(j)
   endif
enddo

m_diag = 0.d0
do i = 1, n
   m_diag(i,i) = e_val(i)
enddo

! cos_tau
do j = 1, n
   do i = 1, n
      if (i==j) then
         cos_tau(i,j) = dcos(dsqrt(e_val(i)))
      else
         cos_tau(i,j) = 0d0
      endif
   enddo
enddo

! sin_tau
do j = 1, n
   do i = 1, n
      if (i==j) then
         sin_tau(i,j) = dsin(dsqrt(e_val(i)))
      else
         sin_tau(i,j) = 0d0
      endif
   enddo
enddo

! Debug, display the cos_tau and sin_tau matrix
!if (debug) then
!   print*, 'cos_tau'
!   do i = 1, n
!      print*, cos_tau(i,:)
!   enddo
!   print*, 'sin_tau'
!   do i = 1, n
!      print*, sin_tau(i,:)
!   enddo
!endif

! tau^-1
do j = 1, n
   do i = 1, n
      if ((i==j) .and. (e_val(i) > 1d-16)) then!0d0)) then !!! Convergence problem can come from here if the threshold is too big/small
         tau_m1(i,j) = 1d0/(dsqrt(e_val(i)))
      else
         tau_m1(i,j) = 0d0
      endif
   enddo
enddo

max_elem = 0d0
do i = 1, n
   if (ABS(tau_m1(i,i)) > ABS(max_elem)) then
      max_elem = tau_m1(i,i)
   endif
enddo
print*,'max elem tau^-1:', max_elem

! Debug
!print*,'eigenvalues:'
!do i = 1, n 
!  print*, e_val(i)
!enddo

!Debug, display tau^-1
!if (debug) then
!   print*, 'tau^-1'
!   do i = 1, n
!      print*,tau_m1(i,:)
!   enddo
!endif

! Rotation matrix 
!     \begin{align*}
!     \textbf{R} = \textbf{W} \cos(\tau) \textbf{W}^{\dagger} + \textbf{W} \tau^{-1} \sin(\tau) \textbf{W}^{\dagger} \textbf{A}
!     \end{align*}
!     \begin{align*}
!     \textbf{Part1} = \textbf{W} \cos(\tau) \textbf{W}^{\dagger}
!     \end{align*}
!     \begin{align*}
!     \textbf{Part2} = \textbf{W} \tau^{-1} \sin(\tau) \textbf{W}^{\dagger} \textbf{A}
!     \end{align*}

!     First:
!     part_1 = dgemm(W, dgemm(cos_tau, W^t))
!     part_1a = dgemm(cos_tau, W^t)
!     part_1 = dgemm(W, part_1a)
!     And:
!     part_2= dgemm(W, dgemm(tau_m1, dgemm(sin_tau, dgemm(W^t, A))))
!     part_2a = dgemm(W^t, A)
!     part_2b = dgemm(sin_tau, part_2a)
!     part_2c = dgemm(tau_m1, part_2b)
!     part_2 = dgemm(W, part_2c)
!     Finally:
!     Rotation matrix, R = part_1+part_2

!     If $R$ is a rotation matrix:
!     $R.R^T=R^T.R=\textbf{1}$

! part_1
call dgemm('N','T',n,n,n,1d0,cos_tau,size(cos_tau,1),W,size(W,1),0d0,part_1a,size(part_1a,1))
call dgemm('N','N',n,n,n,1d0,W,size(W,1),part_1a,size(part_1a,1),0d0,part_1,size(part_1,1))

! part_2
call dgemm('T','N',n,n,n,1d0,W,size(W,1),A,size(A,1),0d0,part_2a,size(part_2a,1))
call dgemm('N','N',n,n,n,1d0,sin_tau,size(sin_tau,1),part_2a,size(part_2a,1),0d0,part_2b,size(part_2b,1))
call dgemm('N','N',n,n,n,1d0,tau_m1,size(tau_m1,1),part_2b,size(part_2b,1),0d0,part_2c,size(part_2c,1))
call dgemm('N','N',n,n,n,1d0,W,size(W,1),part_2c,size(part_2c,1),0d0,part_2,size(part_2,1))

! Rotation matrix R
R = part_1 + part_2

! Matrix check
! R.R^t and R^t.R must be equal to identity matrix
do j = 1, n
   do i=1,n
      if (i==j) then
         RR_t(i,j) = 1d0
      else
         RR_t(i,j) = 0d0
      endif
   enddo
enddo

call dgemm('N','T',n,n,n,1d0,R,size(R,1),R,size(R,1),-1d0,RR_t,size(RR_t,1))

norm = dnrm2(n*n,RR_t,1)
print*, 'Rotation matrix check, norm R.R^T = ', norm 

! Debug
!if (debug) then
!   print*, 'RR_t'
!   do i = 1, n
!      print*, RR_t(i,:)
!   enddo
!endif

! Post conditions

! Check if R.R^T=1
max_elem = 0d0
do j = 1, n
   do i = 1, n 
      if (ABS(RR_t(i,j)) > ABS(max_elem)) then
         max_elem = RR_t(i,j)
      endif
   enddo
enddo

print*, 'Max error in R.R^T:', max_elem
print*, 'e_val(1):', e_val(1)
print*, 'e_val(n):', e_val(n)
print*, 'max elem in A:', max_elem_A

if (ABS(max_elem) > 1d-12) then
  print*, 'WARNING: max error in R.R^T > 1d-12'
  print*, 'Enforce the step cancellation'
  enforce_step_cancellation = .True.
endif

! Matrix elements of R must by non-NaN
do j = 1,n
   do i = 1,LDR
      if (disnan(R(i,j))) then
         info = 666
         print*, 'NaN in rotation matrix'
         call ABORT
      endif
   enddo
enddo

! Display
!if (debug) then
!   print*,'Rotation matrix :'
!   do i = 1, n
!      write(*,'(100(F10.5))') R(i,:)
!   enddo
!endif

! Deallocation, end

deallocate(B)
  deallocate(m_diag,cos_tau,sin_tau,tau_m1)
  deallocate(W,e_val)
  deallocate(part_1,part_1a)
  deallocate(part_2,part_2a,part_2b,part_2c)
  deallocate(RR_t)

  call wall_time(t2)
  t3 = t2-t1
  print*,'Time in rotation matrix:', t3

  print*,'---End rotation_matrix---'

end subroutine
