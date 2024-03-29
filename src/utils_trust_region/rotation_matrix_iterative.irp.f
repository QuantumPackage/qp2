! Rotation matrix with the iterative method

! \begin{align*}
! \textbf{R} = \sum_{k=0}^{\infty} \frac{1}{k!} \textbf{X}^k
! \end{align*}

! !!! Doesn't work !!!


subroutine rotation_matrix_iterative(m,X,R)

  implicit none

  ! in
  integer, intent(in) :: m
  double precision, intent(in) :: X(m,m)
  
  ! out
  double precision, intent(out) :: R(m,m)
 
  ! internal
  double precision :: max_elem, pre_factor
  double precision :: t1,t2,t3
  integer :: k,l,i,j
  logical :: not_converged 
  double precision, allocatable :: RRT(:,:), A(:,:), B(:,:)

  ! Functions
  integer :: factorial

  print*,'---rotation_matrix_iterative---'
  call wall_time(t1)

  allocate(RRT(m,m),A(m,m),B(m,m))

  ! k = 0
  R = 0d0
  do i = 1, m
    R(i,i) = 1d0
  enddo
 
  ! k = 1
  R = R + X
  
  k = 2

  not_converged = .True.
 
  do while (not_converged) 

    pre_factor = 1d0/DBLE(factorial(k))
    if (pre_factor < 1d-15) then
      print*,'pre factor=', pre_factor,'< 1d-15, exit'
      exit
    endif

    A = X
    B = 0d0
    do l = 1, k-1
      call dgemm('N','N',m,m,m,1d0,X,size(X,1),A,size(A,1),0d0,B,size(B,1))
      A = B
    enddo

    !print*,'B'
    !do i = 1, m
    !  print*,B(i,:) * 1d0/DBLE(factorial(k))
    !enddo

    R = R + pre_factor * B

    k = k + 1
    call dgemm('T','N',m,m,m,1d0,R,size(R,1),R,size(R,1),0d0,RRT,size(RRT,1))

    !print*,'R'
    !do i = 1, m
    !  write(*,'(10(ES12.5))') R(i,:)
    !enddo

    do i = 1, m
      RRT(i,i) = RRT(i,i) - 1d0
    enddo

    !print*,'RRT'
    !do i = 1, m
    !  write(*,'(10(ES12.5))') RRT(i,:)
    !enddo

    max_elem = 0d0
    do j = 1, m
      do i = 1, m
        if (dabs(RRT(i,j)) > max_elem) then
          max_elem = dabs(RRT(i,j))
        endif
      enddo
    enddo

    print*, 'Iteration:', k
    print*, 'Max error in R:', max_elem
  
    if (max_elem < 1d-12) then
      not_converged = .False.
    endif

  enddo

  deallocate(RRT,A,B)

  call wall_time(t2)
  t3 = t2 - t1
  print*,'Time in rotation matrix iterative:', t3
  print*,'---End roration_matrix_iterative---'


print*,'Does not work yet, abort'
call abort

end

! Factorial

function factorial(n)

  implicit none

  integer, intent(in) :: n
  integer :: factorial, k

  factorial = 1

  do k = 1, n
    factorial = factorial * k
  enddo
  
end
