! Vect to antisymmetric matrix using mat_to_vec_index

! Vector to antisymmetric matrix transformation using mat_to_vec_index
! subroutine.

! Can be done in OMP (for the first part and with omp critical for the second)


subroutine vec_to_mat_v2(n,m,v_x,m_x)

  BEGIN_DOC
  ! Vector to antisymmetric matrix
  END_DOC
  
  implicit none
  
  integer, intent(in) :: n,m
  double precision, intent(in) :: v_x(n)
  double precision, intent(out) :: m_x(m,m)

  integer :: i,j,k

  ! 1D -> 2D lower diagonal
  m_x = 0d0
  do j = 1, m - 1
    do i = j + 1, m
      call mat_to_vec_index(i,j,k)
      m_x(i,j) = v_x(k)
    enddo
  enddo
  
  ! Antisym
  do i = 1, m - 1
    do j = i + 1, m
      m_x(i,j) = - m_x(j,i) 
    enddo
  enddo

end
