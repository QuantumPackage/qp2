
! Criterion FB (old)

! The criterion is just computed as

! \begin{align*}
! C = - \sum_i^{mo_{num}} (<i|x|i>^2 + <i|y|i>^2 + <i|z|i>^2)
! \end{align*}

! The minus sign is here in order to minimize this criterion

! Output:
! | criterion | double precision | criterion for the Foster-Boys localization |


subroutine criterion_FB_old(criterion)

  implicit none

  BEGIN_DOC
  ! Compute the Foster-Boys localization criterion
  END_DOC

  double precision, intent(out) :: criterion
  integer                       :: i

  ! Criterion (= \sum_i <i|r|i>^2 )
  criterion = 0d0
  do i = 1, mo_num
    criterion = criterion + mo_dipole_x(i,i)**2 + mo_dipole_y(i,i)**2 + mo_dipole_z(i,i)**2
  enddo
  criterion = - criterion

end subroutine

! Criterion FB

subroutine criterion_FB(tmp_list_size, tmp_list, criterion)

  implicit none

  BEGIN_DOC
  ! Compute the Foster-Boys localization criterion
  END_DOC

  integer, intent(in)           :: tmp_list_size, tmp_list(tmp_list_size)
  double precision, intent(out) :: criterion
  integer                       :: i, tmp_i

  ! Criterion (= - \sum_i <i|r|i>^2 )
  criterion = 0d0
  do tmp_i = 1, tmp_list_size
    i = tmp_list(tmp_i)
    criterion = criterion + mo_dipole_x(i,i)**2 + mo_dipole_y(i,i)**2 + mo_dipole_z(i,i)**2
  enddo
  criterion = - criterion

end subroutine

subroutine theta_FB(l, n, m_x, max_elem)

  include 'pi.h'

  BEGIN_DOC
  ! Compute the angles to minimize the Foster-Boys criterion by using pairwise rotations of the MOs
  ! Warning: you must give - the angles to build the rotation matrix...
  END_DOC

  implicit none

  integer, intent(in)           :: n, l(n)
  double precision, intent(out) :: m_x(n,n), max_elem

  integer                       :: i,j, tmp_i, tmp_j
  double precision, allocatable :: cos4theta(:,:), sin4theta(:,:)
  double precision, allocatable :: A(:,:), B(:,:), beta(:,:), gamma(:,:)
  integer                       :: idx_i,idx_j

  allocate(cos4theta(n, n), sin4theta(n, n))
  allocate(A(n,n), B(n,n), beta(n,n), gamma(n,n))

  do tmp_j = 1, n
    j = l(tmp_j)
    do tmp_i = 1, n
      i = l(tmp_i)
      A(tmp_i,tmp_j) = mo_dipole_x(i,j)**2 - 0.25d0 * (mo_dipole_x(i,i) - mo_dipole_x(j,j))**2 &
                     + mo_dipole_y(i,j)**2 - 0.25d0 * (mo_dipole_y(i,i) - mo_dipole_y(j,j))**2 &
                     + mo_dipole_z(i,j)**2 - 0.25d0 * (mo_dipole_z(i,i) - mo_dipole_z(j,j))**2
    enddo
    A(j,j) = 0d0
  enddo

  do tmp_j = 1, n
    j = l(tmp_j)
    do tmp_i = 1, n
      i = l(tmp_i)
      B(tmp_i,tmp_j) = mo_dipole_x(i,j) * (mo_dipole_x(i,i) - mo_dipole_x(j,j)) &
                     + mo_dipole_y(i,j) * (mo_dipole_y(i,i) - mo_dipole_y(j,j)) &
                     + mo_dipole_z(i,j) * (mo_dipole_z(i,i) - mo_dipole_z(j,j))
    enddo
  enddo

  !do tmp_j = 1, n
  !  j = l(tmp_j)
  !  do tmp_i = 1, n
  !    i = l(tmp_i)
  !    beta(tmp_i,tmp_j) =  (mo_dipole_x(i,i) - mo_dipole_x(j,j)) - 4d0 * mo_dipole_x(i,j)**2 &
  !               + (mo_dipole_y(i,i) - mo_dipole_y(j,j)) - 4d0 * mo_dipole_y(i,j)**2 &
  !               + (mo_dipole_z(i,i) - mo_dipole_z(j,j)) - 4d0 * mo_dipole_z(i,j)**2
  !  enddo
  !enddo

  !do tmp_j = 1, n
  !  j = l(tmp_j)
  !  do tmp_i = 1, n
  !    i = l(tmp_i)
  !    gamma(tmp_i,tmp_j) = 4d0 * ( mo_dipole_x(i,j) * (mo_dipole_x(i,i) - mo_dipole_x(j,j)) &
  !                       + mo_dipole_y(i,j) * (mo_dipole_y(i,i) - mo_dipole_y(j,j)) &
  !                       + mo_dipole_z(i,j) * (mo_dipole_z(i,i) - mo_dipole_z(j,j)))
  !  enddo
  !enddo

  !
  !do j = 1, n
  !  do i = 1, n
  !    cos4theta(i,j) = - A(i,j) / dsqrt(A(i,j)**2 + B(i,j)**2)
  !  enddo
  !enddo

  !do j = 1, n
  !  do i = 1, n
  !    sin4theta(i,j) = B(i,j) / dsqrt(A(i,j)**2 + B(i,j)**2)
  !  enddo
  !enddo

  ! Theta
  do j = 1, n
    do i = 1, n
      m_x(i,j) = 0.25d0 * atan2(B(i,j), -A(i,j))
      !m_x(i,j) = 0.25d0 * atan2(sin4theta(i,j), cos4theta(i,j))
    enddo
  enddo

  ! Enforce a perfect antisymmetry
  do j = 1, n-1
    do i = j+1, n
      m_x(j,i) = - m_x(i,j)
    enddo
  enddo
  do i = 1, n
    m_x(i,i) = 0d0
  enddo

  ! Max
  max_elem = 0d0
  do j = 1, n-1
    do i = j+1, n
      if (dabs(m_x(i,j)) > dabs(max_elem)) then
        max_elem = m_x(i,j)
        !idx_i = i
        !idx_j = j
      endif
    enddo
  enddo

  ! Debug
  !print*,''
  !print*,'sin/B'
  !do i = 1, n
  !  write(*,'(100F10.4)') sin4theta(i,:)
  !  !B(i,:)
  !enddo
  !print*,'cos/A'
  !do i = 1, n
  !  write(*,'(100F10.4)') cos4theta(i,:)
  !  !A(i,:)
  !enddo
  !print*,'X'
  !!m_x = 0d0
  !!m_x(idx_i,idx_j) = max_elem
  !!m_x(idx_j,idx_i) = -max_elem
  !do i = 1, n
  !  write(*,'(100F10.4)') m_x(i,:)
  !enddo
  !print*,idx_i,idx_j,max_elem

  max_elem = dabs(max_elem)
  
  deallocate(cos4theta, sin4theta)
  deallocate(A,B,beta,gamma)
  
end

