! Diagonalization of the hessian

! Just a matrix diagonalization using Lapack

! Input:
! | n      | integer          | mo_num*(mo_num-1)/2 |
! | H(n,n) | double precision | hessian             |

! Output:
! | e_val(n) | double precision | eigenvalues of the hessian  |
! | w(n,n)   | double precision | eigenvectors of the hessian |

! Internal:
! | nb_negative_nv | integer          | number of negative eigenvalues                  |
! | lwork          | integer          | for Lapack                                      |
! | work(lwork,n)  | double precision | temporary array for Lapack                      |
! | info           | integer          | if 0 -> ok, else problem in the diagonalization |
! | i,j            | integer          | dummy indexes                                   |


subroutine diagonalization_hessian(n,H,e_val,w)

  include 'constants.h'

  implicit none

  ! Variables

  ! in
  integer, intent(in) :: n
  double precision, intent(in) :: H(n,n)

  ! out
  double precision, intent(out) :: e_val(n), w(n,n)

  ! internal
  double precision, allocatable :: work(:,:)
  integer, allocatable          :: key(:)
  integer                       :: info,lwork
  integer                       :: i,j
  integer                       :: nb_negative_vp
  double precision              :: t1,t2,t3,max_elem

  print*,''
  print*,'---Diagonalization_hessian---'

  call wall_time(t1)

  if (optimization_method == 'full') then
    ! Allocation
    ! For Lapack
    lwork=3*n-1

    allocate(work(lwork,n))

    ! Calculation

    ! Copy the hessian matrix, the eigenvectors will be store in W
    W=H

    ! Diagonalization of the hessian
    call dsyev('V','U',n,W,size(W,1),e_val,work,lwork,info)

    if (info /= 0) then
        print*, 'Error diagonalization : diagonalization_hessian'
        print*, 'info = ', info
        call ABORT
    endif

    if (debug) then
      print *, 'vp Hess:'
      write(*,'(100(F10.5))')  real(e_val(:))
    endif

    ! Number of negative eigenvalues
    max_elem = 0d0
    nb_negative_vp = 0
    do i = 1, n
      if (e_val(i) < 0d0) then
        nb_negative_vp = nb_negative_vp + 1
        if (e_val(i) < max_elem) then
          max_elem = e_val(i)
        endif
        !print*,'e_val < 0 :', e_val(i)
      endif
    enddo
    print*,'Number of negative eigenvalues:', nb_negative_vp
    print*,'Lowest eigenvalue:',max_elem

    !nb_negative_vp = 0
    !do i = 1, n
    !  if (e_val(i) < -thresh_eig) then
    !    nb_negative_vp = nb_negative_vp + 1
    !  endif
    !enddo
    !print*,'Number of negative eigenvalues <', -thresh_eig,':', nb_negative_vp  

    ! Deallocation
    deallocate(work)

  elseif (optimization_method == 'diag') then
    ! Diagonalization of the diagonal hessian by hands
    allocate(key(n))
    
    do i = 1, n
      e_val(i) = H(i,i)
    enddo

    ! Key list for dsort
    do i = 1, n 
      key(i) = i
    enddo

    ! Sort of the eigenvalues
    call dsort(e_val, key, n)

    ! Eigenvectors
    W = 0d0
    do i = 1, n
      j = key(i)
      W(j,i) = 1d0
    enddo

    deallocate(key)
  else
    print*,'Diagonalization_hessian, abort'
    call abort
  endif

  call wall_time(t2)
  t3 = t2 - t1
  print*,'Time in diagonalization_hessian:', t3

  print*,'---End diagonalization_hessian---'

end subroutine
