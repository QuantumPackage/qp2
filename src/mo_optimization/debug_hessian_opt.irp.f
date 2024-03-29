! Debug the hessian

! *Program to check the hessian matrix*

! The program compares the result of the first and last code for the
! hessian. First of all the 4D hessian and after the 2D hessian.

! Provided:
! | mo_num | integer | number of MOs |

! Internal:
! | n                                 | integer          | number of orbitals pairs (p,q) p<q |
! | H(n,n)                            | double precision | Original hessian matrix (2D)       |
! | H2(n,n)                           | double precision | Hessian matrix (2D)                |
! | h_f(mo_num,mo_num,mo_num,mo_num)  | double precision | Original hessian matrix (4D)       |
! | h_f2(mo_num,mo_num,mo_num,mo_num) | double precision | Hessian matrix (4D)                |
! | method                            | integer          | - 1: full hessian                  |
! |                                   |                  | - 2: diagonal hessian              |
! | i,j,p,q,k                         | integer          | indexes                            |
! | threshold                         | double precision | threshold for the errors           |
! | max_error                         | double precision | maximal error in the 4D hessian    |
! | max_error_H                       | double precision | maximal error in the 2D hessian    |
! | nb_error                          | integer          | number of errors in the 4D hessian |
! | nb_error_H                        | integer          | number of errors in the 2D hessian |


program debug_hessian

  implicit none

  ! Variables

  double precision, allocatable :: H(:,:),H2(:,:), h_f(:,:,:,:), h_f2(:,:,:,:)
  integer                       :: n
  integer                       :: i,j,k,l
  double precision              :: max_error, max_error_H
  integer                       :: nb_error, nb_error_H
  double precision              :: threshold
  
  ! Definition of n  
  n = mo_num*(mo_num-1)/2

  PROVIDE mo_two_e_integrals_in_map 

  ! Allocation
  allocate(H(n,n),H2(n,n))  
  allocate(h_f(mo_num,mo_num,mo_num,mo_num),h_f2(mo_num,mo_num,mo_num,mo_num))

  ! Calculation
  
  ! Hessian 
  if (optimization_method == 'full') then 

    print*,'Use the full hessian matrix'
    call hessian_opt(n,H,h_f)
    call first_hessian_opt(n,H2,h_f2)

    ! Difference
    h_f = h_f - h_f2
    H = H - H2
    max_error = 0d0
    nb_error = 0    
    threshold = 1d-12

    do l = 1, mo_num
      do k= 1, mo_num
        do j = 1, mo_num
          do i = 1, mo_num
            if (ABS(h_f(i,j,k,l)) > threshold) then
              print*,h_f(i,j,k,l)
              nb_error = nb_error + 1
              if (ABS(h_f(i,j,k,l)) > ABS(max_error)) then
                max_error = h_f(i,j,k,l)
              endif
            endif
          enddo
        enddo
      enddo
    enddo

   max_error_H = 0d0
   nb_error_H = 0

   do j = 1, n
     do i = 1, n
       if (ABS(H(i,j)) > threshold) then
         print*, H(i,j)
         nb_error_H = nb_error_H + 1

         if (ABS(H(i,j)) > ABS(max_error_H)) then
           max_error_H = H(i,j)
         endif

       endif
     enddo
   enddo 

  elseif (optimization_method == 'diag') then

    print*, 'Use the diagonal hessian matrix'
    call diag_hessian_opt(n,H,h_f)
    call first_diag_hessian_opt(n,H2,h_f2)

    h_f = h_f - h_f2
    max_error = 0d0
    nb_error = 0
    threshold = 1d-12

    do l = 1, mo_num
      do k = 1, mo_num
        do j = 1, mo_num
          do i = 1, mo_num

            if (ABS(h_f(i,j,k,l)) > threshold) then

              print*,h_f(i,j,k,l)
              nb_error = nb_error + 1

              if (ABS(h_f(i,j,k,l)) > ABS(max_error)) then
                max_error = h_f(i,j,k,l)
              endif

            endif

          enddo
        enddo
      enddo
    enddo

    h=H-H2
  
    max_error_H = 0d0
    nb_error_H = 0
 
    do j = 1, n
      do i = 1, n
        if (ABS(H(i,j)) > threshold) then
          print*, H(i,j)
          nb_error_H = nb_error_H + 1
 
          if (ABS(H(i,j)) > ABS(max_error_H)) then
            max_error_H = H(i,j)
          endif
 
        endif
      enddo
    enddo
  
  else
    print*,'Unknown optimization_method, please select full, diag'
    call abort
  endif
  
  print*,''
  if (optimization_method == 'full') then
    print*,'Check the full hessian'
  else
    print*,'Check the diagonal hessian'
  endif
   
  print*,'Threshold :', threshold
  print*,'Nb error :', nb_error
  print*,'Max error :', max_error
  print*,''
  print*,'Nb error_H :', nb_error_H
  print*,'Max error_H :', max_error_H
 
  ! Deallocation
  deallocate(H,H2,h_f,h_f2)

end program
