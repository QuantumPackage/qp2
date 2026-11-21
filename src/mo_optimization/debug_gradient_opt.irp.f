! Debug the gradient

! *Program to check the gradient*

! The program compares the result of the first and last code for the
! gradient.

! Provided:
! | mo_num | integer | number of MOs |

! Internal:
! | n          | integer          | number of orbitals pairs (p,q) p<q |
! | v_grad(n)  | double precision | Original gradient                  |
! | v_grad2(n) | double precision | Gradient                           |
! | i          | integer          | index                              |
! | threshold  | double precision | threshold for the errors           |
! | max_error  | double precision | maximal error in the gradient      |
! | nb_error   | integer          | number of error in the gradient    |


program debug_gradient_opt
  BEGIN_DOC
  ! Debug program for the orbital gradient, for all mo_num orbitals
  ! (n = mo_num*(mo_num-1)/2 rotation pairs). Compares the gradient
  ! of the energy with respect to orbital rotations computed by the
  ! current implementation (gradient_opt) with the original one
  ! (first_gradient_opt), and lists the elements whose difference
  ! exceeds the 1d-12 threshold.
  END_DOC

  implicit none

  ! Variables

  double precision, allocatable :: v_grad(:), v_grad2(:)
  integer                       :: n
  integer                       :: i
  double precision              :: threshold
  double precision              :: max_error, max_elem
  integer                       :: nb_error

  ! Definition of n
  n = mo_num*(mo_num-1)/2

  PROVIDE all_mo_integrals ! Check for suppression

  ! Allocation
  allocate(v_grad(n), v_grad2(n))

  ! Calculation

  call diagonalize_ci

  ! Gradient
  call first_gradient_opt(n,v_grad)
  call gradient_opt(n,v_grad2,max_elem)

  v_grad = v_grad - v_grad2
  nb_error = 0
  max_error = 0d0
  threshold = 1d-12

  do i = 1, n
    if (ABS(v_grad(i)) > threshold) then
       print*,v_grad(i)
       nb_error = nb_error + 1

       if (ABS(v_grad(i)) > max_error) then
         max_error = v_grad(i)
       endif

    endif
  enddo

  print*,''
  print*,'Check the gradient'
  print*,'Threshold :', threshold
  print*,'Nb error :', nb_error
  print*,'Max error :', max_error

  ! Deallocation
  deallocate(v_grad,v_grad2)

end program
