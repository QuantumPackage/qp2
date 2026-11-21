program cas_based_on_top
  implicit none
  BEGIN_DOC
! TODO : Small example to use the different quantities in this plugin
  END_DOC

  !! You force QP2 to read the wave function in the EZFIO folder
  !! It is assumed that all Slater determinants in the wave function 
  !! belongs to an active space defined by core, inactive and active list of orbitals 
  read_wf = .True.
  touch read_wf
!  call routine_test_cas_based_on_top_density
 call routine
end

subroutine routine
  BEGIN_DOC
! Example of the cas_based_on_top program. The PROVIDE
! statements compute (and print) total_cas_on_top_density and
! average_on_top. It then prints accu, the sum over active
! orbital pairs (i,j) of act_2_rdm_ab_mo(i,j,i,j,1).
! The intended sum is over the active-space indices i,j =
! 1..n_act_orb of act_2_rdm_ab_mo(i,j,i,j,1); the current code
! uses the global MO numbers from list_act, a suspected bug
! (author-confirmed 2026-10-31).
  END_DOC
 implicit none
 provide total_cas_on_top_density
 provide average_on_top
 integer :: i,j,ii,jj
 double precision :: accu
 accu = 0.d0
 do ii = 1, n_act_orb
  i = list_act(ii)
  do jj = 1, n_act_orb
   j = list_act(jj)
   accu += act_2_rdm_ab_mo(i,j,i,j,1)
  enddo
 enddo
 print*,'accu = ',accu
end
