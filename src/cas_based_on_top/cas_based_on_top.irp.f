program cas_based_on_top_density
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
 implicit none
 provide total_cas_on_top_density
end
