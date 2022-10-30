 program tc_natorb_bi_ortho
   implicit none
   BEGIN_DOC
 ! TODO : Put the documentation of the program here
   END_DOC
   print *, 'Hello world'
   my_grid_becke = .True.
   my_n_pt_r_grid = 30
   my_n_pt_a_grid = 50
   read_wf = .True.
   touch read_wf
   touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid
   call print_energy_and_mos
   call save_tc_natorb
!   call minimize_tc_orb_angles
 end
 
 subroutine save_tc_natorb 
  implicit none
  print*,'Saving the natorbs '
  provide natorb_tc_leigvec_ao natorb_tc_reigvec_ao
  call ezfio_set_bi_ortho_mos_mo_l_coef(natorb_tc_leigvec_ao)
  call ezfio_set_bi_ortho_mos_mo_r_coef(natorb_tc_reigvec_ao)
  call save_ref_determinant_nstates_1
  call ezfio_set_determinants_read_wf(.False.)
 end
 
 subroutine save_ref_determinant_nstates_1
   implicit none
   use bitmasks
   double precision               :: buffer(1,N_states)
   buffer = 0.d0
   buffer(1,1) = 1.d0
   call save_wavefunction_general(1,1,ref_bitmask,1,buffer)
 end
