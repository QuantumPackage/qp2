program test_dav
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  print *, 'Hello world'
  read_wf = .True.
  touch read_wf
  PROVIDE threshold_davidson nthreads_davidson
  call routine
  call test_dav_dress
end

subroutine routine
 implicit none
  BEGIN_DOC
! Example driver for the two non-dressed Davidson routines
! of this module: h_mat is the provider H_matrix_all_dets
! (|H| on the determinants of psi_det) with nuclear_repulsion
! added to the diagonal; davidson_general is called with
! h_mat, and davidson_general_ext_rout with hcalc_template
! ($v = H u$). Prints the eigenvalues after each call.
  END_DOC
 double precision, allocatable :: u_in(:,:), H_jj(:), energies(:),h_mat(:,:)
 integer :: dim_in,sze,N_st,N_st_diag_in
 logical :: converged 
 integer :: i,j
 external hcalc_template
 N_st = N_states
 N_st_diag_in = N_states_diag
 sze = N_det
 dim_in = sze
 !!!! MARK THAT u_in mut dimensioned with "N_st_diag_in" as a second dimension 
 allocate(u_in(dim_in,N_st_diag_in),H_jj(sze),h_mat(sze,sze),energies(N_st_diag_in))
 u_in = 0.d0
 do i = 1, N_st
  u_in(1,i) = 1.d0
 enddo
 !!! Matrix "h_mat" is the matrix we want to diagonalize with the first routine 
 !!! "davidson_general" 
 do i = 1, sze
  do j = 1, sze
   h_mat(j,i) = H_matrix_all_dets(j,i)
  enddo
  H_jj(i) = H_mat(i,i) + nuclear_repulsion
  h_mat(i,i) = H_mat(i,i) + nuclear_repulsion
 enddo
 provide nthreads_davidson
 call davidson_general(u_in,H_jj,energies,dim_in,sze,N_st,N_st_diag_in,converged,h_mat)
 print*,'energies = ',energies 
 !!! hcalc_template is the routine that computes v = H u 
 !!! and you can use the routine "davidson_general_ext_rout" 
 call davidson_general_ext_rout(u_in,H_jj,energies,sze,N_st,N_st_diag_in,converged,hcalc_template)
 print*,'energies = ',energies 
end


subroutine test_dav_dress
 implicit none
  BEGIN_DOC
! Dressed variant of the test: H_jj is the diagonal of the
! provider H_matrix_all_dets plus nuclear_repulsion, and
! davidson_general_ext_rout_dressed is called with
! hcalc_template ($v = H u$). dressing_state is 0, so the
! zero dressing_vec does not enter H_jj in the called
! routine. Prints the first computed energy.
  END_DOC
 double precision, allocatable :: u_in(:,:), H_jj(:), energies(:)
 integer :: sze,N_st,N_st_diag_in,dressing_state
 logical :: converged 
 integer :: i,j
 external hcalc_template
 double precision, allocatable :: dressing_vec(:)
 integer :: idress
 N_st = N_states
 N_st_diag_in = N_states_diag
 sze = N_det
 dressing_state = 0
 idress = 1
 !!!! MARK THAT u_in mut dimensioned with "N_st_diag_in" as a second dimension 
 allocate(u_in(sze,N_st_diag_in),H_jj(sze),energies(N_st_diag_in))
 allocate(dressing_vec(sze))
 dressing_vec = 0.d0
 u_in = 0.d0
 do i = 1, N_st
  u_in(1,i) = 1.d0
 enddo
 do i = 1, sze
  H_jj(i) = H_matrix_all_dets(i,i) + nuclear_repulsion
 enddo
 print*,'dressing davidson '
 call davidson_general_ext_rout_dressed(u_in,H_jj,energies,sze,N_st,N_st_diag_in,dressing_state,dressing_vec,idress,converged,hcalc_template)
 print*,'energies(1) = ',energies(1)

end
