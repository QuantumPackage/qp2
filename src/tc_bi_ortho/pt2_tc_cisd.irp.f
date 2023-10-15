program pt2_tc_cisd

  BEGIN_DOC
  !
  ! TODO : Reads psi_det in the EZFIO folder and prints out the left- and right-eigenvectors together 
  !        with the energy. Saves the left-right wave functions at the end. 
  !
  END_DOC

  implicit none

  my_grid_becke = .True.
  PROVIDE tc_grid1_a tc_grid1_r
  my_n_pt_r_grid = tc_grid1_r
  my_n_pt_a_grid = tc_grid1_a
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  read_wf = .True.
  touch read_wf

  print*, ' nb of states = ', N_states
  print*, ' nb of det    = ', N_det
  call routine_diag()

  call routine
end

subroutine routine
 implicit none
 integer :: i,h1,p1,h2,p2,s1,s2,degree
 double precision :: h0i,hi0,e00,ei,delta_e
 double precision :: norm,e_corr,coef,e_corr_pos,e_corr_neg,e_corr_abs

 integer                        :: exc(0:2,2,2)
 double precision               :: phase
 double precision :: eh1,ep1,eh2,ep2

 norm = 0.d0 
 e_corr = 0.d0
 e_corr_abs = 0.d0
 e_corr_pos = 0.d0
 e_corr_neg = 0.d0
 call htilde_mu_mat_bi_ortho_tot_slow(psi_det(1,1,1), psi_det(1,1,1), N_int, e00) 
 do i = 2, N_det
  call htilde_mu_mat_bi_ortho_tot_slow(psi_det(1,1,i), psi_det(1,1,1), N_int, hi0) 
  call htilde_mu_mat_bi_ortho_tot_slow(psi_det(1,1,1), psi_det(1,1,i), N_int, h0i) 
  call htilde_mu_mat_bi_ortho_tot_slow(psi_det(1,1,i), psi_det(1,1,i), N_int, ei) 
  call get_excitation_degree(psi_det(1,1,1), psi_det(1,1,i),degree,N_int)
  call get_excitation(psi_det(1,1,1), psi_det(1,1,i),exc,degree,phase,N_int)
  call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
  eh1 = Fock_matrix_tc_diag_mo_tot(h1)
  ep1 = Fock_matrix_tc_diag_mo_tot(p1)
  delta_e = eh1 - ep1
  if (degree==2)then
   eh2 = Fock_matrix_tc_diag_mo_tot(h2)
   ep2 = Fock_matrix_tc_diag_mo_tot(p2)
   delta_e +=  eh2 - ep2
  endif
!  delta_e = e00 - ei
  coef = hi0/delta_e
  norm += coef*coef
  e_corr = coef* h0i
  if(e_corr.lt.0.d0)then
   e_corr_neg += e_corr
  elseif(e_corr.gt.0.d0)then
   e_corr_pos += e_corr
  endif
  e_corr_abs += dabs(e_corr)
 enddo
 print*,'e_corr_abs = ',e_corr_abs
 print*,'e_corr_pos = ',e_corr_pos
 print*,'e_corr_neg = ',e_corr_neg
 print*,'norm       = ',dsqrt(norm)

end

subroutine routine_diag()

  implicit none
  integer          :: i, j, k
  double precision :: dE

  ! provide eigval_right_tc_bi_orth
  ! provide overlap_bi_ortho
  ! provide htilde_matrix_elmt_bi_ortho

  if(N_states .eq. 1) then

    print*,'eigval_right_tc_bi_orth   = ',eigval_right_tc_bi_orth(1)
    print*,'e_tc_left_right           = ',e_tc_left_right
    print*,'e_tilde_bi_orth_00        = ',e_tilde_bi_orth_00
    print*,'e_pt2_tc_bi_orth          = ',e_pt2_tc_bi_orth
    print*,'e_pt2_tc_bi_orth_single   = ',e_pt2_tc_bi_orth_single
    print*,'e_pt2_tc_bi_orth_double   = ',e_pt2_tc_bi_orth_double
    print*,'***'                      
    print*,'e_corr_bi_orth            = ',e_corr_bi_orth
    print*,'e_corr_bi_orth_proj       = ',e_corr_bi_orth_proj
    print*,'e_corr_bi_orth_proj_abs   = ',e_corr_bi_orth_proj_abs
    print*,'e_corr_single_bi_orth     = ',e_corr_single_bi_orth
    print*,'e_corr_double_bi_orth     = ',e_corr_double_bi_orth
    print*,'e_corr_single_bi_orth_abs = ',e_corr_single_bi_orth_abs
    print*,'e_corr_double_bi_orth_abs = ',e_corr_double_bi_orth_abs
    print*,'Left/right eigenvectors'
    do i = 1,N_det
      write(*,'(I5,X,(100(F12.7,X)))')i,leigvec_tc_bi_orth(i,1),reigvec_tc_bi_orth(i,1),leigvec_tc_bi_orth(i,1)*reigvec_tc_bi_orth(i,1)
    enddo

  else

    print*,'eigval_right_tc_bi_orth : '
    do i = 1, N_states
      print*, i, eigval_right_tc_bi_orth(i)
    enddo

    print*,''
    print*,'******************************************************'
    print*,'TC Excitation energies (au)                     (eV)'
    do i = 2, N_states
      dE = eigval_right_tc_bi_orth(i) - eigval_right_tc_bi_orth(1)
      print*, i, dE, dE/0.0367502d0
    enddo
    print*,''

  endif

end



