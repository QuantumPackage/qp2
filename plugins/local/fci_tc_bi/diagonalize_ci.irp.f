
! ---

subroutine diagonalize_CI_tc_bi_ortho(ndet, E_tc, norm, pt2_data, print_pt2)

  BEGIN_DOC
  !  Replace the coefficients of the CI states by the coefficients of the
  !  eigenstates of the CI matrix
  END_DOC

  use selection_types
  implicit none
  integer,          intent(inout) :: ndet       ! number of determinants from before 
  double precision, intent(inout) :: E_tc(N_states), norm(N_states) ! E and norm from previous wave function 
  type(pt2_type)  , intent(in)    :: pt2_data   ! PT2 from previous wave function 
  logical,          intent(in)    :: print_pt2
  integer                         :: i, j,k
  double precision:: pt2_minus,pt2_plus,pt2_tot, pt2_abs,pt1_norm,rpt2_tot
  double precision :: error_pt2_minus, error_pt2_plus, error_pt2_tot, error_pt2_abs

  PROVIDE mo_l_coef mo_r_coef

  print*,'*****'
  print*,'New wave function information'
  print*,'N_det tc               = ',N_det
  do k = 1, N_states
   print*,'************'
   print*,'State ',k
   pt2_plus  = pt2_data % variance(k)
   pt2_minus = pt2_data % pt2(k)
   pt2_abs   = pt2_plus - pt2_minus 
   pt2_tot   = pt2_plus + pt2_minus 
!   error_pt2_minus = pt2_data_err % pt2(k)
!   error_pt2_plus  = pt2_data_err % variance(k)
!   error_pt2_tot = dsqrt(error_pt2_minus**2+error_pt2_plus**2)
!   error_pt2_abs = error_pt2_tot ! same variance because independent variables 
 
   pt1_norm = pt2_data % overlap(k,k)
   rpt2_tot = pt2_tot / (1.d0 + pt1_norm)
 
 
   print*,'norm_ground_left_right_bi_orth = ',norm_ground_left_right_bi_orth(k)
   print*,'eigval_right_tc = ',eigval_right_tc_bi_orth(k)
   print*,'*****'
 
   if(print_pt2) then
     print*,'*****'
     print*,'previous wave function info'
     print*,'norm(before)      = ',norm
     print*,'E(before)         = ',E_tc 
     print*,'PT1 norm          = ',dsqrt(pt1_norm)
     print*,'PT2               = ',pt2_tot
     print*,'rPT2              = ',rpt2_tot
     print*,'|PT2|             = ',pt2_abs 
     print*,'Positive PT2      = ',pt2_plus
     print*,'Negative PT2      = ',pt2_minus
     print*,'E(before) + PT2   = ',E_tc + pt2_tot/norm
     print*,'E(before) +rPT2   = ',E_tc + rpt2_tot/norm
     write(*,'(A28,X,I10,X,100(F16.8,X))')'Ndet,E,E+PT2,E+RPT2,|PT2|=',ndet,E_tc ,E_tc  + pt2_tot/norm,E_tc  + rpt2_tot/norm,pt2_minus, pt2_plus
     print*,'*****'
   endif
   E_tc(k) = eigval_right_tc_bi_orth(k)
   norm(k) = norm_ground_left_right_bi_orth(k)
  enddo

  psi_energy(1:N_states) = eigval_right_tc_bi_orth(1:N_states) - nuclear_repulsion
  psi_s2(1:N_states) = s2_eigvec_tc_bi_orth(1:N_states)

  ndet = N_det
  do j = 1, N_states
    do i = 1, N_det
      psi_l_coef_bi_ortho(i,j) = leigvec_tc_bi_orth(i,j)
      psi_r_coef_bi_ortho(i,j) = reigvec_tc_bi_orth(i,j)
      psi_coef(i,j)            = dabs(psi_l_coef_bi_ortho(i,j) * psi_r_coef_bi_ortho(i,j))   
    enddo
  enddo
  SOFT_TOUCH eigval_left_tc_bi_orth eigval_right_tc_bi_orth leigvec_tc_bi_orth reigvec_tc_bi_orth norm_ground_left_right_bi_orth 
  SOFT_TOUCH psi_l_coef_bi_ortho psi_r_coef_bi_ortho psi_coef psi_energy psi_s2

  call save_tc_bi_ortho_wavefunction()

end

! ---

