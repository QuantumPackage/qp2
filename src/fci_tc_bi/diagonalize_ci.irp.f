
subroutine diagonalize_CI_tc_bi_ortho(ndet, E_tc,norm,pt2_data,print_pt2)
  use selection_types
  implicit none
  integer, intent(inout)              :: ndet      ! number of determinants from before 
  double precision, intent(inout)     :: E_tc,norm ! E and norm from previous wave function 
  type(pt2_type)  , intent(in)        :: pt2_data  ! PT2 from previous wave function 
  logical, intent(in) :: print_pt2
  BEGIN_DOC
!  Replace the coefficients of the CI states by the coefficients of the
!  eigenstates of the CI matrix
  END_DOC
  integer :: i,j
 double precision :: pt2_tmp,pt1_norm,rpt2_tmp,abs_pt2
 pt2_tmp = pt2_data % pt2(1)
 abs_pt2 = pt2_data % variance(1)
 pt1_norm = pt2_data % overlap(1,1)
 rpt2_tmp = pt2_tmp/(1.d0 + pt1_norm)
  print*,'*****'
  print*,'New wave function information'
  print*,'N_det tc               = ',N_det
  print*,'norm_ground_left_right_bi_orth = ',norm_ground_left_right_bi_orth
  print*,'eigval_right_tc = ',eigval_right_tc_bi_orth(1)
  print*,'Ndet, E_tc = ',N_det,eigval_right_tc_bi_orth(1)
  print*,'*****'
  if(print_pt2)then
   print*,'*****'
   print*,'previous wave function info'
   print*,'norm(before)      = ',norm
   print*,'E(before)         = ',E_tc 
   print*,'PT1 norm          = ',dsqrt(pt1_norm)
   print*,'PT2               = ',pt2_tmp
   print*,'rPT2              = ',rpt2_tmp
   print*,'|PT2|             = ',abs_pt2
   print*,'Positive PT2      = ',(pt2_tmp + abs_pt2)*0.5d0
   print*,'Negative PT2      = ',(pt2_tmp - abs_pt2)*0.5d0
   print*,'E(before) + PT2   = ',E_tc + pt2_tmp/norm
   print*,'E(before) +rPT2   = ',E_tc + rpt2_tmp/norm
   write(*,'(A28,X,I10,X,100(F16.8,X))')'Ndet,E,E+PT2,E+RPT2,|PT2|=',ndet,E_tc ,E_tc  + pt2_tmp/norm,E_tc  + rpt2_tmp/norm,abs_pt2
   print*,'*****'
  endif
  E_tc  = eigval_right_tc_bi_orth(1)
  norm  = norm_ground_left_right_bi_orth
  ndet  = N_det
  do j=1,N_states
    do i=1,N_det
      psi_l_coef_bi_ortho(i,j) = leigvec_tc_bi_orth(i,j)
      psi_r_coef_bi_ortho(i,j) = reigvec_tc_bi_orth(i,j)
      psi_coef(i,j) = dabs(psi_l_coef_bi_ortho(i,j) * psi_r_coef_bi_ortho(i,j))   
    enddo
  enddo
  SOFT_TOUCH  eigval_left_tc_bi_orth  eigval_right_tc_bi_orth  leigvec_tc_bi_orth  reigvec_tc_bi_orth norm_ground_left_right_bi_orth psi_coef  psi_l_coef_bi_ortho psi_r_coef_bi_ortho



  call save_tc_bi_ortho_wavefunction
end

subroutine print_CI_dressed(ndet, E_tc,norm,pt2_data,print_pt2)
  use selection_types
  implicit none
  integer, intent(inout)              :: ndet      ! number of determinants from before 
  double precision, intent(inout)     :: E_tc,norm ! E and norm from previous wave function 
  type(pt2_type)  , intent(in)        :: pt2_data  ! PT2 from previous wave function 
  logical, intent(in) :: print_pt2
  BEGIN_DOC
!  Replace the coefficients of the CI states by the coefficients of the
!  eigenstates of the CI matrix
  END_DOC
  integer :: i,j
  print*,'*****'
  print*,'New wave function information'
  print*,'N_det tc               = ',N_det
  print*,'norm_ground_left_right_bi_orth = ',norm_ground_left_right_bi_orth
  print*,'eigval_right_tc = ',eigval_right_tc_bi_orth(1)
  print*,'Ndet, E_tc = ',N_det,eigval_right_tc_bi_orth(1)
  print*,'*****'
  if(print_pt2)then
   print*,'*****'
   print*,'previous wave function info'
   print*,'norm(before)      = ',norm
   print*,'E(before)         = ',E_tc
   print*,'PT1 norm          = ',dsqrt(pt2_data % overlap(1,1))
   print*,'E(before) + PT2   = ',E_tc + (pt2_data % pt2(1))/norm
   print*,'PT2               = ',pt2_data % pt2(1)
   print*,'Ndet, E_tc, E+PT2 = ',ndet,E_tc,E_tc + (pt2_data % pt2(1))/norm,dsqrt(pt2_data % overlap(1,1))
   print*,'*****'
  endif
  E_tc  = eigval_right_tc_bi_orth(1)
  norm  = norm_ground_left_right_bi_orth
  ndet  = N_det
  do j=1,N_states
    do i=1,N_det
      psi_coef(i,j) = reigvec_tc_bi_orth(i,j)
    enddo
  enddo
  SOFT_TOUCH  eigval_left_tc_bi_orth  eigval_right_tc_bi_orth  leigvec_tc_bi_orth  norm_ground_left_right_bi_orth  psi_coef  reigvec_tc_bi_orth 

end

