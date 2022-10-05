program print_tc_bi_ortho
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
!  if(three_body_h_tc)then
!   call provide_all_three_ints_bi_ortho
!  endif
!  call routine
 call write_l_r_wf
end

subroutine write_l_r_wf
 implicit none
 character*(128) :: output
 integer :: i_unit_output,getUnitAndOpen                                                                                 
 output=trim(ezfio_filename)//'.tc_wf'
 i_unit_output = getUnitAndOpen(output,'w')
 integer :: i
 print*,'Writing the left-right wf'
 do i = 1, N_det
  write(i_unit_output,*)i,psi_l_coef_sorted_bi_ortho_left(i),psi_r_coef_sorted_bi_ortho_right(i)
 enddo


end

subroutine routine
 implicit none
 integer :: i,degree
 integer          :: exc(0:2,2,2),h1,p1,s1,h2,p2,s2
 double precision :: hmono,htwoe,hthree,htilde_ij,coef_pt1,e_i0,delta_e,e_pt2
 double precision :: contrib_pt,e_corr,coef,contrib,phase
 double precision :: accu_positive,accu_positive_pt, accu_positive_core,accu_positive_core_pt
 e_pt2 = 0.d0
 accu_positive = 0.D0
 accu_positive_pt = 0.D0
 accu_positive_core = 0.d0
 accu_positive_core_pt = 0.d0
 
 do i = 1, N_det
  call get_excitation_degree(HF_bitmask,psi_det(1,1,i),degree,N_int)
   if(degree == 1 .or. degree == 2)then
    call htilde_mu_mat_bi_ortho(psi_det(1,1,i),HF_bitmask,N_int,hmono,htwoe,hthree,htilde_ij)
    call htilde_mu_mat_bi_ortho(psi_det(1,1,i),psi_det(1,1,i),N_int,hmono,htwoe,hthree,e_i0)
    delta_e = e_tilde_00 - e_i0
    coef_pt1 = htilde_ij / delta_e
 
    call htilde_mu_mat_bi_ortho(HF_bitmask,psi_det(1,1,i),N_int,hmono,htwoe,hthree,htilde_ij)
    contrib_pt = coef_pt1 * htilde_ij
    e_pt2 += contrib_pt
 
    coef = psi_r_coef_bi_ortho(i,1)/psi_r_coef_bi_ortho(1,1)
    contrib = coef * htilde_ij
    e_corr += contrib
    call get_excitation(HF_bitmask,psi_det(1,1,i),exc,degree,phase,N_int)
    call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
    print*,'*********'
    if(degree==1)then
     print*,'s1',s1
     print*,'h1,p1 = ',h1,p1
    else if(degree ==2)then
     print*,'s1',s1
     print*,'h1,p1 = ',h1,p1
     print*,'s2',s2
     print*,'h2,p2 = ',h2,p2
    endif
    print*,'coef_pt1 = ',coef_pt1
    print*,'coef     = ',coef
    print*,'contrib_pt ',contrib_pt
    print*,'contrib  = ',contrib
    if(contrib.gt.0.d0)then
     accu_positive    += contrib
     if(h1==1.or.h2==1)then
      accu_positive_core += contrib
     endif
     if(dabs(contrib).gt.1.d-5)then
      print*,'Found a positive contribution to correlation energy !!'
     endif
    endif
    if(contrib_pt.gt.0.d0)then
     accu_positive_pt += contrib_pt
     if(h2==1.or.h1==1)then
      accu_positive_core_pt += contrib_pt
     endif
    endif
   endif
 enddo
 print*,''
 print*,''
 print*,'Total correlation energy            = ',e_corr
 print*,'Total correlation energy PT         = ',e_pt2
 print*,'Positive contribution to ecorr      = ',accu_positive
 print*,'Positive contribution to ecorr PT   = ',accu_positive_pt
 print*,'Pure core contribution              = ',accu_positive_core
 print*,'Pure core contribution PT           = ',accu_positive_core_pt
end
