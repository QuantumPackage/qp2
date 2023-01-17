program test_ints

  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC

  implicit none

  print *, ' starting test_ints ...'

  my_grid_becke  = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
!  my_n_pt_r_grid = 15 ! small grid for quick debug
!  my_n_pt_a_grid = 26 ! small grid for quick debug
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  my_extra_grid_becke = .True.
  my_n_pt_r_extra_grid = 30
  my_n_pt_a_extra_grid = 50 ! small extra_grid for quick debug
  touch my_extra_grid_becke my_n_pt_r_extra_grid my_n_pt_a_extra_grid

!! OK 
!call routine_int2_u_grad1u_j1b2 
!! OK
!call routine_v_ij_erf_rk_cst_mu_j1b
!! OK 
! call routine_x_v_ij_erf_rk_cst_mu_j1b
!! OK
! call routine_v_ij_u_cst_mu_j1b

!! OK
!call routine_int2_u2_j1b2

!! OK
!call routine_int2_u_grad1u_x_j1b2

!! OK 
! call routine_int2_grad1u2_grad2u2_j1b2
! call routine_int2_u_grad1u_j1b2
! call test_total_grad_lapl
! call test_total_grad_square
! call test_ao_tc_int_chemist
! call test_grid_points_ao
! call test_tc_scf
 !call test_int_gauss

  !call test_fock_3e_uhf_ao()
  !call test_fock_3e_uhf_mo()

  !call test_tc_grad_and_lapl_ao()
  !call test_tc_grad_square_ao()

  call test_two_e_tc_non_hermit_integral()

end

! ---

subroutine test_tc_scf
 implicit none
 integer :: i
! provide int2_u_grad1u_x_j1b2_test
 provide x_v_ij_erf_rk_cst_mu_j1b_test
! provide x_v_ij_erf_rk_cst_mu_j1b_test
! print*,'TC_HF_energy = ',TC_HF_energy
! print*,'grad_non_hermit = ',grad_non_hermit
end

subroutine test_ao_tc_int_chemist
 implicit none
 provide ao_tc_int_chemist
! provide ao_tc_int_chemist_test
! provide tc_grad_square_ao_test
! provide tc_grad_and_lapl_ao_test
end

! ---

subroutine routine_test_j1b
 implicit none
 integer :: i,icount,j
 icount = 0
 do i = 1, List_all_comb_b3_size
  if(dabs(List_all_comb_b3_coef(i)).gt.1.d-10)then
   print*,''
   print*,List_all_comb_b3_expo(i),List_all_comb_b3_coef(i)
   print*,List_all_comb_b3_cent(1:3,i)
   print*,''
   icount += 1
  endif
  
 enddo
 print*,'List_all_comb_b3_coef,icount = ',List_all_comb_b3_size,icount
 do i = 1, ao_num
  do j = 1, ao_num
   do icount = 1, List_comb_thr_b3_size(j,i)
    print*,'',j,i
    print*,List_comb_thr_b3_expo(icount,j,i),List_comb_thr_b3_coef(icount,j,i)
    print*,List_comb_thr_b3_cent(1:3,icount,j,i)
    print*,''
   enddo
!   enddo
  enddo
 enddo
 print*,'max_List_comb_thr_b3_size = ',max_List_comb_thr_b3_size,List_all_comb_b3_size

end

subroutine routine_int2_u_grad1u_j1b2
 implicit none
 integer :: i,j,ipoint,k,l
 double precision :: weight,accu_relat, accu_abs, contrib
 double precision, allocatable :: array(:,:,:,:), array_ref(:,:,:,:)

 allocate(array(ao_num, ao_num, ao_num, ao_num))
 array = 0.d0
 allocate(array_ref(ao_num, ao_num, ao_num, ao_num))
 array_ref = 0.d0
 do ipoint = 1, n_points_final_grid
  weight = final_weight_at_r_vector(ipoint)
  do k = 1, ao_num
   do l = 1, ao_num
    do i = 1, ao_num
     do j = 1, ao_num
      array(j,i,l,k)     += int2_u_grad1u_j1b2_test(j,i,ipoint) * aos_in_r_array(k,ipoint) * aos_in_r_array(l,ipoint) * weight
      array_ref(j,i,l,k) += int2_u_grad1u_j1b2(j,i,ipoint)      * aos_in_r_array(k,ipoint) * aos_in_r_array(l,ipoint) * weight
     enddo
    enddo
   enddo
  enddo
 enddo
 accu_relat = 0.d0
 accu_abs   = 0.d0
  do k = 1, ao_num
   do l = 1, ao_num
    do i = 1, ao_num
     do j = 1, ao_num
      contrib = dabs(array(j,i,l,k) - array_ref(j,i,l,k))
      accu_abs += contrib
      if(dabs(array_ref(j,i,l,k)).gt.1.d-10)then
       accu_relat += contrib/dabs(array_ref(j,i,l,k))
      endif
     enddo
    enddo
   enddo
  enddo
 print*,'accu_abs   = ',accu_abs/dble(ao_num)**4
 print*,'accu_relat = ',accu_relat/dble(ao_num)**4

  

end

subroutine routine_v_ij_erf_rk_cst_mu_j1b
 implicit none
 integer :: i,j,ipoint,k,l
 double precision :: weight,accu_relat, accu_abs, contrib
 double precision, allocatable :: array(:,:,:,:), array_ref(:,:,:,:)
! print*,'ao_overlap_abs = '
! do i = 1, ao_num
!   write(*,'(100(F10.5,X))')ao_overlap_abs(i,:)
! enddo
! print*,'center = '
! do i = 1, ao_num
!   write(*,'(100(F10.5,X))')ao_prod_center(2,i,:)
! enddo
! print*,'sigma = '
! do i = 1, ao_num
!   write(*,'(100(F10.5,X))')ao_prod_sigma(i,:)
! enddo


 allocate(array(ao_num, ao_num, ao_num, ao_num))
 array = 0.d0
 allocate(array_ref(ao_num, ao_num, ao_num, ao_num))
 array_ref = 0.d0
 do ipoint = 1, n_points_final_grid
  weight = final_weight_at_r_vector(ipoint)
  do k = 1, ao_num
   do l = 1, ao_num
    do i = 1, ao_num
     do j = 1, ao_num
      array(j,i,l,k)     += v_ij_erf_rk_cst_mu_j1b_test(j,i,ipoint) * aos_in_r_array(k,ipoint) * aos_in_r_array(l,ipoint) * weight
      array_ref(j,i,l,k) += v_ij_erf_rk_cst_mu_j1b(j,i,ipoint)      * aos_in_r_array(k,ipoint) * aos_in_r_array(l,ipoint) * weight
     enddo
    enddo
   enddo
  enddo
 enddo
 accu_relat = 0.d0
 accu_abs   = 0.d0
  do k = 1, ao_num
   do l = 1, ao_num
    do i = 1, ao_num
     do j = 1, ao_num
      contrib = dabs(array(j,i,l,k) - array_ref(j,i,l,k))
      accu_abs += contrib
      if(dabs(array_ref(j,i,l,k)).gt.1.d-10)then
       accu_relat += contrib/dabs(array_ref(j,i,l,k))
      endif
     enddo
    enddo
   enddo
  enddo
 print*,'accu_abs   = ',accu_abs/dble(ao_num)**4
 print*,'accu_relat = ',accu_relat/dble(ao_num)**4

  

end


subroutine routine_x_v_ij_erf_rk_cst_mu_j1b
 implicit none
 integer :: i,j,ipoint,k,l,m
 double precision :: weight,accu_relat, accu_abs, contrib
 double precision, allocatable :: array(:,:,:,:), array_ref(:,:,:,:)
! print*,'ao_overlap_abs = '
! do i = 1, ao_num
!   write(*,'(100(F10.5,X))')ao_overlap_abs(i,:)
! enddo
! print*,'center = '
! do i = 1, ao_num
!   write(*,'(100(F10.5,X))')ao_prod_center(2,i,:)
! enddo
! print*,'sigma = '
! do i = 1, ao_num
!   write(*,'(100(F10.5,X))')ao_prod_sigma(i,:)
! enddo


 allocate(array(ao_num, ao_num, ao_num, ao_num))
 array = 0.d0
 allocate(array_ref(ao_num, ao_num, ao_num, ao_num))
 array_ref = 0.d0
 do ipoint = 1, n_points_final_grid
  weight = final_weight_at_r_vector(ipoint)
  do k = 1, ao_num
   do l = 1, ao_num
    do i = 1, ao_num
     do j = 1, ao_num
      do m = 1, 3
       array(j,i,l,k)     += x_v_ij_erf_rk_cst_mu_j1b_test(j,i,ipoint,m) * aos_grad_in_r_array_transp(m,k,ipoint) * aos_in_r_array(l,ipoint) * weight
       array_ref(j,i,l,k) += x_v_ij_erf_rk_cst_mu_j1b     (j,i,ipoint,m) * aos_grad_in_r_array_transp(m,k,ipoint) * aos_in_r_array(l,ipoint) * weight
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
 accu_relat = 0.d0
 accu_abs   = 0.d0
  do k = 1, ao_num
   do l = 1, ao_num
    do i = 1, ao_num
     do j = 1, ao_num
      contrib = dabs(array(j,i,l,k) - array_ref(j,i,l,k))
      accu_abs += contrib
      if(dabs(array_ref(j,i,l,k)).gt.1.d-10)then
       accu_relat += contrib/dabs(array_ref(j,i,l,k))
      endif
     enddo
    enddo
   enddo
  enddo
 print*,'accu_abs   = ',accu_abs/dble(ao_num)**4
 print*,'accu_relat = ',accu_relat/dble(ao_num)**4

  

end



subroutine routine_v_ij_u_cst_mu_j1b_test
 implicit none
 integer :: i,j,ipoint,k,l
 double precision :: weight,accu_relat, accu_abs, contrib
 double precision, allocatable :: array(:,:,:,:), array_ref(:,:,:,:)
! print*,'ao_overlap_abs = '
! do i = 1, ao_num
!   write(*,'(100(F10.5,X))')ao_overlap_abs(i,:)
! enddo
! print*,'center = '
! do i = 1, ao_num
!   write(*,'(100(F10.5,X))')ao_prod_center(2,i,:)
! enddo
! print*,'sigma = '
! do i = 1, ao_num
!   write(*,'(100(F10.5,X))')ao_prod_sigma(i,:)
! enddo


 allocate(array(ao_num, ao_num, ao_num, ao_num))
 array = 0.d0
 allocate(array_ref(ao_num, ao_num, ao_num, ao_num))
 array_ref = 0.d0
 do ipoint = 1, n_points_final_grid
  weight = final_weight_at_r_vector(ipoint)
  do k = 1, ao_num
   do l = 1, ao_num
    do i = 1, ao_num
     do j = 1, ao_num
      array(j,i,l,k)     += v_ij_u_cst_mu_j1b_test(j,i,ipoint) * aos_in_r_array(k,ipoint) * aos_in_r_array(l,ipoint) * weight
      array_ref(j,i,l,k) += v_ij_u_cst_mu_j1b(j,i,ipoint)      * aos_in_r_array(k,ipoint) * aos_in_r_array(l,ipoint) * weight
     enddo
    enddo
   enddo
  enddo
 enddo
 accu_relat = 0.d0
 accu_abs   = 0.d0
  do k = 1, ao_num
   do l = 1, ao_num
    do i = 1, ao_num
     do j = 1, ao_num
      contrib = dabs(array(j,i,l,k) - array_ref(j,i,l,k))
      accu_abs += contrib
      if(dabs(array_ref(j,i,l,k)).gt.1.d-10)then
       accu_relat += contrib/dabs(array_ref(j,i,l,k))
      endif
     enddo
    enddo
   enddo
  enddo
 print*,'accu_abs   = ',accu_abs/dble(ao_num)**4
 print*,'accu_relat = ',accu_relat/dble(ao_num)**4

  

end

subroutine routine_int2_grad1u2_grad2u2_j1b2
 implicit none
 integer :: i,j,ipoint,k,l
 integer :: ii , jj
 double precision :: weight,accu_relat, accu_abs, contrib
 double precision, allocatable :: array(:,:,:,:), array_ref(:,:,:,:)
 double precision, allocatable :: ints(:,:,:)
 allocate(ints(ao_num, ao_num, n_points_final_grid))
! do ipoint = 1, n_points_final_grid
!  do i = 1, ao_num
!   do j = 1, ao_num
!    read(33,*)ints(j,i,ipoint)
!   enddo
!  enddo
! enddo

 allocate(array(ao_num, ao_num, ao_num, ao_num))
 array = 0.d0
 allocate(array_ref(ao_num, ao_num, ao_num, ao_num))
 array_ref = 0.d0
 do ipoint = 1, n_points_final_grid
  weight = final_weight_at_r_vector(ipoint)
  do k = 1, ao_num
   do l = 1, ao_num
    do i = 1, ao_num
     do j = 1, ao_num
      array(j,i,l,k)     += int2_grad1u2_grad2u2_j1b2_test(j,i,ipoint) * aos_in_r_array(k,ipoint) * aos_in_r_array(l,ipoint) * weight
!     !array(j,i,l,k)     += int2_grad1u2_grad2u2_j1b2_test(j,i,ipoint) * aos_in_r_array(k,ipoint) * aos_in_r_array(l,ipoint) * weight
!      array_ref(j,i,l,k)     += int2_grad1u2_grad2u2_j1b2_test(j,i,ipoint) * aos_in_r_array(k,ipoint) * aos_in_r_array(l,ipoint) * weight
!     !array(j,i,l,k) += ints(j,i,ipoint)      * aos_in_r_array(k,ipoint) * aos_in_r_array(l,ipoint) * weight
!       array_ref(j,i,l,k) += int2_grad1u2_grad2u2_j1b2(j,i,ipoint)      * aos_in_r_array(k,ipoint) * aos_in_r_array(l,ipoint) * weight
      array_ref(j,i,l,k) += ints(j,i,ipoint)      * aos_in_r_array(k,ipoint) * aos_in_r_array(l,ipoint) * weight
!      if(dabs(int2_grad1u2_grad2u2_j1b2_test(j,i,ipoint)).gt.1.d-6)then
!       if(dabs(int2_grad1u2_grad2u2_j1b2_test(j,i,ipoint) - int2_grad1u2_grad2u2_j1b2_test(j,i,ipoint)).gt.1.d-6)then
!        print*,j,i,ipoint
!        print*,int2_grad1u2_grad2u2_j1b2_test(j,i,ipoint) , int2_grad1u2_grad2u2_j1b2_test(j,i,ipoint), dabs(int2_grad1u2_grad2u2_j1b2_test(j,i,ipoint) - int2_grad1u2_grad2u2_j1b2_test(j,i,ipoint))
!        print*,int2_grad1u2_grad2u2_j1b2_test(i,j,ipoint) , int2_grad1u2_grad2u2_j1b2_test(i,j,ipoint), dabs(int2_grad1u2_grad2u2_j1b2_test(i,j,ipoint) - int2_grad1u2_grad2u2_j1b2_test(i,j,ipoint))
!        stop
!       endif
!      endif
     enddo
    enddo
   enddo
  enddo
 enddo
 double precision :: e_ref, e_new
 accu_relat = 0.d0
 accu_abs   = 0.d0
 e_ref = 0.d0
 e_new = 0.d0
 do ii = 1, elec_alpha_num
  do jj = ii, elec_alpha_num
   do k = 1, ao_num
    do l = 1, ao_num
     do i = 1, ao_num
      do j = 1, ao_num
       e_ref += mo_coef(j,ii) * mo_coef(i,ii) * array_ref(j,i,l,k) * mo_coef(l,jj) * mo_coef(k,jj)
       e_new += mo_coef(j,ii) * mo_coef(i,ii) * array(j,i,l,k) * mo_coef(l,jj) * mo_coef(k,jj)
       contrib = dabs(array(j,i,l,k) - array_ref(j,i,l,k))
       accu_abs += contrib
!       if(dabs(array_ref(j,i,l,k)).gt.1.d-10)then
!        accu_relat += contrib/dabs(array_ref(j,i,l,k))
!       endif
      enddo
     enddo
    enddo
   enddo

  enddo
 enddo
 print*,'e_ref = ',e_ref
 print*,'e_new = ',e_new
! print*,'accu_abs   = ',accu_abs/dble(ao_num)**4
! print*,'accu_relat = ',accu_relat/dble(ao_num)**4

  

end

subroutine routine_int2_u2_j1b2
 implicit none
 integer :: i,j,ipoint,k,l
 double precision :: weight,accu_relat, accu_abs, contrib
 double precision, allocatable :: array(:,:,:,:), array_ref(:,:,:,:)
! print*,'ao_overlap_abs = '
! do i = 1, ao_num
!   write(*,'(100(F10.5,X))')ao_overlap_abs(i,:)
! enddo
! print*,'center = '
! do i = 1, ao_num
!   write(*,'(100(F10.5,X))')ao_prod_center(2,i,:)
! enddo
! print*,'sigma = '
! do i = 1, ao_num
!   write(*,'(100(F10.5,X))')ao_prod_sigma(i,:)
! enddo


 allocate(array(ao_num, ao_num, ao_num, ao_num))
 array = 0.d0
 allocate(array_ref(ao_num, ao_num, ao_num, ao_num))
 array_ref = 0.d0
 do ipoint = 1, n_points_final_grid
  weight = final_weight_at_r_vector(ipoint)
  do k = 1, ao_num
   do l = 1, ao_num
    do i = 1, ao_num
     do j = 1, ao_num
      array(j,i,l,k)     += int2_u2_j1b2_test(j,i,ipoint) * aos_in_r_array(k,ipoint) * aos_in_r_array(l,ipoint) * weight
      array_ref(j,i,l,k) += int2_u2_j1b2(j,i,ipoint)      * aos_in_r_array(k,ipoint) * aos_in_r_array(l,ipoint) * weight
     enddo
    enddo
   enddo
  enddo
 enddo
 accu_relat = 0.d0
 accu_abs   = 0.d0
  do k = 1, ao_num
   do l = 1, ao_num
    do i = 1, ao_num
     do j = 1, ao_num
      contrib = dabs(array(j,i,l,k) - array_ref(j,i,l,k))
      accu_abs += contrib
      if(dabs(array_ref(j,i,l,k)).gt.1.d-10)then
       accu_relat += contrib/dabs(array_ref(j,i,l,k))
      endif
     enddo
    enddo
   enddo
  enddo
 print*,'accu_abs   = ',accu_abs/dble(ao_num)**4
 print*,'accu_relat = ',accu_relat/dble(ao_num)**4

  

end


subroutine routine_int2_u_grad1u_x_j1b2
 implicit none
 integer :: i,j,ipoint,k,l,m
 double precision :: weight,accu_relat, accu_abs, contrib
 double precision, allocatable :: array(:,:,:,:), array_ref(:,:,:,:)
! print*,'ao_overlap_abs = '
! do i = 1, ao_num
!   write(*,'(100(F10.5,X))')ao_overlap_abs(i,:)
! enddo
! print*,'center = '
! do i = 1, ao_num
!   write(*,'(100(F10.5,X))')ao_prod_center(2,i,:)
! enddo
! print*,'sigma = '
! do i = 1, ao_num
!   write(*,'(100(F10.5,X))')ao_prod_sigma(i,:)
! enddo


 allocate(array(ao_num, ao_num, ao_num, ao_num))
 array = 0.d0
 allocate(array_ref(ao_num, ao_num, ao_num, ao_num))
 array_ref = 0.d0
 do ipoint = 1, n_points_final_grid
  weight = final_weight_at_r_vector(ipoint)
  do k = 1, ao_num
   do l = 1, ao_num
    do i = 1, ao_num
     do j = 1, ao_num
      do m = 1, 3
       array(j,i,l,k)     += int2_u_grad1u_x_j1b2_test(j,i,ipoint,m) * aos_grad_in_r_array_transp(m,k,ipoint) * aos_in_r_array(l,ipoint) * weight
       array_ref(j,i,l,k) += int2_u_grad1u_x_j1b2     (j,i,ipoint,m) * aos_grad_in_r_array_transp(m,k,ipoint) * aos_in_r_array(l,ipoint) * weight
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
 accu_relat = 0.d0
 accu_abs   = 0.d0
  do k = 1, ao_num
   do l = 1, ao_num
    do i = 1, ao_num
     do j = 1, ao_num
      contrib = dabs(array(j,i,l,k) - array_ref(j,i,l,k))
      accu_abs += contrib
      if(dabs(array_ref(j,i,l,k)).gt.1.d-10)then
       accu_relat += contrib/dabs(array_ref(j,i,l,k))
      endif
     enddo
    enddo
   enddo
  enddo
 print*,'accu_abs   = ',accu_abs/dble(ao_num)**4
 print*,'accu_relat = ',accu_relat/dble(ao_num)**4

  

end

subroutine routine_v_ij_u_cst_mu_j1b
 implicit none
 integer :: i,j,ipoint,k,l
 double precision :: weight,accu_relat, accu_abs, contrib
 double precision, allocatable :: array(:,:,:,:), array_ref(:,:,:,:)

 allocate(array(ao_num, ao_num, ao_num, ao_num))
 array = 0.d0
 allocate(array_ref(ao_num, ao_num, ao_num, ao_num))
 array_ref = 0.d0
 do ipoint = 1, n_points_final_grid
  weight = final_weight_at_r_vector(ipoint)
  do k = 1, ao_num
   do l = 1, ao_num
    do i = 1, ao_num
     do j = 1, ao_num
      array(j,i,l,k)     += v_ij_u_cst_mu_j1b_test(j,i,ipoint) * aos_in_r_array(k,ipoint) * aos_in_r_array(l,ipoint) * weight
      array_ref(j,i,l,k) += v_ij_u_cst_mu_j1b(j,i,ipoint)      * aos_in_r_array(k,ipoint) * aos_in_r_array(l,ipoint) * weight
     enddo
    enddo
   enddo
  enddo
 enddo
 accu_relat = 0.d0
 accu_abs   = 0.d0
  do k = 1, ao_num
   do l = 1, ao_num
    do i = 1, ao_num
     do j = 1, ao_num
      contrib = dabs(array(j,i,l,k) - array_ref(j,i,l,k))
      accu_abs += contrib
      if(dabs(array_ref(j,i,l,k)).gt.1.d-10)then
       accu_relat += contrib/dabs(array_ref(j,i,l,k))
      endif
     enddo
    enddo
   enddo
  enddo
 print*,'accu_abs   = ',accu_abs/dble(ao_num)**4
 print*,'accu_relat = ',accu_relat/dble(ao_num)**4

end

! ---

subroutine test_fock_3e_uhf_ao()

  implicit none
  integer                       :: i, j
  double precision              :: diff_tot, diff_ij, thr_ih, norm
  double precision, allocatable :: fock_3e_uhf_ao_a_mo(:,:), fock_3e_uhf_ao_b_mo(:,:)

  thr_ih = 1d-7

  PROVIDE fock_a_tot_3e_bi_orth fock_b_tot_3e_bi_orth
  PROVIDE fock_3e_uhf_ao_a fock_3e_uhf_ao_b

  ! ---

  allocate(fock_3e_uhf_ao_a_mo(mo_num,mo_num))
  call ao_to_mo_bi_ortho( fock_3e_uhf_ao_a   , size(fock_3e_uhf_ao_a   , 1) &
                        , fock_3e_uhf_ao_a_mo, size(fock_3e_uhf_ao_a_mo, 1) )

  norm     = 0.d0
  diff_tot = 0.d0
  do i = 1, mo_num
    do j = 1, mo_num

      diff_ij = dabs(fock_3e_uhf_ao_a_mo(j,i) - fock_a_tot_3e_bi_orth(j,i))
      if(diff_ij .gt. thr_ih) then
        print *, ' difference on ', j, i
        print *, ' MANU : ', fock_a_tot_3e_bi_orth(j,i)
        print *, ' UHF  : ', fock_3e_uhf_ao_a_mo  (j,i)
        !stop
      endif

      norm     += dabs(fock_a_tot_3e_bi_orth(j,i))
      diff_tot += diff_ij
    enddo
  enddo
  print *, ' diff on F_a = ', diff_tot / norm
  print *, ' '

  deallocate(fock_3e_uhf_ao_a_mo)

  ! ---

  allocate(fock_3e_uhf_ao_b_mo(mo_num,mo_num))
  call ao_to_mo_bi_ortho( fock_3e_uhf_ao_b   , size(fock_3e_uhf_ao_b   , 1) &
                        , fock_3e_uhf_ao_b_mo, size(fock_3e_uhf_ao_b_mo, 1) )

  norm     = 0.d0
  diff_tot = 0.d0
  do i = 1, mo_num
    do j = 1, mo_num

      diff_ij = dabs(fock_3e_uhf_ao_b_mo(j,i) - fock_b_tot_3e_bi_orth(j,i))
      if(diff_ij .gt. thr_ih) then
        print *, ' difference on ', j, i
        print *, ' MANU : ', fock_b_tot_3e_bi_orth(j,i)
        print *, ' UHF  : ', fock_3e_uhf_ao_b_mo  (j,i)
        !stop
      endif

      norm     += dabs(fock_b_tot_3e_bi_orth(j,i))
      diff_tot += diff_ij
    enddo
  enddo
  print *, ' diff on F_b = ', diff_tot/norm
  print *, ' '

  deallocate(fock_3e_uhf_ao_b_mo)

  ! ---

end subroutine test_fock_3e_uhf_ao()

! ---

subroutine test_fock_3e_uhf_mo()

  implicit none
  integer          :: i, j
  double precision :: diff_tot, diff_ij, thr_ih, norm

  thr_ih = 1d-12

  PROVIDE fock_a_tot_3e_bi_orth fock_b_tot_3e_bi_orth
  PROVIDE fock_3e_uhf_mo_a fock_3e_uhf_mo_b

  ! ---

  norm     = 0.d0
  diff_tot = 0.d0
  do i = 1, mo_num
    do j = 1, mo_num

      diff_ij = dabs(fock_3e_uhf_mo_a(j,i) - fock_a_tot_3e_bi_orth(j,i))
      if(diff_ij .gt. thr_ih) then
        print *, ' difference on ', j, i
        print *, ' MANU : ', fock_a_tot_3e_bi_orth(j,i)
        print *, ' UHF  : ', fock_3e_uhf_mo_a     (j,i)
        !stop
      endif

      norm     += dabs(fock_a_tot_3e_bi_orth(j,i))
      diff_tot += diff_ij
    enddo
  enddo
  print *, ' diff on F_a = ', diff_tot / norm
  print *, '      norm_a = ', norm
  print *, ' '

  ! ---

  norm     = 0.d0
  diff_tot = 0.d0
  do i = 1, mo_num
    do j = 1, mo_num

      diff_ij = dabs(fock_3e_uhf_mo_b(j,i) - fock_b_tot_3e_bi_orth(j,i))
      if(diff_ij .gt. thr_ih) then
        print *, ' difference on ', j, i
        print *, ' MANU : ', fock_b_tot_3e_bi_orth(j,i)
        print *, ' UHF  : ', fock_3e_uhf_mo_b     (j,i)
        !stop
      endif

      norm     += dabs(fock_b_tot_3e_bi_orth(j,i))
      diff_tot += diff_ij
    enddo
  enddo
  print *, ' diff on F_b = ', diff_tot/norm
  print *, '      norm_b = ', norm
  print *, ' '

  ! ---

end subroutine test_fock_3e_uhf_mo

! ---

subroutine test_total_grad_lapl
 implicit none
 integer :: i,j,ipoint,k,l
 double precision :: weight,accu_relat, accu_abs, contrib
 accu_relat = 0.d0
 accu_abs   = 0.d0
  do k = 1, ao_num
   do l = 1, ao_num
    do i = 1, ao_num
     do j = 1, ao_num
      contrib = dabs(tc_grad_and_lapl_ao_test(j,i,l,k) - tc_grad_and_lapl_ao(j,i,l,k))
      accu_abs += contrib
      if(dabs(tc_grad_and_lapl_ao(j,i,l,k)).gt.1.d-10)then
       accu_relat += contrib/dabs(tc_grad_and_lapl_ao(j,i,l,k))
      endif
     enddo
    enddo
   enddo
  enddo
 print*,'accu_abs   = ',accu_abs/dble(ao_num)**4
 print*,'accu_relat = ',accu_relat/dble(ao_num)**4


end

subroutine test_total_grad_square
 implicit none
 integer :: i,j,ipoint,k,l
 double precision :: weight,accu_relat, accu_abs, contrib
 accu_relat = 0.d0
 accu_abs   = 0.d0
  do k = 1, ao_num
   do l = 1, ao_num
    do i = 1, ao_num
     do j = 1, ao_num
      contrib = dabs(tc_grad_square_ao_test(j,i,l,k) - tc_grad_square_ao(j,i,l,k))
      accu_abs += contrib
      if(dabs(tc_grad_square_ao(j,i,l,k)).gt.1.d-10)then
       accu_relat += contrib/dabs(tc_grad_square_ao(j,i,l,k))
      endif
     enddo
    enddo
   enddo
  enddo
 print*,'accu_abs   = ',accu_abs/dble(ao_num)**4
 print*,'accu_relat = ',accu_relat/dble(ao_num)**4


end

subroutine test_grid_points_ao
 implicit none
 integer :: i,j,ipoint,icount,icount_good, icount_bad,icount_full
 double precision :: thr
 thr = 1.d-10
! print*,'max_n_pts_grid_ao_prod = ',max_n_pts_grid_ao_prod
! print*,'n_pts_grid_ao_prod'
 do i = 1, ao_num
  do j = i, ao_num
  icount = 0
  icount_good = 0
  icount_bad = 0
  icount_full = 0
  do ipoint = 1, n_points_final_grid
!   if(dabs(int2_u_grad1u_x_j1b2_test(j,i,ipoint,1)) & 
!    + dabs(int2_u_grad1u_x_j1b2_test(j,i,ipoint,2)) &
!    + dabs(int2_u_grad1u_x_j1b2_test(j,i,ipoint,3)) )
!   if(dabs(int2_u2_j1b2_test(j,i,ipoint)).gt.thr)then
!    icount += 1
!   endif
   if(dabs(v_ij_u_cst_mu_j1b_ng_1_test(j,i,ipoint)).gt.thr*0.1d0)then
    icount_full += 1
   endif
   if(dabs(v_ij_u_cst_mu_j1b_test(j,i,ipoint)).gt.thr)then
    icount += 1
    if(dabs(v_ij_u_cst_mu_j1b_ng_1_test(j,i,ipoint)).gt.thr*0.1d0)then
    icount_good += 1
    else
    print*,j,i,ipoint
    print*,dabs(v_ij_u_cst_mu_j1b_test(j,i,ipoint)),dabs(v_ij_u_cst_mu_j1b_ng_1_test(j,i,ipoint)),dabs(v_ij_u_cst_mu_j1b_ng_1_test(j,i,ipoint))/dabs(v_ij_u_cst_mu_j1b_test(j,i,ipoint))
    icount_bad  += 1
    endif
   endif
!   if(dabs(v_ij_u_cst_mu_j1b_ng_1_test(j,i,ipoint)).gt.thr)then
!   endif
  enddo
   print*,''
   print*,j,i
   print*,icount,icount_full, icount_bad!,n_pts_grid_ao_prod(j,i)
   print*,dble(icount)/dble(n_points_final_grid),dble(icount_full)/dble(n_points_final_grid)
!          dble(n_pts_grid_ao_prod(j,i))/dble(n_points_final_grid)
!   if(icount.gt.n_pts_grid_ao_prod(j,i))then
!    print*,'pb !!'
!   endif
  enddo
 enddo
end

subroutine test_int_gauss
 implicit none
 integer :: i,j
 print*,'center'
 do i = 1, ao_num
  do j = i, ao_num
   print*,j,i
   print*,ao_prod_sigma(j,i),ao_overlap_abs_grid(j,i)
   print*,ao_prod_center(1:3,j,i)
  enddo
 enddo
 print*,''
 double precision :: weight, r(3),integral_1,pi,center(3),f_r,alpha,distance,integral_2
 center = 0.d0
 pi = dacos(-1.d0)
 integral_1 = 0.d0
 integral_2 = 0.d0
 alpha = 0.75d0
 do i = 1,  n_points_final_grid
  ! you get x, y and z of the ith grid point
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
  weight = final_weight_at_r_vector(i)
  distance = dsqrt( (r(1) - center(1))**2 +  (r(2) - center(2))**2 + (r(3) - center(3))**2 )
  f_r = dexp(-alpha * distance*distance)
  ! you add the contribution of the grid point to the integral
  integral_1 += f_r * weight
  integral_2 += f_r * distance * weight
 enddo
 print*,'integral_1      =',integral_1
 print*,'(pi/alpha)**1.5 =',(pi / alpha)**1.5
 print*,'integral_2      =',integral_2
 print*,'(pi/alpha)**1.5 =',2.d0*pi / (alpha)**2


end

! ---

subroutine test_tc_grad_and_lapl_ao()

  implicit none
  integer          :: i, j, k, l
  double precision :: diff_tot, diff, thr_ih, norm

  thr_ih = 1d-10

  PROVIDE tc_grad_and_lapl_ao tc_grad_and_lapl_ao_loop

  norm     = 0.d0
  diff_tot = 0.d0
  do i = 1, ao_num
    do j = 1, ao_num
      do k = 1, ao_num
        do l = 1, ao_num

          diff = dabs(tc_grad_and_lapl_ao_loop(l,k,j,i) - tc_grad_and_lapl_ao(l,k,j,i))
          if(diff .gt. thr_ih) then
            print *, ' difference on ', l, k, j, i
            print *, ' loops : ', tc_grad_and_lapl_ao_loop(l,k,j,i)
            print *, ' lapack: ', tc_grad_and_lapl_ao     (l,k,j,i)
            !stop
          endif

          norm     += dabs(tc_grad_and_lapl_ao_loop(l,k,j,i))
          diff_tot += diff
        enddo
      enddo
    enddo
  enddo

  print *, ' diff tot = ', diff_tot / norm
  print *, '     norm = ', norm
  print *, ' '

  return

end

! ---

subroutine test_tc_grad_square_ao()

  implicit none
  integer          :: i, j, k, l
  double precision :: diff_tot, diff, thr_ih, norm

  thr_ih = 1d-10

  PROVIDE tc_grad_square_ao tc_grad_square_ao_loop

  norm     = 0.d0
  diff_tot = 0.d0
  do i = 1, ao_num
    do j = 1, ao_num
      do k = 1, ao_num
        do l = 1, ao_num

          diff = dabs(tc_grad_square_ao_loop(l,k,j,i) - tc_grad_square_ao(l,k,j,i))
          if(diff .gt. thr_ih) then
            print *, ' difference on ', l, k, j, i
            print *, ' loops : ', tc_grad_square_ao_loop(l,k,j,i)
            print *, ' lapack: ', tc_grad_square_ao     (l,k,j,i)
            !stop
          endif

          norm     += dabs(tc_grad_square_ao_loop(l,k,j,i))
          diff_tot += diff
        enddo
      enddo
    enddo
  enddo

  print *, ' diff tot = ', diff_tot / norm
  print *, '     norm = ', norm
  print *, ' '

  return

end

! ---

subroutine test_two_e_tc_non_hermit_integral()

  implicit none
  integer          :: i, j
  double precision :: diff_tot, diff, thr_ih, norm

  thr_ih = 1d-10

  PROVIDE two_e_tc_non_hermit_integral_beta two_e_tc_non_hermit_integral_alpha
  PROVIDE two_e_tc_non_hermit_integral_seq_beta two_e_tc_non_hermit_integral_seq_alpha

  ! ---

  norm     = 0.d0
  diff_tot = 0.d0
  do i = 1, ao_num
    do j = 1, ao_num

      diff = dabs(two_e_tc_non_hermit_integral_seq_alpha(j,i) - two_e_tc_non_hermit_integral_alpha(j,i))
      if(diff .gt. thr_ih) then
        print *, ' difference on ', j, i
        print *, ' seq         : ', two_e_tc_non_hermit_integral_seq_alpha(j,i)
        print *, ' //          : ', two_e_tc_non_hermit_integral_alpha    (j,i)
        !stop
      endif

      norm     += dabs(two_e_tc_non_hermit_integral_seq_alpha(j,i))
      diff_tot += diff
    enddo
  enddo

  print *, ' diff tot a = ', diff_tot / norm
  print *, '     norm a = ', norm
  print *, ' '

  ! ---

  norm     = 0.d0
  diff_tot = 0.d0
  do i = 1, ao_num
    do j = 1, ao_num

      diff = dabs(two_e_tc_non_hermit_integral_seq_beta(j,i) - two_e_tc_non_hermit_integral_beta(j,i))
      if(diff .gt. thr_ih) then
        print *, ' difference on ', j, i
        print *, ' seq         : ', two_e_tc_non_hermit_integral_seq_beta(j,i)
        print *, ' //          : ', two_e_tc_non_hermit_integral_beta    (j,i)
        !stop
      endif

      norm     += dabs(two_e_tc_non_hermit_integral_seq_beta(j,i))
      diff_tot += diff
    enddo
  enddo

  print *, ' diff tot b = ', diff_tot / norm
  print *, '     norm b = ', norm
  print *, ' '

  ! ---

  return

end

! ---

>>>>>>> 92a4e33f8a21717cab0c0e4f8412ed6903afb04a
