subroutine dm_dft_alpha_beta_at_r(r,dm_a,dm_b)
 implicit none
 BEGIN_DOC
! input: r(1) ==> r(1) = x, r(2) = y, r(3) = z
! output : dm_a = alpha density evaluated at r(3)
! output : dm_b = beta  density evaluated at r(3)
 END_DOC
 double precision, intent(in) :: r(3)
 double precision, intent(out) :: dm_a(N_states),dm_b(N_states)
 integer :: istate
 double precision  :: aos_array(ao_num),aos_array_bis(ao_num),u_dot_v
 call give_all_aos_at_r(r,aos_array)
 do istate = 1, N_states
  aos_array_bis = aos_array
  ! alpha density
  call dgemv('N',ao_num,ao_num,1.d0,one_e_dm_alpha_ao_for_dft(1,1,istate),ao_num,aos_array,1,0.d0,aos_array_bis,1)
  dm_a(istate) = u_dot_v(aos_array,aos_array_bis,ao_num)
  ! beta density
  aos_array_bis = aos_array
  call dgemv('N',ao_num,ao_num,1.d0,one_e_dm_beta_ao_for_dft(1,1,istate),ao_num,aos_array,1,0.d0,aos_array_bis,1)
  dm_b(istate) = u_dot_v(aos_array,aos_array_bis,ao_num)
 enddo
end


subroutine dm_dft_alpha_beta_and_all_aos_at_r(r,dm_a,dm_b,aos_array)
 BEGIN_DOC
! input: r(1) ==> r(1) = x, r(2) = y, r(3) = z
! output : dm_a = alpha density evaluated at r
! output : dm_b = beta  density evaluated at r
! output : aos_array(i) = ao(i) evaluated at r
 END_DOC
 implicit none
 double precision, intent(in) :: r(3)
 double precision, intent(out) :: dm_a(N_states),dm_b(N_states)
 double precision, intent(out) :: aos_array(ao_num)
 integer :: istate
 double precision  :: aos_array_bis(ao_num),u_dot_v
 call give_all_aos_at_r(r,aos_array)
 do istate = 1, N_states
  aos_array_bis = aos_array
  ! alpha density
  call dsymv('U',ao_num,1.d0,one_e_dm_alpha_ao_for_dft(1,1,istate),size(one_e_dm_alpha_ao_for_dft,1),aos_array,1,0.d0,aos_array_bis,1)
  dm_a(istate) = u_dot_v(aos_array,aos_array_bis,ao_num)
  ! beta density
  aos_array_bis = aos_array
  call dsymv('U',ao_num,1.d0,one_e_dm_beta_ao_for_dft(1,1,istate),size(one_e_dm_beta_ao_for_dft,1),aos_array,1,0.d0,aos_array_bis,1)
  dm_b(istate) = u_dot_v(aos_array,aos_array_bis,ao_num)
 enddo
end



 subroutine density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,dm_a,dm_b, grad_dm_a, grad_dm_b, aos_array, grad_aos_array)
 implicit none
 BEGIN_DOC
! input:
!
! * r(1) ==> r(1) = x, r(2) = y, r(3) = z
!
! output:
!
! * dm_a = alpha density evaluated at r
! * dm_b = beta  density evaluated at r
! * aos_array(i) = ao(i) evaluated at r
! * grad_dm_a(1) = X gradient of the alpha density evaluated in r
! * grad_dm_a(1) = X gradient of the beta  density evaluated in r
! * grad_aos_array(1) = X gradient of the aos(i) evaluated at r
!
 END_DOC
 double precision, intent(in)  :: r(3)
 double precision, intent(out) :: dm_a(N_states),dm_b(N_states)
 double precision, intent(out) :: grad_dm_a(3,N_states),grad_dm_b(3,N_states)
 double precision, intent(out) :: grad_aos_array(3,ao_num)
 integer :: i,j,istate
 double precision  :: aos_array(ao_num),aos_array_bis(ao_num),u_dot_v
 double precision  :: aos_grad_array(ao_num,3), aos_grad_array_bis(ao_num,3)

 call give_all_aos_and_grad_at_r(r,aos_array,grad_aos_array)
 do i = 1, ao_num
  do j = 1, 3
   aos_grad_array(i,j) =  grad_aos_array(j,i)
  enddo
 enddo

 do istate = 1, N_states
   ! alpha density
   ! aos_array_bis = \rho_ao * aos_array
   call dsymv('U',ao_num,1.d0,one_e_dm_alpha_ao_for_dft(1,1,istate),size(one_e_dm_alpha_ao_for_dft,1),aos_array,1,0.d0,aos_array_bis,1)
   dm_a(istate) = u_dot_v(aos_array,aos_array_bis,ao_num)

   ! grad_dm(1) = \sum_i aos_grad_array(i,1) * aos_array_bis(i)
   grad_dm_a(1,istate) = u_dot_v(aos_grad_array(1,1),aos_array_bis,ao_num)
   grad_dm_a(2,istate) = u_dot_v(aos_grad_array(1,2),aos_array_bis,ao_num)
   grad_dm_a(3,istate) = u_dot_v(aos_grad_array(1,3),aos_array_bis,ao_num)
   ! aos_grad_array_bis = \rho_ao * aos_grad_array

   ! beta density
   call dsymv('U',ao_num,1.d0,one_e_dm_beta_ao_for_dft(1,1,istate),size(one_e_dm_beta_ao_for_dft,1),aos_array,1,0.d0,aos_array_bis,1)
   dm_b(istate) = u_dot_v(aos_array,aos_array_bis,ao_num)

   ! grad_dm(1) = \sum_i aos_grad_array(i,1) * aos_array_bis(i)
   grad_dm_b(1,istate) = u_dot_v(aos_grad_array(1,1),aos_array_bis,ao_num)
   grad_dm_b(2,istate) = u_dot_v(aos_grad_array(1,2),aos_array_bis,ao_num)
   grad_dm_b(3,istate) = u_dot_v(aos_grad_array(1,3),aos_array_bis,ao_num)
   ! aos_grad_array_bis = \rho_ao * aos_grad_array
 enddo
   grad_dm_a *= 2.d0
   grad_dm_b *= 2.d0
 end

subroutine dm_dft_alpha_beta_no_core_at_r(r,dm_a,dm_b)
 implicit none
 BEGIN_DOC
! input: r(1) ==> r(1) = x, r(2) = y, r(3) = z
! output : dm_a = alpha density evaluated at r(3) without the core orbitals 
! output : dm_b = beta  density evaluated at r(3) without the core orbitals 
 END_DOC
 double precision, intent(in) :: r(3)
 double precision, intent(out) :: dm_a(N_states),dm_b(N_states)
 integer :: istate
 double precision  :: aos_array(ao_num),aos_array_bis(ao_num),u_dot_v
 call give_all_aos_at_r(r,aos_array)
 do istate = 1, N_states
  aos_array_bis = aos_array
  ! alpha density
  call dgemv('N',ao_num,ao_num,1.d0,one_e_dm_alpha_ao_for_dft_no_core(1,1,istate),ao_num,aos_array,1,0.d0,aos_array_bis,1)
  dm_a(istate) = u_dot_v(aos_array,aos_array_bis,ao_num)
  ! beta density
  aos_array_bis = aos_array
  call dgemv('N',ao_num,ao_num,1.d0,one_e_dm_beta_ao_for_dft_no_core(1,1,istate),ao_num,aos_array,1,0.d0,aos_array_bis,1)
  dm_b(istate) = u_dot_v(aos_array,aos_array_bis,ao_num)
 enddo
end

 subroutine dens_grad_a_b_no_core_and_aos_grad_aos_at_r(r,dm_a,dm_b, grad_dm_a, grad_dm_b, aos_array, grad_aos_array)
 implicit none
 BEGIN_DOC
! input:
!
! * r(1) ==> r(1) = x, r(2) = y, r(3) = z
!
! output:
!
! * dm_a = alpha density evaluated at r without the core orbitals 
! * dm_b = beta  density evaluated at r without the core orbitals 
! * aos_array(i) = ao(i) evaluated at r without the core orbitals 
! * grad_dm_a(1) = X gradient of the alpha density evaluated in r without the core orbitals 
! * grad_dm_a(1) = X gradient of the beta  density evaluated in r without the core orbitals 
! * grad_aos_array(1) = X gradient of the aos(i) evaluated at r
!
 END_DOC
 double precision, intent(in)  :: r(3)
 double precision, intent(out) :: dm_a(N_states),dm_b(N_states)
 double precision, intent(out) :: grad_dm_a(3,N_states),grad_dm_b(3,N_states)
 double precision, intent(out) :: grad_aos_array(3,ao_num)
 integer :: i,j,istate
 double precision  :: aos_array(ao_num),aos_array_bis(ao_num),u_dot_v
 double precision  :: aos_grad_array(ao_num,3), aos_grad_array_bis(ao_num,3)

 call give_all_aos_and_grad_at_r(r,aos_array,grad_aos_array)
 do i = 1, ao_num
  do j = 1, 3
   aos_grad_array(i,j) =  grad_aos_array(j,i)
  enddo
 enddo

 do istate = 1, N_states
   ! alpha density
   ! aos_array_bis = \rho_ao * aos_array
   call dsymv('U',ao_num,1.d0,one_e_dm_alpha_ao_for_dft_no_core(1,1,istate),size(one_e_dm_alpha_ao_for_dft_no_core,1),aos_array,1,0.d0,aos_array_bis,1)
   dm_a(istate) = u_dot_v(aos_array,aos_array_bis,ao_num)

   ! grad_dm(1) = \sum_i aos_grad_array(i,1) * aos_array_bis(i)
   grad_dm_a(1,istate) = u_dot_v(aos_grad_array(1,1),aos_array_bis,ao_num)
   grad_dm_a(2,istate) = u_dot_v(aos_grad_array(1,2),aos_array_bis,ao_num)
   grad_dm_a(3,istate) = u_dot_v(aos_grad_array(1,3),aos_array_bis,ao_num)
   ! aos_grad_array_bis = \rho_ao * aos_grad_array

   ! beta density
   call dsymv('U',ao_num,1.d0,one_e_dm_beta_ao_for_dft_no_core(1,1,istate),size(one_e_dm_beta_ao_for_dft_no_core,1),aos_array,1,0.d0,aos_array_bis,1)
   dm_b(istate) = u_dot_v(aos_array,aos_array_bis,ao_num)

   ! grad_dm(1) = \sum_i aos_grad_array(i,1) * aos_array_bis(i)
   grad_dm_b(1,istate) = u_dot_v(aos_grad_array(1,1),aos_array_bis,ao_num)
   grad_dm_b(2,istate) = u_dot_v(aos_grad_array(1,2),aos_array_bis,ao_num)
   grad_dm_b(3,istate) = u_dot_v(aos_grad_array(1,3),aos_array_bis,ao_num)
   ! aos_grad_array_bis = \rho_ao * aos_grad_array
 enddo
   grad_dm_a *= 2.d0
   grad_dm_b *= 2.d0
 end



 BEGIN_PROVIDER [double precision, one_e_dm_alpha_in_r, (n_points_integration_angular,n_points_radial_grid,nucl_num,N_states) ]
&BEGIN_PROVIDER [double precision, one_e_dm_beta_in_r, (n_points_integration_angular,n_points_radial_grid,nucl_num,N_states) ]
 implicit none
 integer :: i,j,k,l,m,istate
 double precision :: contrib
 double precision :: r(3)
 double precision :: aos_array(ao_num),mos_array(mo_num)
  do j = 1, nucl_num
   do k = 1, n_points_radial_grid -1
    do l = 1, n_points_integration_angular
     do istate = 1, N_States
      one_e_dm_alpha_in_r(l,k,j,istate) = 0.d0
      one_e_dm_beta_in_r(l,k,j,istate) = 0.d0
     enddo
     r(1) = grid_points_per_atom(1,l,k,j)
     r(2) = grid_points_per_atom(2,l,k,j)
     r(3) = grid_points_per_atom(3,l,k,j)

     double precision :: dm_a(N_states),dm_b(N_states)
     call dm_dft_alpha_beta_at_r(r,dm_a,dm_b)
     do istate=1,N_states
      one_e_dm_alpha_in_r(l,k,j,istate) = dm_a(istate)
      one_e_dm_beta_in_r(l,k,j,istate) = dm_b(istate)
     enddo

    enddo
   enddo
  enddo

END_PROVIDER


 BEGIN_PROVIDER [double precision, one_e_dm_alpha_at_r, (n_points_final_grid,N_states) ]
&BEGIN_PROVIDER [double precision, one_e_dm_beta_at_r, (n_points_final_grid,N_states) ]
&BEGIN_PROVIDER [double precision, elec_beta_num_grid_becke , (N_states) ]
&BEGIN_PROVIDER [double precision, elec_alpha_num_grid_becke , (N_states) ]
 implicit none
 BEGIN_DOC
! one_e_dm_alpha_at_r(i,istate) = n_alpha(r_i,istate)
! one_e_dm_beta_at_r(i,istate) =  n_beta(r_i,istate)
! where r_i is the ith point of the grid and istate is the state number
 END_DOC
 integer :: i,istate
 double precision :: r(3)
 double precision, allocatable :: dm_a(:),dm_b(:)
 allocate(dm_a(N_states),dm_b(N_states))
 do istate = 1, N_states
  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)
   r(2) = final_grid_points(2,i)
   r(3) = final_grid_points(3,i)
   call dm_dft_alpha_beta_at_r(r,dm_a,dm_b)
   one_e_dm_alpha_at_r(i,istate) = dm_a(istate)
   one_e_dm_beta_at_r(i,istate) = dm_b(istate)
  enddo
 enddo

END_PROVIDER


 BEGIN_PROVIDER [double precision, one_e_dm_and_grad_alpha_in_r, (4,n_points_final_grid,N_states) ]
&BEGIN_PROVIDER [double precision, one_e_dm_and_grad_beta_in_r,  (4,n_points_final_grid,N_states) ]
&BEGIN_PROVIDER [double precision, one_e_grad_2_dm_alpha_at_r, (n_points_final_grid,N_states) ]
&BEGIN_PROVIDER [double precision, one_e_grad_2_dm_beta_at_r, (n_points_final_grid,N_states) ]
 BEGIN_DOC
! one_e_dm_and_grad_alpha_in_r(1,i,i_state) = d\dx n_alpha(r_i,istate)
! one_e_dm_and_grad_alpha_in_r(2,i,i_state) = d\dy n_alpha(r_i,istate)
! one_e_dm_and_grad_alpha_in_r(3,i,i_state) = d\dz n_alpha(r_i,istate)
! one_e_dm_and_grad_alpha_in_r(4,i,i_state) = n_alpha(r_i,istate)
! one_e_grad_2_dm_alpha_at_r(i,istate)      = d\dx n_alpha(r_i,istate)^2 + d\dy n_alpha(r_i,istate)^2 + d\dz n_alpha(r_i,istate)^2
! where r_i is the ith point of the grid and istate is the state number
 END_DOC
 implicit none
 integer :: i,j,k,l,m,istate
 double precision :: contrib
 double precision :: r(3)
 double precision, allocatable :: aos_array(:),grad_aos_array(:,:)
 double precision, allocatable :: dm_a(:),dm_b(:), dm_a_grad(:,:), dm_b_grad(:,:)
 allocate(dm_a(N_states),dm_b(N_states), dm_a_grad(3,N_states), dm_b_grad(3,N_states))
 allocate(aos_array(ao_num),grad_aos_array(3,ao_num))
 do istate = 1, N_states
  do i = 1, n_points_final_grid
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
 !!!! Works also with the ao basis
   call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,dm_a,dm_b,  dm_a_grad, dm_b_grad, aos_array, grad_aos_array)
   one_e_dm_and_grad_alpha_in_r(1,i,istate)  =  dm_a_grad(1,istate)
   one_e_dm_and_grad_alpha_in_r(2,i,istate)  =  dm_a_grad(2,istate)
   one_e_dm_and_grad_alpha_in_r(3,i,istate)  =  dm_a_grad(3,istate)
   one_e_dm_and_grad_alpha_in_r(4,i,istate)  =  dm_a(istate)
   one_e_grad_2_dm_alpha_at_r(i,istate) = dm_a_grad(1,istate) * dm_a_grad(1,istate) + dm_a_grad(2,istate) * dm_a_grad(2,istate) + dm_a_grad(3,istate) * dm_a_grad(3,istate)

   one_e_dm_and_grad_beta_in_r(1,i,istate)  =  dm_b_grad(1,istate)
   one_e_dm_and_grad_beta_in_r(2,i,istate)  =  dm_b_grad(2,istate)
   one_e_dm_and_grad_beta_in_r(3,i,istate)  =  dm_b_grad(3,istate)
   one_e_dm_and_grad_beta_in_r(4,i,istate)  =  dm_b(istate)
   one_e_grad_2_dm_beta_at_r(i,istate) = dm_b_grad(1,istate) * dm_b_grad(1,istate) + dm_b_grad(2,istate) * dm_b_grad(2,istate) + dm_b_grad(3,istate) * dm_b_grad(3,istate)
  enddo
 enddo

END_PROVIDER


 BEGIN_PROVIDER [double precision, one_e_dm_no_core_and_grad_alpha_in_r, (4,n_points_final_grid,N_states) ]
&BEGIN_PROVIDER [double precision, one_e_dm_no_core_and_grad_beta_in_r,  (4,n_points_final_grid,N_states) ]
 BEGIN_DOC
! one_e_dm_no_core_and_grad_alpha_in_r(1,i,i_state) = d\dx n_alpha(r_i,istate) without core orbitals 
! one_e_dm_no_core_and_grad_alpha_in_r(2,i,i_state) = d\dy n_alpha(r_i,istate) without core orbitals 
! one_e_dm_no_core_and_grad_alpha_in_r(3,i,i_state) = d\dz n_alpha(r_i,istate) without core orbitals 
! one_e_dm_no_core_and_grad_alpha_in_r(4,i,i_state) = n_alpha(r_i,istate) without core orbitals 
! where r_i is the ith point of the grid and istate is the state number
 END_DOC
 implicit none
 integer :: i,j,k,l,m,istate
 double precision :: contrib
 double precision :: r(3)
 double precision, allocatable :: aos_array(:),grad_aos_array(:,:)
 double precision, allocatable :: dm_a(:),dm_b(:), dm_a_grad(:,:), dm_b_grad(:,:)
 allocate(dm_a(N_states),dm_b(N_states), dm_a_grad(3,N_states), dm_b_grad(3,N_states))
 allocate(aos_array(ao_num),grad_aos_array(3,ao_num))
 do istate = 1, N_states
  do i = 1, n_points_final_grid
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
 !!!! Works also with the ao basis
   call dens_grad_a_b_no_core_and_aos_grad_aos_at_r(r,dm_a,dm_b,  dm_a_grad, dm_b_grad, aos_array, grad_aos_array)
   one_e_dm_no_core_and_grad_alpha_in_r(1,i,istate)  =  dm_a_grad(1,istate)
   one_e_dm_no_core_and_grad_alpha_in_r(2,i,istate)  =  dm_a_grad(2,istate)
   one_e_dm_no_core_and_grad_alpha_in_r(3,i,istate)  =  dm_a_grad(3,istate)
   one_e_dm_no_core_and_grad_alpha_in_r(4,i,istate)  =  dm_a(istate)

   one_e_dm_no_core_and_grad_beta_in_r(1,i,istate)  =  dm_b_grad(1,istate)
   one_e_dm_no_core_and_grad_beta_in_r(2,i,istate)  =  dm_b_grad(2,istate)
   one_e_dm_no_core_and_grad_beta_in_r(3,i,istate)  =  dm_b_grad(3,istate)
   one_e_dm_no_core_and_grad_beta_in_r(4,i,istate)  =  dm_b(istate)
  enddo
 enddo

END_PROVIDER


