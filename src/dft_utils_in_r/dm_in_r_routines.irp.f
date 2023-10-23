
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

 subroutine density_and_grad_alpha_beta(r,dm_a,dm_b, grad_dm_a, grad_dm_b)
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
! * grad_dm_a(1) = X gradient of the alpha density evaluated in r
! * grad_dm_a(1) = X gradient of the beta  density evaluated in r
!
 END_DOC
 double precision, intent(in)  :: r(3)
 double precision, intent(out) :: dm_a(N_states),dm_b(N_states)
 double precision, intent(out) :: grad_dm_a(3,N_states),grad_dm_b(3,N_states)
 double precision              :: grad_aos_array(3,ao_num)
 integer :: i,j,istate
 double precision  :: aos_array(ao_num),aos_array_bis(ao_num),u_dot_v
 double precision  :: aos_grad_array(ao_num,3), aos_grad_array_bis(ao_num,3)

 call give_all_aos_and_grad_at_r(r,aos_array,grad_aos_array)
 do i = 1, ao_num
  do j = 1, 3
   aos_grad_array(i,j) =  grad_aos_array(j,i)
  enddo
 enddo

 ! TODO : build the vector of chi_i(r) chi_j(r) and conscequently grad_i(r) grad_j(r) 
 !      : the same for gamma_ij and big dot product 
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



 subroutine density_and_grad_lapl_alpha_beta_and_all_aos_and_grad_aos_at_r(r,dm_a,dm_b, grad_dm_a, grad_dm_b, lapl_dm_a, lapl_dm_b, aos_array, grad_aos_array, lapl_aos_array)
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
 double precision, intent(out) :: lapl_dm_a(3,N_states),lapl_dm_b(3,N_states)
 double precision, intent(out) :: grad_aos_array(3,ao_num)
 double precision, intent(out) :: lapl_aos_array(3,ao_num)
 integer :: i,j,istate
 double precision  :: aos_array(ao_num),aos_array_bis(ao_num),u_dot_v
 double precision  :: aos_grad_array(ao_num,3), aos_grad_array_bis(ao_num,3)
 double precision  :: aos_lapl_array(ao_num,3)

 call give_all_aos_and_grad_and_lapl_at_r(r,aos_array,grad_aos_array,lapl_aos_array)
 do i = 1, ao_num
  do j = 1, 3
   aos_grad_array(i,j) =  grad_aos_array(j,i)
   aos_lapl_array(i,j) =  lapl_aos_array(j,i)
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

   ! lapl_dm(1) = \sum_i aos_lapl_array(i,1) * aos_array_bis(i)
   lapl_dm_a(1,istate) = 2.d0 * u_dot_v(aos_lapl_array(1,1),aos_array_bis,ao_num)
   lapl_dm_a(2,istate) = 2.d0 * u_dot_v(aos_lapl_array(1,2),aos_array_bis,ao_num)
   lapl_dm_a(3,istate) = 2.d0 * u_dot_v(aos_lapl_array(1,3),aos_array_bis,ao_num)

   ! aos_grad_array_bis(1) = \rho_ao * aos_grad_array(1)
   call dsymv('U',ao_num,1.d0,one_e_dm_alpha_ao_for_dft(1,1,istate),size(one_e_dm_alpha_ao_for_dft,1),aos_grad_array(1,1),1,0.d0,aos_grad_array_bis(1,1),1)
   call dsymv('U',ao_num,1.d0,one_e_dm_alpha_ao_for_dft(1,1,istate),size(one_e_dm_alpha_ao_for_dft,1),aos_grad_array(1,2),1,0.d0,aos_grad_array_bis(1,2),1)
   call dsymv('U',ao_num,1.d0,one_e_dm_alpha_ao_for_dft(1,1,istate),size(one_e_dm_alpha_ao_for_dft,1),aos_grad_array(1,3),1,0.d0,aos_grad_array_bis(1,3),1)
   ! lapl_dm(1) += \sum_i aos_grad_array(i,1) * aos_grad_array_bis(i)
   lapl_dm_a(1,istate) += 2.d0 * u_dot_v(aos_grad_array(1,1),aos_grad_array_bis,ao_num)
   lapl_dm_a(2,istate) += 2.d0 * u_dot_v(aos_grad_array(1,2),aos_grad_array_bis,ao_num)
   lapl_dm_a(3,istate) += 2.d0 * u_dot_v(aos_grad_array(1,3),aos_grad_array_bis,ao_num)
 

   ! beta density
   call dsymv('U',ao_num,1.d0,one_e_dm_beta_ao_for_dft(1,1,istate),size(one_e_dm_beta_ao_for_dft,1),aos_array,1,0.d0,aos_array_bis,1)
   dm_b(istate) = u_dot_v(aos_array,aos_array_bis,ao_num)

   ! grad_dm(1) = \sum_i aos_grad_array(i,1) * aos_array_bis(i)
   grad_dm_b(1,istate) = u_dot_v(aos_grad_array(1,1),aos_array_bis,ao_num)
   grad_dm_b(2,istate) = u_dot_v(aos_grad_array(1,2),aos_array_bis,ao_num)
   grad_dm_b(3,istate) = u_dot_v(aos_grad_array(1,3),aos_array_bis,ao_num)

   ! lapl_dm(1) = \sum_i aos_lapl_array(i,1) * aos_array_bis(i)
   lapl_dm_b(1,istate) = 2.d0 * u_dot_v(aos_lapl_array(1,1),aos_array_bis,ao_num)
   lapl_dm_b(2,istate) = 2.d0 * u_dot_v(aos_lapl_array(1,2),aos_array_bis,ao_num)
   lapl_dm_b(3,istate) = 2.d0 * u_dot_v(aos_lapl_array(1,3),aos_array_bis,ao_num)

   ! aos_grad_array_bis(1) = \rho_ao * aos_grad_array(1)
   call dsymv('U',ao_num,1.d0,one_e_dm_alpha_ao_for_dft(1,1,istate),size(one_e_dm_alpha_ao_for_dft,1),aos_grad_array(1,1),1,0.d0,aos_grad_array_bis(1,1),1)
   call dsymv('U',ao_num,1.d0,one_e_dm_alpha_ao_for_dft(1,1,istate),size(one_e_dm_alpha_ao_for_dft,1),aos_grad_array(1,2),1,0.d0,aos_grad_array_bis(1,2),1)
   call dsymv('U',ao_num,1.d0,one_e_dm_alpha_ao_for_dft(1,1,istate),size(one_e_dm_alpha_ao_for_dft,1),aos_grad_array(1,3),1,0.d0,aos_grad_array_bis(1,3),1)
   ! lapl_dm(1) += \sum_i aos_grad_array(i,1) * aos_grad_array_bis(i)
   lapl_dm_b(1,istate) += 2.d0 * u_dot_v(aos_grad_array(1,1),aos_grad_array_bis,ao_num)
   lapl_dm_b(2,istate) += 2.d0 * u_dot_v(aos_grad_array(1,2),aos_grad_array_bis,ao_num)
   lapl_dm_b(3,istate) += 2.d0 * u_dot_v(aos_grad_array(1,3),aos_grad_array_bis,ao_num)
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


