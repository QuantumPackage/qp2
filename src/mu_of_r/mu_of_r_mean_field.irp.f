BEGIN_PROVIDER [ double precision, two_e_int_mf, (elec_beta_num,elec_alpha_num,elec_beta_num,elec_alpha_num)]
 implicit none
 integer :: i,j,k,l 
 double precision :: get_two_e_integral
 do i = 1, elec_alpha_num
  do j = 1, elec_beta_num
   do k = 1, elec_alpha_num
    do l = 1, elec_beta_num
     two_e_int_mf(l,k,j,i) = get_two_e_integral(l,k,j,i,mo_integrals_map) 
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 

subroutine get_f_mf_ab(r,f_mf_ab,two_bod_dens, dm_a, dm_b)
 implicit none
 double precision, intent(in) :: r(3)
 double precision, intent(out):: f_mf_ab,two_bod_dens, dm_a, dm_b
 double precision, allocatable :: mos_array_r(:),mos_array_a(:), mos_array_b(:)
 integer :: i,j,k,l
 allocate(mos_array_r(mo_num), mos_array_a(elec_alpha_num), mos_array_b(elec_alpha_num))
 call give_all_mos_at_r(r,mos_array_r) 
 do i = 1, elec_alpha_num
  mos_array_a(i) = mos_array_r(i)
 enddo
 do i = 1, elec_beta_num
  mos_array_b(i) = mos_array_r(i)
 enddo

 dm_a = 0.d0
 do i = 1, elec_alpha_num
  dm_a += mos_array_a(i) * mos_array_a(i) 
 enddo

 dm_b = 0.d0
 do i = 1, elec_beta_num
  dm_b += mos_array_b(i) * mos_array_b(i) 
 enddo
 two_bod_dens = dm_a * dm_b
 
 f_mf_ab = 0.d0
 do i = 1, elec_alpha_num
  do j = 1, elec_beta_num
   do k = 1, elec_alpha_num
    do l = 1, elec_beta_num
     f_mf_ab += two_e_int_mf(l,k,j,i) * mos_array_a(i) * mos_array_a(k) * mos_array_b(j) * mos_array_b(l)
    enddo
   enddo
  enddo
 enddo
 ! multiply by two to adapt to the N(N-1) normalization condition of the active two-rdm
 f_mf_ab *= 2.d0 
 two_bod_dens *= 2.d0

end

subroutine get_grad_f_mf_ab(r,grad_f_mf_ab, grad_two_bod_dens,f_mf_ab,two_bod_dens, dm_a, dm_b,grad_dm_a, grad_dm_b)
 implicit none
 double precision, intent(in)  :: r(3)
 double precision, intent(out) :: f_mf_ab, two_bod_dens
 double precision, intent(out) :: grad_two_bod_dens(3), grad_f_mf_ab(3)
 double precision, intent(out) :: dm_a, dm_b, grad_dm_a(3), grad_dm_b(3)

 double precision, allocatable :: mos_array_r(:), mos_grad_array_r(:,:)
 double precision, allocatable :: mos_array_a(:), mos_array_b(:)
 double precision, allocatable :: mos_grad_array_a(:,:), mos_grad_array_b(:,:)
 double precision :: mo_i, mo_j, mo_k, mo_l
 double precision :: grad_mo_i(3), grad_mo_j(3), grad_mo_k(3), grad_mo_l(3)
 
 integer :: i,j,k,l
 allocate(mos_array_r(mo_num),mos_grad_array_r(3,mo_num))
 allocate(mos_array_a(elec_alpha_num), mos_array_b(elec_beta_num))
 allocate(mos_grad_array_a(3,elec_alpha_num), mos_grad_array_b(3,elec_beta_num))
 call give_all_mos_and_grad_at_r(r,mos_array_r,mos_grad_array_r)
 do i = 1, elec_alpha_num
  mos_array_a(i) = mos_array_r(i)
  mos_grad_array_a(1:3,i) = mos_grad_array_r(1:3,i)
 enddo
 do i = 1, elec_beta_num
  mos_array_b(i) = mos_array_r(i)
  mos_grad_array_b(1:3,i) = mos_grad_array_r(1:3,i)
 enddo

 ! ALPHA DENSITY AND GRADIENT 
 dm_a = 0.d0
 grad_dm_a = 0.d0
 do i = 1, elec_alpha_num
  dm_a += mos_array_a(i) * mos_array_a(i) 
  grad_dm_a(1:3) += 2.d0 * mos_array_a(i) * mos_grad_array_a(1:3,i)  
 enddo

 ! BETA DENSITY AND GRADIENT 
 dm_b = 0.d0
 grad_dm_b = 0.d0
 do i = 1, elec_beta_num
  dm_b += mos_array_b(i) * mos_array_b(i) 
  grad_dm_b(1:3) += 2.d0 * mos_array_b(i) * mos_grad_array_b(1:3,i)  
 enddo
 ! TWO-BODY DENSITY AND GRADIENT 
 two_bod_dens = dm_a * dm_b
 grad_two_bod_dens(1:3) = dm_a * grad_dm_b(1:3) + dm_b * grad_dm_a(1:3)

 ! F_MF and GRADIENT 
 grad_f_mf_ab = 0.d0
 f_mf_ab  = 0.d0
 do i = 1, elec_alpha_num
  mo_i = mos_array_a(i)
  grad_mo_i(1:3) = mos_grad_array_a(1:3,i)
  do j = 1, elec_beta_num
   mo_j = mos_array_b(j)
   grad_mo_j(1:3) = mos_grad_array_b(1:3,j)
   do k = 1, elec_alpha_num
    mo_k = mos_array_a(k)
    grad_mo_k(1:3) = mos_grad_array_a(1:3,k)
    do l = 1, elec_beta_num
     mo_l = mos_array_b(l)
     grad_mo_l(1:3) = mos_grad_array_b(1:3,l)
     f_mf_ab += two_e_int_mf(l,k,j,i) * mo_i * mo_j * mo_k * mo_l
     grad_f_mf_ab(1:3) += two_e_int_mf(l,k,j,i) * & 
     (mo_i * mo_j * mo_k * grad_mo_l(1:3) + mo_i * mo_j * grad_mo_k(1:3) * mo_l & 
     +mo_i * grad_mo_j(1:3) * mo_k * mo_l + grad_mo_i(1:3) * mo_j * mo_k * mo_l)
    enddo
   enddo
  enddo
 enddo

 f_mf_ab *= 2.d0 
 two_bod_dens *= 2.d0
 grad_f_mf_ab *= 2.D0
 grad_two_bod_dens *= 2.d0
end
