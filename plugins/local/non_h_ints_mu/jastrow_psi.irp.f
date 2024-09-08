BEGIN_PROVIDER [ double precision, c_ij_ab_jastrow, (mo_num, mo_num, elec_alpha_num, elec_beta_num)]
 implicit none
 integer :: iunit, getUnitAndOpen
 c_ij_ab_jastrow = 0.d0
 iunit = getUnitAndOpen(trim(ezfio_work_dir)//'c_ij_ab', 'R')                                                     
 read(iunit) c_ij_ab_jastrow
 close(iunit)
 print*,'c_ij_ab_jastrow = '
 integer :: i,j,a,b
 do i = 1, elec_beta_num ! r2
  do j = 1, elec_alpha_num ! r1
   do a = elec_beta_num+1, mo_num ! r2
    do b = elec_alpha_num+1, mo_num ! r1 
!     print*,b,a,j,i
     print*,c_ij_ab_jastrow(b,a,j,i),b,a,j,i
     if(dabs(c_ij_ab_jastrow(b,a,j,i)).lt.1.d-12)then
      c_ij_ab_jastrow(b,a,j,i) = 0.d0
     endif
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 

double precision function jastrow_psi(r1,r2)
 implicit none
 double precision, intent(in) :: r1(3), r2(3)
 integer :: i,j,a,b
 double precision, allocatable :: mos_array_r1(:), mos_array_r2(:)
 allocate(mos_array_r1(mo_num), mos_array_r2(mo_num))
 call give_all_mos_at_r(r1,mos_array_r1)
 call give_all_mos_at_r(r2,mos_array_r2)
 double precision :: eps,coef, numerator,denominator
 double precision :: phi_i_phi_j
 eps = a_boys
 jastrow_psi= 0.d0
 do i = 1, elec_beta_num ! r1
  do j = 1, elec_alpha_num ! r2
   phi_i_phi_j = mos_array_r1(i) * mos_array_r2(j) + eps
   denominator = 1.d0/phi_i_phi_j
   do a = elec_beta_num+1, mo_num ! r1
    do b = elec_alpha_num+1, mo_num ! r2 
     coef = c_ij_ab_jastrow(b,a,j,i)
     numerator = mos_array_r2(b) * mos_array_r1(a)
     jastrow_psi += coef * numerator*denominator
    enddo
   enddo
  enddo
 enddo
end

subroutine get_grad_r1_jastrow_psi(r1,r2,grad_j_psi_r1,jast)
 implicit none
 double precision, intent(in) :: r1(3), r2(3)
 double precision, intent(out):: grad_j_psi_r1(3),jast
 integer :: i,j,a,b
 double precision, allocatable :: mos_array_r1(:), mos_array_r2(:)
 double precision, allocatable :: mos_grad_array_r1(:,:),mos_grad_array_r2(:,:)
 double precision :: num_j, denom_j, num_j_grad(3), denom_j_grad(3),delta,coef
 double precision :: inv_denom_j
 allocate(mos_array_r1(mo_num), mos_array_r2(mo_num))
 allocate(mos_grad_array_r1(3,mo_num), mos_grad_array_r2(3,mo_num))
 delta = a_boys
 call give_all_mos_and_grad_at_r(r1,mos_array_r1,mos_grad_array_r1)
 call give_all_mos_and_grad_at_r(r2,mos_array_r2,mos_grad_array_r2)
 grad_j_psi_r1 = 0.d0
 jast = 0.d0
 do i = 1, elec_beta_num ! r1
  do j = 1, elec_alpha_num ! r2
   call denom_jpsi(i,j,delta,mos_array_r1,mos_grad_array_r1,mos_array_r2,denom_j, denom_j_grad)
   inv_denom_j = 1.d0/denom_j
   do a = elec_beta_num+1, mo_num ! r1
    do b = elec_alpha_num+1, mo_num ! r2 
     call numerator_psi(a,b,mos_array_r1,mos_grad_array_r1,mos_array_r2,num_j, num_j_grad)
     coef = c_ij_ab_jastrow(b,a,j,i)
     jast += coef * num_j * inv_denom_j
     grad_j_psi_r1 += coef * (num_j_grad * denom_j - num_j * denom_j_grad) * inv_denom_j * inv_denom_j
    enddo
   enddo
  enddo
 enddo
 if(jast.lt.-1.d0.or.dabs(jast).gt.1.d0)then
  print*,'pb ! '
  print*,jast
  print*,dsqrt(r1(1)**2+r1(2)**2+r1(3)**2),dsqrt(r2(1)**2+r2(2)**2+r2(3)**2)
  print*,r1
!  print*,mos_array_r1(1:2)
  print*,r2
!  print*,mos_array_r2(1:2)
  stop
 endif
 if(log_jpsi)then
  grad_j_psi_r1 = grad_j_psi_r1/(1.d0 + jast)
 endif

end


subroutine denom_jpsi(i,j,delta,mos_array_r1,mos_grad_array_r1,mos_array_r2,denom, grad_denom)
 implicit none
 integer, intent(in)           :: i,j
 double precision, intent(in)  :: mos_array_r1(mo_num),mos_grad_array_r1(3,mo_num),mos_array_r2(mo_num),delta
 double precision, intent(out) :: denom, grad_denom(3)
 double precision :: coef,phi_i_phi_j,inv_phi_i_phi_j,inv_phi_i_phi_j_2
 phi_i_phi_j = mos_array_r1(i) * mos_array_r2(j)
 if(phi_i_phi_j /= 0.d0)then
  inv_phi_i_phi_j   = 1.d0/phi_i_phi_j
  inv_phi_i_phi_j_2 = 1.d0/(phi_i_phi_j * phi_i_phi_j)
 else
  inv_phi_i_phi_j   = huge(1.0)
  inv_phi_i_phi_j_2 = huge(1.d0)
 endif
 denom         = phi_i_phi_j + delta * inv_phi_i_phi_j
 grad_denom(:) = (1.d0 - delta*inv_phi_i_phi_j_2) * mos_array_r2(j) * mos_grad_array_r1(:,i)
end

subroutine numerator_psi(a,b,mos_array_r1,mos_grad_array_r1,mos_array_r2,num, grad_num)
 implicit none
 integer, intent(in)           :: a,b
 double precision, intent(in)  :: mos_array_r1(mo_num),mos_grad_array_r1(3,mo_num),mos_array_r2(mo_num)
 double precision, intent(out) :: num, grad_num(3)
 num = mos_array_r1(a) * mos_array_r2(b)
 grad_num(:) = mos_array_r2(b) * mos_grad_array_r1(:,a)
end
