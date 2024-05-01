program plot_j
 implicit none
 double precision :: r1(3),rI(3),r2(3)
 double precision :: r12,dx,xmax, j_1e,j_2e,j_een,j_tot
 double precision :: j_mu_F_x_j
 integer :: i,nx,m,i_charge,sm_j

 character*(128) :: output
 integer :: i_unit_output_He_sm_7,i_unit_output_Ne_sm_7
 integer :: i_unit_output_He_sm_17,i_unit_output_Ne_sm_17
 integer :: getUnitAndOpen
 output='J_SM_7_He'
 i_unit_output_He_sm_7 = getUnitAndOpen(output,'w')
 output='J_SM_7_Ne'
 i_unit_output_Ne_sm_7 = getUnitAndOpen(output,'w')

 output='J_SM_17_He'
 i_unit_output_He_sm_17 = getUnitAndOpen(output,'w')
 output='J_SM_17_Ne'
 i_unit_output_Ne_sm_17 = getUnitAndOpen(output,'w')

 rI = 0.d0
 r1 = 0.d0
 r2 = 0.d0
 r1(1) = 1.5d0
 xmax = 20.d0
 r2(1) = -xmax*0.5d0
 nx = 1000
 dx = xmax/dble(nx)
 do i = 1, nx
  r12 = 0.d0
  do m = 1, 3
   r12 += (r1(m) - r2(m))*(r1(m) - r2(m))
  enddo
  r12 = dsqrt(r12)
  double precision :: jmu,env_nucl,jmu_env,jmu_scaled, jmu_scaled_env
  double precision :: b_I,d_I,r_inucl,r_jnucl,r_ij
  b_I = 1.D0
  d_I = 1.D0
  call get_rescaled_variables_j_sm_90(r1,r2,rI,b_I,d_I,r_inucl,r_jnucl,r_ij)
  jmu=j_mu_F_x_j(r12)
  jmu_scaled=j_mu_F_x_j(r_ij)
  jmu_env = jmu * env_nucl(r1) * env_nucl(r2)
!  jmu_scaled_env= jmu_scaled * (1.d0 - env_coef(1) * dexp(-env_expo(1)*r_inucl**2)) * (1.d0 - env_coef(1) * dexp(-env_expo(1)*r_jnucl**2))
  jmu_scaled_env= jmu_scaled * env_nucl(r1) * env_nucl(r2)
  ! He 
  i_charge = 2
  ! SM 7 Jastrow 
  sm_j = 7
  call get_full_sm_90_jastrow(r1,r2,rI,sm_j,i_charge, j_1e,j_2e,j_een,j_tot)
  write(i_unit_output_He_sm_7,'(100(F16.10,X))')r2(1),r12,j_mu_F_x_j(r12), j_1e,j_2e,j_een,j_tot,jmu_env,jmu_scaled,jmu_scaled_env
  ! SM 17 Jastrow 
  sm_j = 17
  call get_full_sm_90_jastrow(r1,r2,rI,sm_j,i_charge, j_1e,j_2e,j_een,j_tot)
  write(i_unit_output_He_sm_17,'(100(F16.10,X))')r2(1),r12,j_mu_F_x_j(r12), j_1e,j_2e,j_een,j_tot,jmu_env,jmu_scaled,jmu_scaled_env
  ! Ne 
  i_charge = 10
  ! SM 7 Jastrow 
  sm_j = 7
  call get_full_sm_90_jastrow(r1,r2,rI,sm_j,i_charge, j_1e,j_2e,j_een,j_tot)
  write(i_unit_output_Ne_sm_7,'(100(F16.10,X))')r2(1),r12,j_mu_F_x_j(r12), j_1e,j_2e,j_een,j_tot,jmu_env,jmu_scaled,jmu_scaled_env
  ! SM 17 Jastrow 
  sm_j = 17
  call get_full_sm_90_jastrow(r1,r2,rI,sm_j,i_charge, j_1e,j_2e,j_een,j_tot)
  write(i_unit_output_Ne_sm_17,'(100(F16.10,X))')r2(1),r12,j_mu_F_x_j(r12), j_1e,j_2e,j_een,j_tot,jmu_env,jmu_scaled,jmu_scaled_env
  r2(1) += dx
 enddo
 
end
