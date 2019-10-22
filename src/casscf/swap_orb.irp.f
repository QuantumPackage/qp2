 BEGIN_PROVIDER [double precision, super_ci_dm, (mo_num,mo_num)]
 implicit none 
 super_ci_dm = 0.d0
 integer :: i,j
 integer :: iorb,jorb
 integer :: a,aorb
 double precision :: c0,ci
 c0 = SXeigenvec(1,1)
 ! equation B3.a of the annex B of Roos et. al. Chemical Physics 48 (1980) 157-173
 do i = 1, n_core_inact_orb
  iorb = list_core_inact(i)
  super_ci_dm(iorb,iorb) = 2.d0
  do j = 1, n_core_inact_orb
   jorb = list_core_inact(j)
   ! loop over the core/inact 
   do a = 1, n_virt_orb
    aorb = list_virt(a)
    super_ci_dm(jorb,iorb) += lowest_super_ci_coef_mo(aorb,iorb) * lowest_super_ci_coef_mo(aorb,jorb)
   enddo
  enddo
 enddo

 END_PROVIDER 
 
 BEGIN_PROVIDER [double precision, lowest_super_ci_coef_mo, (mo_num,mo_num)]
 implicit none
 integer :: i,j,iorb,jorb
 integer :: a, aorb,t, torb
 double precision :: sqrt2

 sqrt2 = 1.d0/dsqrt(2.d0)
 do i = 1, nMonoEx
  iorb = excit(1,i)
  jorb = excit(2,i)
  super_ci_coef_mo(iorb,jorb) = SXeigenvec(i+1,1)
  super_ci_coef_mo(jorb,iorb) = SXeigenvec(i+1,1)
 enddo

 ! a_{ia} of the equation B.2 of Roos et. al. Chemical Physics 48 (1980) 157-173
 do i = 1, n_core_inact_orb
  iorb = list_core_inact(i)
  do a = 1, n_virt_orb
   aorb = list_virt(a)
   super_ci_coef_mo(aorb,iorb) *= sqrt2
   super_ci_coef_mo(iorb,aorb) *= sqrt2
  enddo
 enddo
 
 ! a_{it} of the equation B.2 of Roos et. al. Chemical Physics 48 (1980) 157-173
 do i = 1, n_core_inact_orb
  iorb = list_core_inact(i)
  do t = 1, n_act_orb
   torb = list_act(t)
   super_ci_coef_mo(torb,iorb) *= (2.d0 - occ_act(t))**(-0.5d0)
   super_ci_coef_mo(iorb,torb) *= (2.d0 - occ_act(t))**(-0.5d0)
  enddo
 enddo

 ! a_{ta} of the equation B.2 of Roos et. al. Chemical Physics 48 (1980) 157-173
 do a = 1, n_virt_orb
  aorb = list_virt(a)
  do t = 1, n_act_orb
   torb = list_act(t)
   super_ci_coef_mo(torb,aorb) *= occ_act(t)**(-0.5d0)
   super_ci_coef_mo(aorb,torb) *= occ_act(t)**(-0.5d0)
  enddo
 enddo

 END_PROVIDER 

