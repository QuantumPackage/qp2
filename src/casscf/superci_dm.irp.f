 BEGIN_PROVIDER [double precision, super_ci_dm, (mo_num,mo_num)]
 implicit none 
 BEGIN_DOC
! density matrix of the super CI matrix, in the basis of NATURAL ORBITALS OF THE CASCI WF 
! 
! This is obtained from annex B of Roos et. al. Chemical Physics 48 (1980) 157-173
!
! WARNING ::: in the equation B3.d there is a TYPO with a forgotten MINUS SIGN (see variable mat_tmp_dm_super_ci ) 
 END_DOC
 super_ci_dm = 0.d0
 integer :: i,j,iorb,jorb
 integer :: a,aorb,b,borb
 integer :: t,torb,v,vorb,u,uorb,x,xorb
 double precision :: c0,ci
 c0 = SXeigenvec(1,1)
 ! equation B3.a of the annex B of Roos et. al. Chemical Physics 48 (1980) 157-173
 ! loop over the core/inact 
 do i = 1, n_core_inact_orb
  iorb = list_core_inact(i)
  super_ci_dm(iorb,iorb) = 2.d0 ! first term of B3.a
  ! loop over the core/inact 
  do j = 1, n_core_inact_orb
   jorb = list_core_inact(j)
   ! loop over the virtual 
   do a = 1, n_virt_orb
    aorb = list_virt(a)
    super_ci_dm(jorb,iorb) += -2.d0 * lowest_super_ci_coef_mo(aorb,iorb) * lowest_super_ci_coef_mo(aorb,jorb) ! second term in B3.a
   enddo
   do t = 1, n_act_orb
    torb = list_act(t)
    ! thrid term of the B3.a 
    super_ci_dm(jorb,iorb) += - lowest_super_ci_coef_mo(iorb,torb) * lowest_super_ci_coef_mo(jorb,torb) * (2.d0 - occ_act(t))
   enddo
  enddo
 enddo

 ! equation B3.b of the annex B of Roos et. al. Chemical Physics 48 (1980) 157-173
 do i = 1, n_core_inact_orb
  iorb = list_core_inact(i)
  do t = 1, n_act_orb
   torb = list_act(t)
   super_ci_dm(iorb,torb) = c0 * lowest_super_ci_coef_mo(torb,iorb) * (2.d0 - occ_act(t))
   super_ci_dm(torb,iorb) = c0 * lowest_super_ci_coef_mo(torb,iorb) * (2.d0 - occ_act(t))
   do a = 1, n_virt_orb
    aorb = list_virt(a)
    super_ci_dm(iorb,torb) += - lowest_super_ci_coef_mo(aorb,iorb) * lowest_super_ci_coef_mo(aorb,torb) * occ_act(t)
    super_ci_dm(torb,iorb) += - lowest_super_ci_coef_mo(aorb,iorb) * lowest_super_ci_coef_mo(aorb,torb) * occ_act(t)
   enddo
  enddo
 enddo

 ! equation B3.c of the annex B of Roos et. al. Chemical Physics 48 (1980) 157-173
 do i = 1, n_core_inact_orb
  iorb = list_core_inact(i)
  do a = 1, n_virt_orb
   aorb = list_virt(a)
   super_ci_dm(aorb,iorb) = 2.d0 * c0 * lowest_super_ci_coef_mo(aorb,iorb)
   super_ci_dm(iorb,aorb) = 2.d0 * c0 * lowest_super_ci_coef_mo(aorb,iorb)
  enddo
 enddo

 ! equation B3.d of the annex B of Roos et. al. Chemical Physics 48 (1980) 157-173
 do t = 1, n_act_orb
  torb = list_act(t)
  super_ci_dm(torb,torb) = occ_act(t) ! first term of equation B3.d 
  do x = 1, n_act_orb
   xorb = list_act(x) 
   super_ci_dm(torb,torb) += - occ_act(x) * occ_act(t)* mat_tmp_dm_super_ci(x,x) ! second term involving the ONE-rdm 
  enddo
  do u = 1, n_act_orb
   uorb = list_act(u)

   ! second term of equation B3.d 
   do x = 1, n_act_orb
    xorb = list_act(x) 
    do v = 1, n_act_orb
     vorb = list_act(v) 
     super_ci_dm(torb,uorb) +=  2.d0 * P0tuvx_no(v,x,t,u) * mat_tmp_dm_super_ci(v,x) ! second term involving the TWO-rdm 
    enddo
   enddo
   
   ! third term of equation B3.d 
   do i = 1, n_core_inact_orb
    iorb = list_core_inact(i)
    super_ci_dm(torb,uorb) += lowest_super_ci_coef_mo(iorb,torb) * lowest_super_ci_coef_mo(iorb,uorb) * (2.d0 - occ_act(t) - occ_act(u)) 
   enddo

  enddo
 enddo

 ! equation B3.e of the annex B of Roos et. al. Chemical Physics 48 (1980) 157-173
 do t = 1, n_act_orb
  torb = list_act(t)
  do a = 1, n_virt_orb
   aorb = list_virt(a)
   super_ci_dm(aorb,torb) += c0 * lowest_super_ci_coef_mo(aorb,torb) * occ_act(t)
   super_ci_dm(torb,aorb) += c0 * lowest_super_ci_coef_mo(aorb,torb) * occ_act(t)
   do i = 1, n_core_inact_orb
    iorb = list_core_inact(i)
    super_ci_dm(aorb,torb) += lowest_super_ci_coef_mo(iorb,aorb) * lowest_super_ci_coef_mo(iorb,torb) * (2.d0 - occ_act(t))
    super_ci_dm(torb,aorb) += lowest_super_ci_coef_mo(iorb,aorb) * lowest_super_ci_coef_mo(iorb,torb) * (2.d0 - occ_act(t))
   enddo
  enddo
 enddo

 ! equation B3.f of the annex B of Roos et. al. Chemical Physics 48 (1980) 157-173
 do a = 1, n_virt_orb
  aorb = list_virt(a)
  do b = 1, n_virt_orb
   borb= list_virt(b)

   ! First term of equation B3.f 
   do i = 1, n_core_inact_orb
    iorb = list_core_inact(i)
    super_ci_dm(borb,aorb) += 2.d0 * lowest_super_ci_coef_mo(iorb,aorb) * lowest_super_ci_coef_mo(iorb,borb) 
   enddo

   ! Second term of equation B3.f
   do t = 1, n_act_orb
    torb = list_act(t)
    super_ci_dm(borb,aorb) += lowest_super_ci_coef_mo(torb,aorb) * lowest_super_ci_coef_mo(torb,borb) * occ_act(t)
   enddo
  enddo
 enddo

 END_PROVIDER 

 BEGIN_PROVIDER [double precision, superci_natorb, (ao_num,mo_num)
&BEGIN_PROVIDER [double precision, superci_nat_occ, (mo_num)
 implicit none
 call general_mo_coef_new_as_svd_vectors_of_mo_matrix_eig(super_ci_dm,mo_num,mo_num,mo_num,NatOrbsFCI,superci_nat_occ,superci_natorb)

END_PROVIDER 

 BEGIN_PROVIDER [double precision, mat_tmp_dm_super_ci, (n_act_orb,n_act_orb)]
 implicit none
 BEGIN_DOC
 ! computation of the term in [ ] in the equation B3.d of Roos et. al. Chemical Physics 48 (1980) 157-173
 !
 ! !!!!! WARNING !!!!!! there is a TYPO: a MINUS SIGN SHOULD APPEAR in that term 
 END_DOC
 integer :: a,aorb,i,iorb
 integer :: x,xorb,v,vorb
 mat_tmp_dm_super_ci = 0.d0 
 do v = 1, n_act_orb
  vorb = list_act(v)
  do x = 1, n_act_orb
   xorb = list_act(x)
   do a = 1, n_virt_orb
    aorb = list_virt(a)
    mat_tmp_dm_super_ci(x,v) += lowest_super_ci_coef_mo(aorb,vorb) * lowest_super_ci_coef_mo(aorb,xorb)
   enddo
   do i = 1, n_core_inact_orb
    iorb = list_core_inact(i)
    ! MARK THE MINUS SIGN HERE !!!!!!!!!!! BECAUSE OF TYPO IN THE ORIGINAL PAPER 
    mat_tmp_dm_super_ci(x,v) -= lowest_super_ci_coef_mo(iorb,vorb) * lowest_super_ci_coef_mo(iorb,xorb)
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
  lowest_super_ci_coef_mo(iorb,jorb) = SXeigenvec(i+1,1)
  lowest_super_ci_coef_mo(jorb,iorb) = SXeigenvec(i+1,1)
 enddo

 ! a_{it} of the equation B.2 of Roos et. al. Chemical Physics 48 (1980) 157-173
 do i = 1, n_core_inact_orb
  iorb = list_core_inact(i)
  do t = 1, n_act_orb
   torb = list_act(t)
   lowest_super_ci_coef_mo(torb,iorb) *= (2.d0 - occ_act(t))**(-0.5d0)
   lowest_super_ci_coef_mo(iorb,torb) *= (2.d0 - occ_act(t))**(-0.5d0)
  enddo
 enddo

 ! a_{ia} of the equation B.2 of Roos et. al. Chemical Physics 48 (1980) 157-173
 do i = 1, n_core_inact_orb
  iorb = list_core_inact(i)
  do a = 1, n_virt_orb
   aorb = list_virt(a)
   lowest_super_ci_coef_mo(aorb,iorb) *= sqrt2
   lowest_super_ci_coef_mo(iorb,aorb) *= sqrt2
  enddo
 enddo
 
 ! a_{ta} of the equation B.2 of Roos et. al. Chemical Physics 48 (1980) 157-173
 do a = 1, n_virt_orb
  aorb = list_virt(a)
  do t = 1, n_act_orb
   torb = list_act(t)
   lowest_super_ci_coef_mo(torb,aorb) *= occ_act(t)**(-0.5d0)
   lowest_super_ci_coef_mo(aorb,torb) *= occ_act(t)**(-0.5d0)
  enddo
 enddo

 END_PROVIDER 

