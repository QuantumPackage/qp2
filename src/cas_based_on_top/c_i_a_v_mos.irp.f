
 BEGIN_PROVIDER[double precision, core_mos_in_r_array, (n_core_orb,n_points_final_grid)]
&BEGIN_PROVIDER[double precision, core_mos_in_r_array_transp,(n_points_final_grid,n_core_orb)]
 implicit none
 BEGIN_DOC
! all COREE  MOs on the grid points, arranged in two different ways
 END_DOC
 integer :: i,j,k
 do i = 1, n_core_orb
  j = list_core(i) 
  do k = 1, n_points_final_grid
   core_mos_in_r_array_transp(k,i) = mos_in_r_array_transp(k,j)
  enddo
 enddo

 do k = 1, n_points_final_grid
  do i = 1, n_core_orb
   core_mos_in_r_array(i,k) = core_mos_in_r_array_transp(k,i)
  enddo
 enddo

END_PROVIDER 


 BEGIN_PROVIDER[double precision, inact_mos_in_r_array, (n_inact_orb,n_points_final_grid)]
&BEGIN_PROVIDER[double precision, inact_mos_in_r_array_transp,(n_points_final_grid,n_inact_orb)]
 implicit none
 BEGIN_DOC
! all INACTIVE MOs on the grid points, arranged in two different ways
 END_DOC
 integer :: i,j,k
 do i = 1, n_inact_orb
  j = list_inact(i) 
  do k = 1, n_points_final_grid
   inact_mos_in_r_array_transp(k,i) = mos_in_r_array_transp(k,j)
  enddo
 enddo

 do k = 1, n_points_final_grid
  do i = 1, n_inact_orb
   inact_mos_in_r_array(i,k) = inact_mos_in_r_array_transp(k,i)
  enddo
 enddo

END_PROVIDER 

 BEGIN_PROVIDER[double precision, act_mos_in_r_array, (n_act_orb,n_points_final_grid)]
&BEGIN_PROVIDER[double precision, act_mos_in_r_array_transp,(n_points_final_grid,n_act_orb)]
 implicit none
 BEGIN_DOC
! all ACTIVE MOs on the grid points, arranged in two different ways
 END_DOC
 integer :: i,j,k
 do i = 1, n_act_orb
  j = list_act(i) 
  do k = 1, n_points_final_grid
   act_mos_in_r_array_transp(k,i) = mos_in_r_array_transp(k,j)
  enddo
 enddo

 do k = 1, n_points_final_grid
  do i = 1, n_act_orb
   act_mos_in_r_array(i,k) = act_mos_in_r_array_transp(k,i)
  enddo
 enddo

END_PROVIDER 


 BEGIN_PROVIDER[double precision, virt_mos_in_r_array, (n_virt_orb,n_points_final_grid)]
&BEGIN_PROVIDER[double precision, virt_mos_in_r_array_transp,(n_points_final_grid,n_virt_orb)]
 implicit none
 BEGIN_DOC
! all VIRTUAL MOs on the grid points, arranged in two different ways
 END_DOC
 integer :: i,j,k
 do i = 1, n_virt_orb
  j = list_virt(i) 
  do k = 1, n_points_final_grid
   virt_mos_in_r_array_transp(k,i) = mos_in_r_array_transp(k,j)
  enddo
 enddo

 do k = 1, n_points_final_grid
  do i = 1, n_virt_orb
   virt_mos_in_r_array(i,k) = virt_mos_in_r_array_transp(k,i)
  enddo
 enddo

END_PROVIDER 

 BEGIN_PROVIDER[double precision, core_inact_act_mos_in_r_array, (n_core_inact_act_orb,n_points_final_grid)]
&BEGIN_PROVIDER[double precision, core_inact_act_mos_in_r_array_transp,(n_points_final_grid,n_core_inact_act_orb)]
 implicit none
 integer :: i,j,k
 do i = 1, n_core_inact_act_orb
  j = list_core_inact_act(i) 
  do k = 1, n_points_final_grid
   core_inact_act_mos_in_r_array_transp(k,i) = mos_in_r_array_transp(k,j)
  enddo
 enddo

 do k = 1, n_points_final_grid
  do i = 1, n_core_inact_act_orb
   core_inact_act_mos_in_r_array(i,k) = core_inact_act_mos_in_r_array_transp(k,i)
  enddo
 enddo

END_PROVIDER 

 BEGIN_PROVIDER[double precision, core_inact_act_mos_grad_in_r_array, (3,n_core_inact_act_orb,n_points_final_grid)]
 implicit none
 integer :: i,j,k,l
 do i = 1, n_core_inact_act_orb
  j = list_core_inact_act(i) 
  do k = 1, n_points_final_grid
   do l = 1, 3
    core_inact_act_mos_grad_in_r_array(l,i,k) = mos_grad_in_r_array(j,k,l)
   enddo
  enddo
 enddo

END_PROVIDER 


