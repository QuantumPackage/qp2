BEGIN_PROVIDER [ double precision, ao_extra_center]
 implicit none
 ao_extra_center = 0.01d0
END_PROVIDER 

 BEGIN_PROVIDER [ integer, n_func_tot]
 implicit none
 BEGIN_DOC
 ! n_func_tot :: total number of functions in the fitted basis set
 ! 
 ! returned in an uncontracted way
 END_DOC
 integer :: i,prefact
 n_func_tot = 0
   print*,'n_func_tot '
 do i = 1, ao_num
   if(ao_l(i) == 0)then
    prefact = 1 ! s functions 
   else
    ! p functions are fitted with 2 functions
    ! d functions are fitted with 4 functions etc ...
    prefact=2*ao_l(i) 
   endif
   n_func_tot += prefact * ao_prim_num(i) 
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ integer, n_prim_tot_orig]
 implicit none
 integer :: i
 n_prim_tot_orig = 0
 do i = 1, ao_num
  n_prim_tot_orig += ao_prim_num(i) 
 enddo
END_PROVIDER 


BEGIN_PROVIDER [ logical, lmax_too_big]
 implicit none
 if (ao_l_max.gt.1)then
  lmax_too_big = .True.
 else
  lmax_too_big = .False.
 endif
 if(lmax_too_big)then
  print*,'STOPPING !! lmax is larger than 1 !'
  print*,'Cannot yet fit with 1s functions ...'
  stop
 endif
END_PROVIDER 

 BEGIN_PROVIDER [ integer, n_2p_func_orig]
&BEGIN_PROVIDER [ integer, n_2p_func_tot]
 implicit none
 integer :: i
 BEGIN_DOC
 ! n_2p_func_orig :: number of 2p functions in the original basis
 !
 ! n_2p_func_tot :: total number of p functions in the fitted basis
 END_DOC
 n_2p_func_orig= 0
 n_2p_func_tot = 0
 do i = 1, ao_num
  if(ao_l(i)==1)then
   n_2p_func_orig+= 1
   n_2p_func_tot += ao_prim_num(i) * 2
  endif
 enddo
 print*,'n_2p_func_tot = ',n_2p_func_tot
END_PROVIDER 

BEGIN_PROVIDER [ integer, list_2p_functions, (n_2p_func_orig)]
 implicit none
 BEGIN_DOC
 ! list of 2p functions in the original basis
 END_DOC
 integer :: i,j
 j=0
 do i = 1, ao_num
  if(ao_l(i)==1)then
   j+=1
   list_2p_functions(j) = i
  endif
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ integer, extra_fictious_nucl]
 implicit none
 extra_fictious_nucl = n_2p_func_tot
END_PROVIDER 

BEGIN_PROVIDER [ integer, new_nucl_num]
 implicit none
 new_nucl_num = nucl_num + n_2p_func_tot
 print*,'new_nucl_num = ',new_nucl_num
END_PROVIDER 

 BEGIN_PROVIDER [ character*(32), new_nucl_label_1s , (new_nucl_num) ]
&BEGIN_PROVIDER [ integer, list_real_nucl, (nucl_num) ]
&BEGIN_PROVIDER [ integer, list_fict_nucl, (extra_fictious_nucl) ]
 implicit none
 integer :: i,j
 do i = 1, nucl_num 
  new_nucl_label_1s(i) = nucl_label(i)
  list_real_nucl(i) = i
 enddo
 j=0
 do i = nucl_num+1,new_nucl_num
  j+=1
  new_nucl_label_1s(i) = "X"
  list_fict_nucl(j) = i
 enddo
END_PROVIDER 
 
 BEGIN_PROVIDER [ double precision,  new_nucl_coord_1s, (new_nucl_num,3)]
 implicit none
 integer :: i,j
 do i = 1, new_nucl_num
  new_nucl_coord_1s(i,1:3) = new_nucl_coord_1s_transp(1:3,i)
 enddo
 END_PROVIDER 

 BEGIN_PROVIDER [ double precision, new_nucl_coord_1s_transp, (3,new_nucl_num)]
&BEGIN_PROVIDER [ double precision, new_nucl_charge_1s, (new_nucl_num)]
&BEGIN_PROVIDER [ integer, extra_nucl_real_fictious_list_prov, (extra_fictious_nucl)]
 implicit none
 BEGIN_DOC
! the real atoms are located in the first nucl_num entries 
! 
! then the fictious atoms are located after
 END_DOC
 integer :: i,ii,j,i_ao,k,n_ao
 integer :: return_xyz_int,power(3),good_i
 new_nucl_charge_1s = 0.d0
 do i = 1, nucl_num
  new_nucl_coord_1s_transp(1:3,i) = nucl_coord_transp(1:3,i)
  new_nucl_charge_1s(i) = nucl_charge(i)
 enddo
 k = nucl_num
 do i = 1, nucl_num
  do ii = 1, Nucl_N_Aos(i)
   i_ao = nucl_aos_transposed(ii,i)
   if(ao_l(i_ao)==1)then
    ! split the function into 2 s functions 
    ! one is centered in R_x + d 
    power(1:3) = ao_power(i_ao,1:3)
    good_i = return_xyz_int(power)
    do j = 1, ao_prim_num(i_ao)
     k+=1
     new_nucl_coord_1s_transp(1:3,k)= nucl_coord_transp(1:3,i)
     new_nucl_coord_1s_transp(good_i,k)+= ao_extra_center
     new_nucl_charge_1s(k) = 0.d0
     extra_nucl_real_fictious_list_prov(k-nucl_num)=i
     k+=1
     ! one is centered in R_x - d 
     new_nucl_coord_1s_transp(1:3,k)= nucl_coord_transp(1:3,i)
     new_nucl_coord_1s_transp(good_i,k)-= ao_extra_center
     new_nucl_charge_1s(k) = 0.d0
     extra_nucl_real_fictious_list_prov(k-nucl_num)=i
    enddo
   else if(ao_l(i_ao).gt.1)then
    print*,'WARNING ! Lmax value not implemented yet !'
    print*,'stopping ...'
    stop
   endif
  enddo
 enddo

END_PROVIDER 

 BEGIN_PROVIDER [ integer, new_n_AOs_max]
 implicit none
 new_n_AOs_max = ao_prim_num_max * n_AOs_max
 
 END_PROVIDER 


 BEGIN_PROVIDER [ integer, new_Nucl_N_Aos, (new_nucl_num)]
&BEGIN_PROVIDER [ integer, new_nucl_aos_transposed, (new_n_AOs_max,new_nucl_num) ]
&BEGIN_PROVIDER [ double precision, new_ao_expo_1s , (n_func_tot) ]
&BEGIN_PROVIDER [ integer, new_ao_nucl_1s, (n_func_tot)]
 implicit none
 integer :: i,j,ii,i_ao,n_func,n_func_total,n_nucl
 double precision :: coef
 n_func_total = 0
 do i = 1, nucl_num
  n_func = 0
  do ii = 1, Nucl_N_Aos(i)
   i_ao = nucl_aos_transposed(ii,i)
   if(ao_l(i_ao)==0)then
    do j = 1, ao_prim_num(i_ao)
     coef= ao_expo(i_ao,j)
     n_func_total += 1
     n_func +=1
     new_nucl_aos_transposed(n_func,i) = n_func_total
     new_ao_expo_1s(n_func_total) = coef
     new_ao_nucl_1s(n_func_total) = i
    enddo
   endif
  enddo
  new_Nucl_N_Aos(i) = n_func
 enddo
 n_nucl=nucl_num
 do i = 1, nucl_num
  do ii = 1, Nucl_N_Aos(i)
   i_ao = nucl_aos_transposed(ii,i)
   if(ao_l(i_ao)==1)then
    do j = 1, ao_prim_num(i_ao)
     coef= ao_expo(i_ao,j)
     n_func_total+=1
     n_nucl +=1
     new_nucl_aos_transposed(1,n_nucl) = n_func_total
     new_ao_expo_1s(n_func_total) = coef
     new_Nucl_N_Aos(n_nucl)=1
     new_ao_nucl_1s(n_func_total) = n_nucl

     n_func_total+=1
     n_nucl +=1
     new_nucl_aos_transposed(1,n_nucl) = n_func_total
     new_ao_expo_1s(n_func_total) = coef
     new_Nucl_N_Aos(n_nucl)=1
     new_ao_nucl_1s(n_func_total) = n_nucl
    enddo
   endif
  enddo
 enddo

END_PROVIDER 

BEGIN_PROVIDER [ double precision, new_ao_coef_1s , (n_func_tot) ]
  implicit none
 integer :: i
  BEGIN_DOC
! Primitive coefficients, read from input. Those should not be used directly, as the MOs are expressed on the basis of **normalized** AOs.
  END_DOC
  do i = 1, n_func_tot
   new_ao_coef_1s(i) = 1.d0
  enddo
END_PROVIDER 

BEGIN_PROVIDER [ integer, new_ao_prim_num_1s, (n_func_tot)]
 implicit none
 integer :: i
 do i = 1, n_func_tot
  new_ao_prim_num_1s(i) = 1
 enddo
END_PROVIDER 

BEGIN_PROVIDER [integer, new_ao_power_1s, (n_func_tot,3)]
 implicit none
 new_ao_power_1s = 0
END_PROVIDER 

integer function return_xyz_int(power)
 implicit none
 integer,intent(in) :: power(3)
 if(power(1) == 1 .and. power(2) ==0 .and. power(3) ==0)then
  return_xyz_int = 1
 else if (power(2) == 1 .and. power(1) ==0 .and. power(3) ==0)then
  return_xyz_int = 2
 else if (power(3) == 1 .and. power(1) ==0 .and. power(2) ==0)then
  return_xyz_int = 3
 else 
  return_xyz_int = -1000
 endif
end
