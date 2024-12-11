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

BEGIN_PROVIDER [ integer, new_nucl_num]
 implicit none
 new_nucl_num = n_func_tot
END_PROVIDER 

BEGIN_PROVIDER [ double precision, new_ao_expo_1s , (n_func_tot) ]
 implicit none
 integer :: i,j,ii,i_ao,k
 k = 0
 do i = 1, nucl_num
  do ii = 1, Nucl_N_Aos(i)
   i_ao = nucl_aos_transposed(ii,i)
   if(ao_l(i_ao)==0)then
    do j = 1, ao_prim_num(i_ao)
     k+=1
     new_ao_expo_1s(k)= ao_expo(i_ao,j)
    enddo
   else if(ao_l(i_ao)==1)then
    ! for 'p' functions 
    ! you replace the function by 2 functions 's' 
    do j = 1, ao_prim_num(i_ao)
     k+=1
     new_ao_expo_1s(k)= ao_expo(i_ao,j)
     k+=1
     new_ao_expo_1s(k)= ao_expo(i_ao,j)
    enddo
   else 
    print*,'WARNING ! Lmax value not implemented yet !'
    print*,'stopping ...'
    stop
   endif
  enddo
 enddo
 if(k.ne.n_func_tot)then
  print*,'pb !!! k NE n_func_tot !!'
  print*,k,n_func_tot
  stop
 endif

END_PROVIDER 

BEGIN_PROVIDER [ double precision, new_nucl_coord, (3,new_nucl_num)]
 implicit none
 integer :: i,ii,j,i_ao,k
 k = 0
 do i = 1, nucl_num
  do ii = 1, Nucl_N_Aos(i)
   i_ao = nucl_aos_transposed(ii,i)
   if(ao_l(i_ao)==0)then
    do j = 1, ao_prim_num(i_ao)
     k+=1
     new_nucl_coord(1:3,k)=nucl_coord_transp(1:3,i)
    enddo
   else if(ao_l(i_ao)==1)then
    ! split the function into 2 s functions 
    ! one is centered in R_x + d 
    do j = 1, ao_prim_num(i_ao)
     k+=1
     new_nucl_coord(2:3,k)= nucl_coord_transp(2:3,i)
     new_nucl_coord(1,k)= nucl_coord_transp(1,i)+ao_extra_center
     k+=1
     ! one is centered in R_x - d 
     new_nucl_coord(2:3,k)= nucl_coord_transp(2:3,i)
     new_nucl_coord(1,k)= nucl_coord_transp(1,i)-ao_extra_center
    enddo
   else 
    print*,'WARNING ! Lmax value not implemented yet !'
    print*,'stopping ...'
    stop
   endif
  enddo
 enddo
 if(k.ne.n_func_tot)then
  print*,'pb !!! k NE n_func_tot !!'
  print*,k,n_func_tot
  stop
 endif

END_PROVIDER 
