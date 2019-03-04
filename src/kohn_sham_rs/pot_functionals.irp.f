 BEGIN_PROVIDER [double precision, ao_potential_alpha_xc, (ao_num, ao_num)]
&BEGIN_PROVIDER [double precision, ao_potential_beta_xc, (ao_num, ao_num)]
 implicit none
 integer :: i,j,k,l
 ao_potential_alpha_xc = 0.d0
 ao_potential_beta_xc = 0.d0
 !if(same_xc_func)then
 ! print*,'coucou '
 ! do i = 1, ao_num
 !  do j = 1, ao_num
 !   ao_potential_alpha_xc(j,i) =  potential_xc_alpha_ao(j,i,1) 
 !   ao_potential_beta_xc(j,i)  =  potential_xc_beta_ao(j,i,1)  
 !  enddo
 ! enddo
 !else
   print*,'coucou 2'
   do i = 1, ao_num
    do j = 1, ao_num
     ao_potential_alpha_xc(j,i) =  potential_c_alpha_ao(j,i,1) + potential_x_alpha_ao(j,i,1)
     ao_potential_beta_xc(j,i)  =  potential_c_beta_ao(j,i,1)  + potential_x_beta_ao(j,i,1)
    enddo
   enddo
 !endif

END_PROVIDER

BEGIN_PROVIDER [double precision, e_exchange_dft]
 implicit none
  print*,'energy_x = ',energy_x
  e_exchange_dft = energy_x(1)

END_PROVIDER

BEGIN_PROVIDER [double precision, e_correlation_dft]
 implicit none
  print*,'energy_c = ',energy_c
  e_correlation_dft = energy_c(1)

END_PROVIDER
