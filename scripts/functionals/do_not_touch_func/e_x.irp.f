
 BEGIN_PROVIDER[double precision, energy_x_new_functional, (N_states) ]
&BEGIN_PROVIDER[double precision, energy_c_new_functional, (N_states) ]
 implicit none
 BEGIN_DOC
! energy_x_new_functional = define here your functional 
! energy_c_new_functional = define here your functional 
 END_DOC
 energy_x_new_functional = 0.d0
 energy_c_new_functional = 0.d0
 if(trim(new_exchange_functional)=="your_new_keyword")then
  energy_x_new_functional = 0.d0 ! replace by your new provider  
 else if(new_exchange_functional.EQ."None")then
  energy_x_new_functional = 0.d0 ! replace by your new provider  
 else 
   print*, 'Exchange functional required does not exist ...'
   print*,'new_exchange_functional',new_exchange_functional
   stop
 endif

 if(trim(new_correlation_functional)=="your_new_keyword")then
  energy_c_new_functional = 0.d0  ! replace by your new provider  
 else if(new_correlation_functional.EQ."None")then
  energy_c_new_functional = 0.d0  ! replace by your new provider  
 else 
   print*, 'Correlation functional required does not exist ...'
   print*,'new_correlation_functional',new_correlation_functional
   stop
 endif

END_PROVIDER 


