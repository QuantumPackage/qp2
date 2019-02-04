
 BEGIN_PROVIDER[double precision, energy_x_new_functional, (N_states) ]
&BEGIN_PROVIDER[double precision, energy_c_new_functional, (N_states) ]
 implicit none
 BEGIN_DOC
! energy_x_new_functional = define here your functional 
! energy_c_new_functional = define here your functional 
 END_DOC
  energy_x_new_functional = 0.d0 ! replace by your new provider  
  energy_c_new_functional = 0.d0 ! replace by your new provider  

END_PROVIDER 


