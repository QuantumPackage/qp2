
 BEGIN_PROVIDER [double precision, potential_new_functional_x_alpha_ao,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_new_functional_x_beta_ao,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_new_functional_c_alpha_ao,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_new_functional_c_beta_ao,(ao_num,ao_num,N_states)]
 implicit none
 BEGIN_DOC
! define here your exchange/correlation potentials for alpha/beta electrons 
 END_DOC
   potential_new_functional_x_alpha_ao = 0.d0 ! replace by your new provider 
   potential_new_functional_x_beta_ao  = 0.d0 ! replace by your new provider 

   potential_new_functional_c_alpha_ao = 0.d0 ! replace by your new provider 
   potential_new_functional_c_beta_ao  = 0.d0 ! replace by your new provider 


END_PROVIDER

