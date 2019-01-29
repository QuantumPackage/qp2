
 BEGIN_PROVIDER [double precision, potential_new_functional_x_alpha_ao,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_new_functional_x_beta_ao,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_new_functional_c_alpha_ao,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_new_functional_c_beta_ao,(ao_num,ao_num,N_states)]
 implicit none
 BEGIN_DOC
! define here your exchange/correlation potentials for alpha/beta electrons 
 END_DOC
 potential_new_functional_x_alpha_ao = 0.d0
 potential_new_functional_c_alpha_ao = 0.d0
 potential_new_functional_x_beta_ao = 0.d0
 potential_new_functional_c_beta_ao = 0.d0
  if(trim(new_exchange_functional)=="your_new_keyword")then
   potential_new_functional_x_alpha_ao = 0.d0 ! replace by your new provider 
   potential_new_functional_x_beta_ao  = 0.d0 ! replace by your new provider 
  else if(new_exchange_functional.EQ."None")then
   potential_new_functional_x_alpha_ao = 0.d0
   potential_new_functional_x_beta_ao = 0.d0
  else
   print*, 'Exchange functional required does not exist ...'
   print*,'new_exchange_functional',new_exchange_functional
   stop
  endif

  if(trim(new_correlation_functional)=="your_new_keyword")then
   potential_new_functional_c_alpha_ao = 0.d0 ! replace by your new provider 
   potential_new_functional_c_beta_ao  = 0.d0 ! replace by your new provider 
  else if(new_correlation_functional.EQ."None")then
   potential_new_functional_c_alpha_ao = 0.d0
   potential_new_functional_c_beta_ao = 0.d0
  else
   print*, 'Correlation functional required does not ecist ...'
   print*,'new_correlation_functional',new_correlation_functional
   stop
  endif


END_PROVIDER

