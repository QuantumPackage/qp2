BEGIN_PROVIDER[double precision, energy_x_none, (N_states) ]
  implicit none
  BEGIN_DOC
  ! null exchange energy 
  END_DOC
  energy_x_none = 0.d0
END_PROVIDER 


BEGIN_PROVIDER[double precision, energy_c_none, (N_states) ]
  implicit none
  BEGIN_DOC
  ! null correlation energy 
  END_DOC
  energy_c_none = 0.d0
END_PROVIDER 


BEGIN_PROVIDER [double precision, potential_x_alpha_ao_none,(ao_num,ao_num,N_states)]
 implicit none
 BEGIN_DOC
 ! Potential for a null exchange-correlation functional
 END_DOC
 potential_x_alpha_ao_none = 0.d0 
END_PROVIDER


BEGIN_PROVIDER [double precision, potential_x_beta_ao_none,(ao_num,ao_num,N_states)]
 implicit none
 BEGIN_DOC
 ! Potential for a null exchange-correlation functional
 END_DOC
 potential_x_beta_ao_none  = 0.d0 
END_PROVIDER


BEGIN_PROVIDER [double precision, potential_c_alpha_ao_none,(ao_num,ao_num,N_states)]
 implicit none
 BEGIN_DOC
 ! Potential for a null exchange-correlation functional
 END_DOC
 potential_c_alpha_ao_none = 0.d0 
END_PROVIDER


BEGIN_PROVIDER [double precision, potential_c_beta_ao_none,(ao_num,ao_num,N_states)]
 implicit none
 BEGIN_DOC
 ! Potential for a null exchange-correlation functional
 END_DOC
 potential_c_beta_ao_none  = 0.d0 
END_PROVIDER


BEGIN_PROVIDER [double precision, potential_xc_alpha_ao_none ,(ao_num,ao_num,N_states)]
 implicit none
 BEGIN_DOC
 ! Potential for a null exchange-correlation functional
 END_DOC
 potential_xc_alpha_ao_none  = 0.d0 
END_PROVIDER


BEGIN_PROVIDER [double precision, potential_xc_beta_ao_none ,(ao_num,ao_num,N_states)]
 implicit none
 BEGIN_DOC
 ! Potential for a null exchange-correlation functional
 END_DOC
 potential_xc_beta_ao_none  = 0.d0 
END_PROVIDER

