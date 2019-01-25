
  BEGIN_PROVIDER [double precision, energy_x, (N_states)]
 &BEGIN_PROVIDER [double precision, energy_c, (N_states)]
 implicit none
 BEGIN_DOC
 ! correlation and exchange energies general providers.
 END_DOC
  if(trim(exchange_functional)=="short_range_LDA")then
   energy_x = energy_sr_x_LDA
   energy_x = energy_sr_x_LDA
  else if(exchange_functional.EQ."short_range_PBE")then
   energy_x = energy_sr_x_PBE
   energy_x = energy_sr_x_PBE
  else if(exchange_functional.EQ."None")then
   energy_x = 0.d0
   energy_x = 0.d0
  else
   print*, 'Exchange functional required does not exist ...'
   print*,'exchange_functional',exchange_functional
   stop
  endif

  if(trim(correlation_functional)=="short_range_LDA")then
   energy_c = energy_sr_c_LDA
   energy_c = energy_sr_c_LDA
  else if(correlation_functional.EQ."short_range_PBE")then
   energy_c = energy_sr_c_PBE
   energy_c = energy_sr_c_PBE
  else if(correlation_functional.EQ."None")then
   energy_c = 0.d0
   energy_c = 0.d0
  else
   print*, 'Correlation functional required does not ecist ...'
   print*,'correlation_functional',correlation_functional
   stop
  endif

END_PROVIDER
