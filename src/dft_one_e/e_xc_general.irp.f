
BEGIN_PROVIDER [double precision, energy_x, (N_states)]
 implicit none
 BEGIN_DOC
 ! correlation energies general providers.
 END_DOC
  if(trim(exchange_functional)=="short_range_LDA")then
   energy_x = energy_sr_x_LDA
  else if(exchange_functional.EQ."short_range_PBE")then
   energy_x = energy_sr_x_PBE
  else if(exchange_functional.EQ."PBE")then
   energy_x = energy_sr_x_PBE
  else if(exchange_functional.EQ."LDA")then
   energy_x = energy_sr_x_LDA
  else if(exchange_functional.EQ."my_functional")then
   energy_x = energy_x_new_functional
  else
   print*, 'Exchange functional required does not exist ...'
   print*,'exchange_functional',exchange_functional
   stop
  endif

BEGIN_SHELL [ /usr/bin/env python ]
import os
functionals = map(lambda x : x.replace(".irp.f",""), os.listdir("../dft_utils_one_e/functionals/"))

prefix = ""
for f in functionals:
  print """
  %sif (trim(exchange_functional) == '%s') then
    energy_x = energy_x_%s"""%(prefix, f, f)
  prefix = "else " 
print "endif"

END_SHELL

END_PROVIDER




BEGIN_PROVIDER [double precision, energy_c, (N_states)]
 implicit none
 BEGIN_DOC
 ! correlation and exchange energies general providers.
 END_DOC
  if(trim(correlation_functional)=="short_range_LDA")then
   energy_c = energy_sr_c_LDA
  else if(correlation_functional.EQ."short_range_PBE")then
   energy_c = energy_sr_c_PBE
  else if(correlation_functional.EQ."PBE")then
   energy_c = energy_sr_c_PBE
  else if(correlation_functional.EQ."LDA")then
   energy_c = energy_sr_c_LDA
  else if(correlation_functional.EQ."my_functional")then
   energy_c = energy_c_new_functional
  else
   print*, 'Correlation functional required does not ecist ...'
   print*,'correlation_functional',correlation_functional
   stop
  endif

BEGIN_SHELL [ /usr/bin/env python ]
import os
functionals = map(lambda x : x.replace(".irp.f",""), os.listdir("../dft_utils_one_e/functionals/"))

prefix = ""
for f in functionals:
  print """
  %sif (trim(correlation_functional) == '%s') then
    energy_c = energy_c_%s"""%(prefix, f, f)
  prefix = "else " 
print "endif"

END_SHELL

END_PROVIDER
