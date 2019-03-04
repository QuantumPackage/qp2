
BEGIN_PROVIDER [double precision, energy_x, (N_states)]
 implicit none
 BEGIN_DOC
 ! correlation energies general providers.
 END_DOC
  if(trim(exchange_functional)=="short_range_lda")then
   energy_x = energy_sr_x_lda
  else if(exchange_functional.EQ."short_range_pbe")then
   energy_x = energy_sr_x_pbe
  else if(exchange_functional.EQ."pbe")then
   energy_x = energy_sr_x_pbe
  else if(exchange_functional.EQ."lda")then
   energy_x = energy_sr_x_lda
  else
   print*, 'Exchange functional required does not exist ...'
   print*,'exchange_functional',exchange_functional
   stop
  endif

BEGIN_SHELL [ /usr/bin/env python ]
import os 
import glob 
import sys
qproot=os.environ['QP_ROOT']
funcdir='../functionals/'
os.chdir(funcdir)
functionals = map(lambda x : x.replace(".irp.f",""), glob.glob("*.irp.f"))
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
  if(trim(correlation_functional)=="short_range_lda")then
   energy_c = energy_sr_c_lda
  else if(correlation_functional.EQ."short_range_pbe")then
   energy_c = energy_sr_c_pbe
  else if(correlation_functional.EQ."pbe")then
   energy_c = energy_sr_c_pbe
  else if(correlation_functional.EQ."lda")then
   energy_c = energy_sr_c_lda
  else
   print*, 'Correlation functional required does not ecist ...'
   print*,'correlation_functional',correlation_functional
   stop
  endif

BEGIN_SHELL [ /usr/bin/env python ]
import os 
import glob 
import sys
qproot=os.environ['QP_ROOT']
funcdir='../functionals/'
os.chdir(funcdir)
functionals = map(lambda x : x.replace(".irp.f",""), glob.glob("*.irp.f"))
prefix = ""
for f in functionals:
  print """
  %sif (trim(exchange_functional) == '%s') then
    energy_c = energy_c_%s"""%(prefix, f, f)
  prefix = "else " 
print "endif"

END_SHELL

END_PROVIDER
