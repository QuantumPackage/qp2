BEGIN_PROVIDER [double precision, energy_x, (N_states)]
  implicit none
  BEGIN_DOC
  ! correlation energies general providers.
  END_DOC
 
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
print """
  else
   print*, 'exchange functional required does not exist ...'
   print*,'exchange_functional ',exchange_functional
   stop"""
print "endif"

 END_SHELL
 
 
 END_PROVIDER
 
 
 
 
 BEGIN_PROVIDER [double precision, energy_c, (N_states)]
  implicit none
  BEGIN_DOC
  ! correlation and exchange energies general providers.
  END_DOC
 
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
  %sif (trim(correlation_functional) == '%s') then
    energy_c = energy_c_%s"""%(prefix, f, f)
  prefix = "else "

print """
  else
   print*, 'Correlation functional required does not exist ...'
   print*,'correlation_functional ',correlation_functional
   stop"""
print "endif"

 END_SHELL
 
 END_PROVIDER

