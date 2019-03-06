 BEGIN_PROVIDER [double precision, potential_x_alpha_ao,(ao_num,ao_num,N_states)]
 &BEGIN_PROVIDER [double precision, potential_x_beta_ao ,(ao_num,ao_num,N_states)]
  implicit none
  BEGIN_DOC
 ! general providers for the alpha/beta exchange potentials on the AO basis
  END_DOC
 
 BEGIN_SHELL [ /usr/bin/env python ]
import os
import glob
from qp_path import QP_SRC
funcdir=QP_SRC+'/functionals/'
os.chdir(funcdir)
functionals = map(lambda x : x.replace(".irp.f",""), glob.glob("*.irp.f"))

prefix = ""
for f in functionals:
  print """
  %sif (trim(exchange_functional) == '%s') then
    potential_x_alpha_ao = potential_x_alpha_ao_%s
    potential_x_beta_ao  = potential_x_beta_ao_%s"""%(prefix, f, f, f)
  prefix = "else "
print """
  else
   print*, 'exchange functional required does not exist ...'
   print*,'exchange_functional ',exchange_functional
   stop"""
print "endif"

 END_SHELL
 
 
 END_PROVIDER
 
 
 
  BEGIN_PROVIDER [double precision, potential_c_alpha_ao,(ao_num,ao_num,N_states)]
 &BEGIN_PROVIDER [double precision, potential_c_beta_ao,(ao_num,ao_num,N_states)]
  implicit none
  BEGIN_DOC
 ! general providers for the alpha/beta correlation potentials on the AO basis
  END_DOC
 
 BEGIN_SHELL [ /usr/bin/env python ]
import os
import glob
from qp_path import QP_SRC
funcdir=QP_SRC+'/functionals/'
os.chdir(funcdir)
functionals = map(lambda x : x.replace(".irp.f",""), glob.glob("*.irp.f"))

prefix = ""
for f in functionals:
  print """
  %sif (trim(correlation_functional) == '%s') then
    potential_c_alpha_ao = potential_c_alpha_ao_%s
    potential_c_beta_ao  = potential_c_beta_ao_%s"""%(prefix, f, f, f)
  prefix = "else "

print """
  else
   print*, 'Correlation functional required does not exist ...'
   print*,'correlation_functional ',correlation_functional
   stop"""
print "endif"

 END_SHELL
 
 END_PROVIDER
 
 
 
 
 
  BEGIN_PROVIDER [double precision, potential_x_alpha_mo,(mo_num,mo_num,N_states)]
 &BEGIN_PROVIDER [double precision, potential_x_beta_mo ,(mo_num,mo_num,N_states)]
  implicit none
  BEGIN_DOC
 ! general providers for the alpha/beta exchange potentials on the MO basis
  END_DOC
  integer :: istate
  do istate = 1, N_states
     call ao_to_mo(                                                   &
         potential_x_alpha_ao(1,1,istate),                                 &
         size(potential_x_alpha_ao,1),                                &
         potential_x_alpha_mo(1,1,istate),                                 &
         size(potential_x_alpha_mo,1)                                 &
         )
 
     call ao_to_mo(                                                   &
         potential_x_beta_ao(1,1,istate),                                  &
         size(potential_x_beta_ao,1),                                 &
         potential_x_beta_mo(1,1,istate),                                  &
         size(potential_x_beta_mo,1)                                  &
         )
  enddo
 
 END_PROVIDER
 
  BEGIN_PROVIDER [double precision, potential_c_alpha_mo,(mo_num,mo_num,N_states)]
 &BEGIN_PROVIDER [double precision, potential_c_beta_mo, (mo_num,mo_num,N_states)]
  implicit none
  BEGIN_DOC
 ! general providers for the alpha/beta correlation potentials on the MO basis
  END_DOC
  integer :: istate
  do istate = 1, N_states
     call ao_to_mo(                                                   &
         potential_c_alpha_ao(1,1,istate),                                 &
         size(potential_c_alpha_ao,1),                                &
         potential_c_alpha_mo(1,1,istate),                                 &
         size(potential_c_alpha_mo,1)                                 &
         )
 
     call ao_to_mo(                                                   &
         potential_c_beta_ao(1,1,istate),                                  &
         size(potential_c_beta_ao,1),                                 &
         potential_c_beta_mo(1,1,istate),                                  &
         size(potential_c_beta_mo,1)                                  &
         )
  enddo
 
 END_PROVIDER


  BEGIN_PROVIDER [double precision, Trace_v_xc, (N_states)]
 &BEGIN_PROVIDER [double precision, Trace_v_H, (N_states)]
 &BEGIN_PROVIDER [double precision, Trace_v_Hxc, (N_states)]
  implicit none
  integer :: i,j,istate
  double precision :: dm
  BEGIN_DOC
 ! Trace_v_xc  = \sum_{i,j} (rho_{ij}_\alpha v^{xc}_{ij}^\alpha  + rho_{ij}_\beta v^{xc}_{ij}^\beta)
 ! Trace_v_Hxc = \sum_{i,j} v^{H}_{ij} (rho_{ij}_\alpha + rho_{ij}_\beta)
 ! Trace_v_Hxc = \sum_{i,j} rho_{ij} v^{Hxc}_{ij}
  END_DOC
  do istate = 1, N_states
   Trace_v_xc(istate) = 0.d0
   Trace_v_H(istate) = 0.d0
   do i = 1, mo_num
    do j = 1, mo_num
      Trace_v_xc(istate) += (potential_x_alpha_mo(j,i,istate) + potential_c_alpha_mo(j,i,istate)) * one_e_dm_mo_alpha_for_dft(j,i,istate)
      Trace_v_xc(istate) += (potential_x_beta_mo(j,i,istate)  + potential_c_beta_mo(j,i,istate) ) * one_e_dm_mo_beta_for_dft(j,i,istate)
      dm = one_e_dm_mo_alpha_for_dft(j,i,istate) + one_e_dm_mo_beta_for_dft(j,i,istate)
      Trace_v_H(istate) += dm * short_range_Hartree_operator(j,i,istate)
    enddo
   enddo
   Trace_v_Hxc(istate) = Trace_v_xc(istate) + Trace_v_H(istate)
  enddo
 
 END_PROVIDER
 
  BEGIN_PROVIDER [double precision, Trace_v_xc_new, (N_states)]
  implicit none
  integer :: i,j,istate
  double precision :: dm
  BEGIN_DOC
 ! Trace_v_xc  = \sum_{i,j} (rho_{ij}_\alpha v^{xc}_{ij}^\alpha  + rho_{ij}_\beta v^{xc}_{ij}^\beta)
  END_DOC
  do istate = 1, N_states
   Trace_v_xc_new(istate) = 0.d0
   do i = 1, mo_num
    do j = 1, mo_num
      Trace_v_xc_new(istate) += (potential_xc_alpha_mo(j,i,istate) ) * one_e_dm_mo_alpha_for_dft(j,i,istate)
      Trace_v_xc_new(istate) += (potential_xc_beta_mo(j,i,istate)  ) * one_e_dm_mo_beta_for_dft(j,i,istate)
    enddo
   enddo
  enddo
 
 END_PROVIDER
 
  BEGIN_PROVIDER [double precision, potential_xc_alpha_mo,(mo_num,mo_num,N_states)]
 &BEGIN_PROVIDER [double precision, potential_xc_beta_mo,(mo_num,mo_num,N_states)]
  implicit none
  integer :: istate
 
  do istate = 1, N_states
     call ao_to_mo(                                                   &
         potential_xc_alpha_ao(1,1,istate),                                 &
         size(potential_xc_alpha_ao,1),                                &
         potential_xc_alpha_mo(1,1,istate),                                 &
         size(potential_xc_alpha_mo,1)                                 &
         )
 
     call ao_to_mo(                                                   &
         potential_xc_beta_ao(1,1,istate),                                  &
         size(potential_xc_beta_ao,1),                                 &
         potential_xc_beta_mo(1,1,istate),                                  &
         size(potential_xc_beta_mo,1)                                  &
         )
  enddo
 
 END_PROVIDER
 
 
  BEGIN_PROVIDER [double precision, potential_xc_alpha_ao,(ao_num,ao_num,N_states)]
 &BEGIN_PROVIDER [double precision, potential_xc_beta_ao,(ao_num,ao_num,N_states)]
  implicit none
  BEGIN_DOC
 ! general providers for the alpha/beta exchange/correlation potentials on the AO basis
  END_DOC
 
 BEGIN_SHELL [ /usr/bin/env python ]
import os
import glob
from qp_path import QP_SRC
funcdir=QP_SRC+'/functionals/'
os.chdir(funcdir)
functionals = map(lambda x : x.replace(".irp.f",""), glob.glob("*.irp.f"))

prefix = ""
for f in functionals:
  print """
  %sif (trim(exchange_functional) == '%s') then
    potential_xc_alpha_ao = potential_xc_alpha_ao_%s
    potential_xc_beta_ao  = potential_xc_beta_ao_%s"""%(prefix, f, f, f)
  prefix = "else "
print """
  else
   print*, 'exchange functional required does not exist ...'
   print*,'exchange_functional ',exchange_functional
   stop"""
print "endif"

END_SHELL
 
 END_PROVIDER

