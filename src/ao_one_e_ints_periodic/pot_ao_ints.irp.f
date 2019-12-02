BEGIN_PROVIDER [ complex*16, ao_integrals_n_e, (ao_num,ao_num)]
  BEGIN_DOC
  !  Nucleus-electron interaction, in the |AO| basis set.
  !
  !  :math:`\langle \chi_i | -\sum_A \frac{1}{|r-R_A|} | \chi_j \rangle`
  END_DOC

   if (read_ao_integrals_e_n) then
    call read_one_e_integrals_complex('ao_ne_integral', ao_integrals_n_e,      &
            size(ao_integrals_n_e,1), size(ao_integrals_n_e,2))
     print *,  'AO N-e integrals read from disk'
   else
     print *, 'complex AO N-e integrals must be provided'
     stop
   endif
END_PROVIDER

BEGIN_PROVIDER [ double precision, ao_integrals_n_e_per_atom, (ao_num,ao_num,nucl_num)]
  BEGIN_DOC
! Nucleus-electron interaction in the |AO| basis set, per atom A.
!
! :math:`\langle \chi_i | -\frac{1}{|r-R_A|} | \chi_j \rangle`
  END_DOC
  print *, 'ao_nucl_elec_integral_per_atom not implemented for k-points'
  stop

END_PROVIDER

