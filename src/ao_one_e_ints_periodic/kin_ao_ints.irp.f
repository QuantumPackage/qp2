BEGIN_PROVIDER [complex*16, ao_kinetic_integrals, (ao_num,ao_num)]
  implicit none
  BEGIN_DOC
  ! array of the priminitve basis kinetic integrals
  !  \langle \chi_i |\hat{T}| \chi_j \rangle
  END_DOC
  
  if (read_ao_integrals_kinetic) then
    call read_one_e_integrals_complex('ao_kinetic_integrals', ao_kinetic_integrals,&
        size(ao_kinetic_integrals,1), size(ao_kinetic_integrals,2))
    print *,  'AO kinetic integrals read from disk'
  else
    print *, 'complex AO kinetic integrals must be provided'
    stop
  endif
END_PROVIDER


