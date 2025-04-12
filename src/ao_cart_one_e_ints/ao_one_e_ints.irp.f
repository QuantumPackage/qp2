 BEGIN_PROVIDER [ double precision, ao_cart_one_e_integrals,(ao_cart_num,ao_cart_num)]
&BEGIN_PROVIDER [ double precision, ao_cart_one_e_integrals_diag,(ao_cart_num)]
  implicit none
  integer :: i,j,n,l
  BEGIN_DOC
 ! One-electron Hamiltonian in the |AO| basis.
  END_DOC

  IF (read_ao_cart_one_e_integrals) THEN
     call ezfio_get_ao_cart_one_e_ints_ao_cart_one_e_integrals(ao_cart_one_e_integrals)
  ELSE
     ao_cart_one_e_integrals = ao_cart_integrals_n_e + ao_cart_kinetic_integrals

  ENDIF

  DO j = 1, ao_cart_num
    ao_cart_one_e_integrals_diag(j) = ao_cart_one_e_integrals(j,j)
  ENDDO

  IF (write_ao_cart_one_e_integrals) THEN
       call ezfio_set_ao_cart_one_e_ints_ao_cart_one_e_integrals(ao_cart_one_e_integrals)
       print *,  'AO one-e integrals written to disk'
  ENDIF

END_PROVIDER

BEGIN_PROVIDER [ double precision, ao_cart_one_e_integrals_imag,(ao_cart_num,ao_cart_num)]
  implicit none
  integer :: i,j,n,l
  BEGIN_DOC
 ! One-electron Hamiltonian in the |AO| basis.
  END_DOC

  IF (read_ao_cart_one_e_integrals) THEN
     call ezfio_get_ao_cart_one_e_ints_ao_cart_one_e_integrals(ao_cart_one_e_integrals_imag)
  ELSE
     print *,  irp_here, ': Not yet implemented'
     stop -1
  ENDIF

  IF (write_ao_cart_one_e_integrals) THEN
       call ezfio_set_ao_cart_one_e_ints_ao_cart_one_e_integrals(ao_cart_one_e_integrals_imag)
       print *,  'AO one-e integrals written to disk'
  ENDIF

END_PROVIDER

