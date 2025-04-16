
BEGIN_PROVIDER [ double precision, ao_integrals_pt_chrg, (ao_num,ao_num)]

  BEGIN_DOC
  !  Point charge-electron interaction, in the |AO| basis set.
  !
  !  :math:`\langle \chi_i | -\sum_charge charge * \frac{1}{|r-R_charge|} | \chi_j \rangle`
  !
  ! Notice the minus sign convention as it is supposed to be for electrons. 
  END_DOC

  implicit none
  call ao_cart_to_ao_basis(ao_cart_integrals_pt_chrg, ao_cart_num, ao_integrals_pt_chrg, ao_num)

END_PROVIDER
