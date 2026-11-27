use sizet

 BEGIN_PROVIDER [double precision, mo_deriv_1_x , (mo_num,mo_num)]
&BEGIN_PROVIDER [double precision, mo_deriv_1_y , (mo_num,mo_num)]
&BEGIN_PROVIDER [double precision, mo_deriv_1_z , (mo_num,mo_num)]
 BEGIN_DOC
 ! array of the integrals of MO_i * d/dx MO_j
 ! array of the integrals of MO_i * d/dy MO_j
 ! array of the integrals of MO_i * d/dz MO_j
 END_DOC
 implicit none

  call ao_to_mo(                                                     &
      ao_deriv_1_x,                                                   &
      int(size(ao_deriv_1_x,1), size_t),                              &
      mo_deriv_1_x,                                                   &
      int(size(mo_deriv_1_x,1), size_t)                               &
      )
  call ao_to_mo(                                                     &
      ao_deriv_1_y,                                                   &
      int(size(ao_deriv_1_y,1), size_t),                              &
      mo_deriv_1_y,                                                   &
      int(size(mo_deriv_1_y,1), size_t)                               &
      )
  call ao_to_mo(                                                     &
      ao_deriv_1_z,                                                   &
      int(size(ao_deriv_1_z,1), size_t),                              &
      mo_deriv_1_z,                                                   &
      int(size(mo_deriv_1_z,1), size_t)                               &
      )

END_PROVIDER
