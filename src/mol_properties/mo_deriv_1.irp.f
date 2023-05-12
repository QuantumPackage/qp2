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
      size(ao_deriv_1_x,1),                                           &
      mo_deriv_1_x,                                                   &
      size(mo_deriv_1_x,1)                                            &
      )
  call ao_to_mo(                                                     &
      ao_deriv_1_y,                                                   &
      size(ao_deriv_1_y,1),                                           &
      mo_deriv_1_y,                                                   &
      size(mo_deriv_1_y,1)                                            &
      )
  call ao_to_mo(                                                     &
      ao_deriv_1_z,                                                   &
      size(ao_deriv_1_z,1),                                           &
      mo_deriv_1_z,                                                   &
      size(mo_deriv_1_z,1)                                            &
      )

END_PROVIDER
