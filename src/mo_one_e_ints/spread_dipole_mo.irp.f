 BEGIN_PROVIDER [double precision, mo_dipole_x , (mo_num,mo_num)]
&BEGIN_PROVIDER [double precision, mo_dipole_y , (mo_num,mo_num)]
&BEGIN_PROVIDER [double precision, mo_dipole_z , (mo_num,mo_num)]
 BEGIN_DOC
 ! array of the integrals of MO_i * x MO_j
 ! array of the integrals of MO_i * y MO_j
 ! array of the integrals of MO_i * z MO_j
 END_DOC
 implicit none

  call ao_to_mo(                                                     &
      ao_dipole_x,                                                   &
      size(ao_dipole_x,1),                                           &
      mo_dipole_x,                                                   &
      size(mo_dipole_x,1)                                            &
      )
  call ao_to_mo(                                                     &
      ao_dipole_y,                                                   &
      size(ao_dipole_y,1),                                           &
      mo_dipole_y,                                                   &
      size(mo_dipole_y,1)                                            &
      )
  call ao_to_mo(                                                     &
      ao_dipole_z,                                                   &
      size(ao_dipole_z,1),                                           &
      mo_dipole_z,                                                   &
      size(mo_dipole_z,1)                                            &
      )

END_PROVIDER

 BEGIN_PROVIDER [double precision, mo_spread_x , (mo_num,mo_num)]
&BEGIN_PROVIDER [double precision, mo_spread_y , (mo_num,mo_num)]
&BEGIN_PROVIDER [double precision, mo_spread_z , (mo_num,mo_num)]
 BEGIN_DOC
 ! array of the integrals of MO_i * x^2 MO_j
 ! array of the integrals of MO_i * y^2 MO_j
 ! array of the integrals of MO_i * z^2 MO_j
 END_DOC
 implicit none
  call ao_to_mo(                                                     &
      ao_spread_x,                                                   &
      size(ao_spread_x,1),                                           &
      mo_spread_x,                                                   &
      size(mo_spread_x,1)                                            &
      )
  call ao_to_mo(                                                     &
      ao_spread_y,                                                   &
      size(ao_spread_y,1),                                           &
      mo_spread_y,                                                   &
      size(mo_spread_y,1)                                            &
      )
  call ao_to_mo(                                                     &
      ao_spread_z,                                                   &
      size(ao_spread_z,1),                                           &
      mo_spread_z,                                                   &
      size(mo_spread_z,1)                                            &
      )
END_PROVIDER

