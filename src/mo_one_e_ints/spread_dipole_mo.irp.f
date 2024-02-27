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

 BEGIN_PROVIDER [double precision, mo_spread_centered_x, (mo_num, mo_num) ]
&BEGIN_PROVIDER [double precision, mo_spread_centered_y, (mo_num, mo_num) ]
&BEGIN_PROVIDER [double precision, mo_spread_centered_z, (mo_num, mo_num) ]
 BEGIN_DOC
 ! array of the integrals of MO_i * (x^2 - <MO_i|x|MO_j>^2) MO_j = MO_i x^2 MO_j - (MO_i x MO_j)^2
 ! array of the integrals of MO_i * (y^2 - <MO_i|y|MO_j>^2) MO_j = MO_i y^2 MO_j - (MO_i y MO_j)^2
 ! array of the integrals of MO_i * (z^2 - <MO_i|z|MO_j>^2) MO_j = MO_i z^2 MO_j - (MO_i z MO_j)^2
 END_DOC
 implicit none
 integer :: i,j
 do i = 1, mo_num
  do j = 1, mo_num
   mo_spread_centered_x(j,i) = mo_spread_x(j,i) - mo_dipole_x(j,i)**2
   mo_spread_centered_y(j,i) = mo_spread_y(j,i) - mo_dipole_y(j,i)**2
   mo_spread_centered_z(j,i) = mo_spread_z(j,i) - mo_dipole_z(j,i)**2
  enddo
 enddo
END_PROVIDER
