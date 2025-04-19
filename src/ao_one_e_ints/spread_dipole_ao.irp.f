 BEGIN_PROVIDER [ double precision, ao_spread_x, (ao_num,ao_num)]
 implicit none
 BEGIN_DOC
 ! * array of the integrals of AO_i * x^2 AO_j
 END_DOC
   call ao_cart_to_ao_basis(ao_cart_spread_x, ao_cart_num, ao_spread_x, ao_num)
 END_PROVIDER
 BEGIN_PROVIDER [ double precision, ao_spread_y, (ao_num,ao_num)]
 implicit none
 BEGIN_DOC
 ! * array of the integrals of AO_i * y^2 AO_j
 END_DOC
   call ao_cart_to_ao_basis(ao_cart_spread_y, ao_cart_num, ao_spread_y, ao_num)
 END_PROVIDER
 BEGIN_PROVIDER [ double precision, ao_spread_z, (ao_num,ao_num)]
 implicit none
 BEGIN_DOC
 ! * array of the integrals of AO_i * z^2 AO_j
 END_DOC
   call ao_cart_to_ao_basis(ao_cart_spread_z, ao_cart_num, ao_spread_z, ao_num)
 END_PROVIDER


 BEGIN_PROVIDER [ double precision, ao_dipole_x, (ao_num,ao_num)]
 BEGIN_DOC
 ! * array of the integrals of AO_i * x AO_j
 END_DOC
 implicit none
   call ao_cart_to_ao_basis(ao_cart_dipole_x, ao_cart_num, ao_dipole_x, ao_num)
 END_PROVIDER
 BEGIN_PROVIDER [ double precision, ao_dipole_y, (ao_num,ao_num)]
 BEGIN_DOC
 ! * array of the integrals of AO_i * y AO_j
 END_DOC
 implicit none
   call ao_cart_to_ao_basis(ao_cart_dipole_y, ao_cart_num, ao_dipole_y, ao_num)
 END_PROVIDER
 BEGIN_PROVIDER [ double precision, ao_dipole_z, (ao_num,ao_num)]
 BEGIN_DOC
 ! * array of the integrals of AO_i * z AO_j
 END_DOC
 implicit none
   call ao_cart_to_ao_basis(ao_cart_dipole_z, ao_cart_num, ao_dipole_z, ao_num)
 END_PROVIDER


  BEGIN_PROVIDER [ double precision, ao_deriv_1_x, (ao_num,ao_num)]
 BEGIN_DOC
 ! * array of the integrals of AO_i * d/dx  AO_j
 END_DOC
 implicit none
   call ao_cart_to_ao_basis(ao_cart_deriv_1_x, ao_cart_num, ao_deriv_1_x, ao_num)
 END_PROVIDER
  BEGIN_PROVIDER [ double precision, ao_deriv_1_y, (ao_num,ao_num)]
 BEGIN_DOC
 ! * array of the integrals of AO_i * d/dy  AO_j
 END_DOC
 implicit none
   call ao_cart_to_ao_basis(ao_cart_deriv_1_y, ao_cart_num, ao_deriv_1_y, ao_num)
 END_PROVIDER
  BEGIN_PROVIDER [ double precision, ao_deriv_1_z, (ao_num,ao_num)]
 BEGIN_DOC
 ! * array of the integrals of AO_i * d/dz  AO_j
 END_DOC
 implicit none
   call ao_cart_to_ao_basis(ao_cart_deriv_1_z, ao_cart_num, ao_deriv_1_z, ao_num)
 END_PROVIDER

