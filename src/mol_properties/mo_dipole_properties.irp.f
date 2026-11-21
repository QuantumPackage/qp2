 BEGIN_PROVIDER [double precision, mo_prop_dipole_x , (mo_num,mo_num)]
&BEGIN_PROVIDER [double precision, mo_prop_dipole_y , (mo_num,mo_num)]
&BEGIN_PROVIDER [double precision, mo_prop_dipole_z , (mo_num,mo_num)]
 implicit none
 BEGIN_DOC
 ! mo_prop_dipole_x(i,j) = <MO_i | x | MO_j>, the MO dipole
 ! matrix element between MO i and MO j; y and z components in
 ! mo_prop_dipole_y and mo_prop_dipole_z. If dft_grid_mo_dipole
 ! is true, the values are those of mo_dipole_x_becke (Becke
 ! grid quadrature); otherwise those of mo_dipole_x (AO dipole
 ! integrals transformed to the MO basis with ao_to_mo).
 END_DOC
 if(dft_grid_mo_dipole)then
   mo_prop_dipole_x=mo_dipole_x_becke
   mo_prop_dipole_y=mo_dipole_y_becke
   mo_prop_dipole_z=mo_dipole_z_becke
 else
  mo_prop_dipole_x=mo_dipole_x 
  mo_prop_dipole_y=mo_dipole_y 
  mo_prop_dipole_z=mo_dipole_z 
 endif
END_PROVIDER 
