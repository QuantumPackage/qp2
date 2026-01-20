 BEGIN_PROVIDER [double precision, mo_prop_dipole_x , (mo_num,mo_num)]
&BEGIN_PROVIDER [double precision, mo_prop_dipole_y , (mo_num,mo_num)]
&BEGIN_PROVIDER [double precision, mo_prop_dipole_z , (mo_num,mo_num)]
 implicit none
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
