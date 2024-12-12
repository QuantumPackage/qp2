program fit_1s_basis
 implicit none
 provide lmax_too_big
 integer :: i,j
 print*,'////////////////////////////////////////////////////'
 print*,'////////////////////////////////////////////////////'
 print*,'Fitting the original basis set on uncontracted s only functions '
 print*,'WARNING :: works for now with only P functions at most !!'
 print*,'WARNING :: otherwise it will stop '
 print*,'Writting the results in the extra_nuclei and ao_extra_basis folders of EZFIO'
 print*,'New number of atomic functions : '
 print*,'n_func_tot = ',n_func_tot
 call ezfio_set_extra_nuclei_extra_nucl_num(new_nucl_num)
 call ezfio_set_extra_nuclei_extra_nucl_charge(new_nucl_charge_1s)
 call ezfio_set_extra_nuclei_extra_nucl_coord(new_nucl_coord_1s)
 call ezfio_set_extra_nuclei_extra_nucl_label(new_nucl_label_1s)
 call ezfio_set_ao_extra_basis_ao_extra_num(n_func_tot)
 call ezfio_set_ao_extra_basis_ao_extra_nucl(new_ao_nucl_1s)
 call ezfio_set_ao_extra_basis_ao_extra_prim_num(new_ao_prim_num_1s)
 call ezfio_set_ao_extra_basis_ao_extra_coef(new_ao_coef_1s)
 call ezfio_set_ao_extra_basis_ao_extra_expo(new_ao_expo_1s)
 call ezfio_set_ao_extra_basis_ao_extra_power(new_ao_power_1s)
end

