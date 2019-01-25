BEGIN_PROVIDER [double precision, ao_ortho_lowdin_nucl_elec_integrals, (mo_num,mo_num)]
 implicit none
 integer :: i1,j1,i,j
 double precision :: c_i1,c_j1

 ao_ortho_lowdin_nucl_elec_integrals = 0.d0
 !$OMP PARALLEL DO DEFAULT(none) &
 !$OMP PRIVATE(i,j,i1,j1,c_j1,c_i1) &
 !$OMP SHARED(mo_num,ao_num,ao_ortho_lowdin_coef, &
 !$OMP   ao_ortho_lowdin_nucl_elec_integrals, ao_integrals_n_e)
 do i = 1, mo_num
   do j = 1, mo_num
    do i1 = 1,ao_num
     c_i1 = ao_ortho_lowdin_coef(i1,i)
     do j1 = 1,ao_num
       c_j1 = c_i1*ao_ortho_lowdin_coef(j1,j)
       ao_ortho_lowdin_nucl_elec_integrals(j,i) = ao_ortho_lowdin_nucl_elec_integrals(j,i) + &
                                           c_j1 * ao_integrals_n_e(j1,i1)
     enddo
    enddo
   enddo
 enddo
 !$OMP END PARALLEL DO
END_PROVIDER

