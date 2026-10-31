BEGIN_PROVIDER [double precision, ao_canonical_nucl_elec_ints, (mo_num,mo_num)]
 implicit none
 BEGIN_DOC
 ! Nuclear attraction (nucleus-electron) integrals in the canonical
 ! orthonormalized AO basis (ao_ortho_canonical_coef).
 ! Element (j,i) is the matrix element between the jth and the ith
 ! ao_ortho_canonical orbital, obtained by transforming
 ! ao_integrals_n_e (the nucleus-electron integrals in the AO basis,
 ! including the pseudopotential integrals) with the coefficient
 ! matrix.
 END_DOC
 integer :: i1,j1,i,j
 double precision :: c_i1,c_j1

 ao_canonical_nucl_elec_ints = 0.d0
 !$OMP PARALLEL DO DEFAULT(none) &
 !$OMP PRIVATE(i,j,i1,j1,c_j1,c_i1) &
 !$OMP SHARED(mo_num,ao_num,ao_ortho_canonical_coef, &
 !$OMP   ao_canonical_nucl_elec_ints, ao_integrals_n_e)
 do i = 1, mo_num
   do j = 1, mo_num
    do i1 = 1,ao_num
     c_i1 = ao_ortho_canonical_coef(i1,i)
     do j1 = 1,ao_num
       c_j1 = c_i1*ao_ortho_canonical_coef(j1,j)
       ao_canonical_nucl_elec_ints(j,i) = ao_canonical_nucl_elec_ints(j,i) + &
                                           c_j1 * ao_integrals_n_e(j1,i1)
     enddo
    enddo
   enddo
 enddo
 !$OMP END PARALLEL DO
END_PROVIDER

