!!! TODO:: optimize when "ao_extra_only_1s" is True


double precision function v_extra_nucl_extra_ao(i_ao,j_ao)
 implicit none
  BEGIN_DOC
  !
  ! Computes the following integral :
  ! $\int_{-\infty}^{infty} dr \chi_i(r) \chi_j(r) v_ne^{extra}(r)$.
  !
  !
  ! where BOTH $\chi_i(r)$ AND $\chi_j(r)$ belongs to the EXTRA basis 
  !
  ! and v_ne^{extra}(r) is the Coulomb potential coming from the EXTRA nuclei
  END_DOC
 integer, intent(in) ::i_ao,j_ao
 double precision :: mu_in,charge,coord(3)
 double precision :: NAI_pol_mult_erf_ao_extra
 mu_in = 1.d10
 integer :: i
 v_extra_nucl_extra_ao = 0.d0
 do i = 1, extra_nucl_num
  charge = extra_nucl_charge(i)
  coord(1:3) = extra_nucl_coord_transp(1:3,i)
  v_extra_nucl_extra_ao -= charge * NAI_pol_mult_erf_ao_extra(i_ao, j_ao, mu_in, coord)
 enddo
end

double precision function v_extra_nucl_ao(i_ao,j_ao)
 implicit none
  BEGIN_DOC
  !
  ! Computes the following integral :
  ! $\int_{-\infty}^{infty} dr \chi_i(r) \chi_j(r) v_ne(r)$.
  !
  !
  ! where BOTH $\chi_i(r)$ AND $\chi_j(r)$ belongs to the REGULAR basis 
  !
  ! and v_ne(r) is the Coulomb potential coming from the EXTRA nuclei
  END_DOC
 integer, intent(in) ::i_ao,j_ao
 integer :: i
 double precision :: mu_in, coord(3),charge, integral
 double precision :: NAI_pol_mult_erf_ao
 mu_in = 1.d+10
 do i = 1, extra_nucl_num
  coord(1:3) = extra_nucl_coord_transp(1:3,i)
  charge = extra_nucl_charge(i)
  v_extra_nucl_ao += -NAI_pol_mult_erf_ao(i_ao, j_ao, mu_in, coord) * charge
 enddo
end


double precision function v_nucl_extra_ao(i_ao,j_ao)
 implicit none
  BEGIN_DOC
  !
  ! Computes the following integral :
  ! $\int_{-\infty}^{infty} dr \chi_i(r) \chi_j(r) v_ne(r)$.
  !
  !
  ! where BOTH $\chi_i(r)$ AND $\chi_j(r)$ belongs to the EXTRA basis 
  !
  ! and v_ne(r) is the Coulomb potential coming from the REGULAR nuclei
  END_DOC
 integer, intent(in) ::i_ao,j_ao
 double precision :: mu_in,charge,coord(3)
 double precision :: NAI_pol_mult_erf_ao_extra
 mu_in = 1.d10
 integer :: i
 v_nucl_extra_ao = 0.d0
 do i = 1, nucl_num
  charge = nucl_charge(i)
  coord(1:3) = nucl_coord_transp(1:3,i)
  v_nucl_extra_ao -= charge * NAI_pol_mult_erf_ao_extra(i_ao, j_ao, mu_in, coord)
 enddo
end


double precision function v_extra_nucl_mixed_ao(i_ao,j_ao)
 implicit none
  BEGIN_DOC
  !
  ! Computes the following integral :
  ! $\int_{-\infty}^{infty} dr \chi_i(r) \chi_j(r) v_ne^{extra}(r)$.
  !
  !
  ! where $\chi_i(r)$ belongs to the EXTRA basis and $\chi_j(r)$ to the REGULAR basis
  !
  ! and v_ne^{extra}(r) is the Coulomb potential coming from the EXTRA nuclei
  END_DOC
 integer, intent(in) ::i_ao,j_ao
 double precision :: mu_in,charge,coord(3)
 double precision :: NAI_pol_mult_erf_ao_extra_mixed
 mu_in = 1.d10
 integer :: i
 v_extra_nucl_mixed_ao = 0.d0
 do i = 1, extra_nucl_num
  charge = extra_nucl_charge(i)
  coord(1:3) = extra_nucl_coord_transp(1:3,i)
  v_extra_nucl_mixed_ao -= charge * NAI_pol_mult_erf_ao_extra_mixed(i_ao, j_ao, mu_in, coord)
 enddo
end

double precision function v_nucl_mixed_ao(i_ao,j_ao)
 implicit none
  BEGIN_DOC
  !
  ! Computes the following integral :
  ! $\int_{-\infty}^{infty} dr \chi_i(r) \chi_j(r) v_ne(r)$.
  !
  !
  ! where $\chi_i(r)$ belongs to the EXTRA basis and $\chi_j(r)$ to the REGULAR basis
  !
  ! and v_ne(r) is the Coulomb potential coming from the REGULAR nuclei
  END_DOC
 integer, intent(in) ::i_ao,j_ao
 double precision :: mu_in,charge,coord(3)
 double precision :: NAI_pol_mult_erf_ao_extra_mixed
 mu_in = 1.d10
 integer :: i
 v_nucl_mixed_ao = 0.d0
 do i = 1, nucl_num
  charge = nucl_charge(i)
  coord(1:3) = nucl_coord_transp(1:3,i)
  v_nucl_mixed_ao -= charge * NAI_pol_mult_erf_ao_extra_mixed(i_ao, j_ao, mu_in, coord)
 enddo
end

