BEGIN_PROVIDER [ double precision, pot_vne_A_extra_basis, (ao_extra_num,ao_extra_num)]
 implicit none
  BEGIN_DOC
  !
  ! Computes the following integral :
  ! $\sum_{R in the USUAL nuclei} -Z <chi_i|1/|r-R||chi_j>$ 
  !
  !
  ! where $\chi_i(r)$ AND $\chi_j(r)$ belongs to the EXTRA basis 
  END_DOC
  integer :: mu,nu
  double precision :: v_nucl_extra_ao
  pot_vne_A_extra_basis = 0.d0
  do mu = 1, ao_extra_num
   do nu = 1, ao_extra_num
    pot_vne_A_extra_basis(nu,mu)= v_nucl_extra_ao(mu,nu)
   enddo
  enddo

END_PROVIDER 


BEGIN_PROVIDER [ double precision, pot_vne_extra_basis, (ao_num,ao_num)]
 implicit none
  BEGIN_DOC
  !
  ! Computes the following integral :
  ! $\sum_{R in EXTRA nuclei} -Z <chi_i|1/|r-R||chi_j>$ 
  !
  !
  ! where $\chi_i(r)$ AND $\chi_j(r)$ belongs to the USUAL basis 
  END_DOC
  integer :: mu,nu,k_nucl
  double precision :: mu_in, R_nucl(3),charge_nucl, integral
  double precision :: NAI_pol_mult_erf_ao
  mu_in = 10.d0**10
  pot_vne_extra_basis = 0.d0
  do mu = 1, ao_num
   do nu = 1, ao_num
    do k_nucl = 1, extra_nucl_num
     R_nucl(1:3) = extra_nucl_coord_transp(1:3,k_nucl)
     charge_nucl = extra_nucl_charge(k_nucl)
     integral = NAI_pol_mult_erf_ao(mu, nu, mu_in, R_nucl)
     pot_vne_extra_basis(nu,mu) += -integral * charge_nucl
    enddo
   enddo
  enddo

END_PROVIDER 



double precision function NAI_pol_mult_erf_ao_extra(i_ao, j_ao, mu_in, C_center)

  BEGIN_DOC
  !
  ! Computes the following integral :
  ! $\int_{-\infty}^{infty} dr \chi_i(r) \chi_j(r) \frac{\erf(\mu |r - R_C|)}{|r - R_C|}$.
  !
  !
  ! where $\chi_i(r)$ AND $\chi_j(r)$ belongs to the extra basis 
  END_DOC

  implicit none
  integer,          intent(in)   :: i_ao, j_ao
  double precision, intent(in)   :: mu_in, C_center(3)

  integer                        :: i, j, num_A, num_B, power_A(3), power_B(3), n_pt_in
  double precision               :: A_center(3), B_center(3), integral, alpha, beta

  double precision               :: NAI_pol_mult_erf

  num_A         = ao_extra_nucl(i_ao)
  power_A(1:3)  = ao_extra_power(i_ao,1:3)
  A_center(1:3) = extra_nucl_coord(num_A,1:3)
  num_B         = ao_extra_nucl(j_ao)
  power_B(1:3)  = ao_extra_power(j_ao,1:3)
  B_center(1:3) = extra_nucl_coord(num_B,1:3)

  n_pt_in = n_pt_max_extra_basis_integrals

  NAI_pol_mult_erf_ao_extra = 0.d0
  do i = 1, ao_extra_prim_num(i_ao)
    alpha = ao_extra_expo_ordered_transp(i,i_ao)
    do j = 1, ao_extra_prim_num(j_ao)
      beta = ao_extra_expo_ordered_transp(j,j_ao)

      integral = NAI_pol_mult_erf(A_center, B_center, power_A, power_B, alpha, beta, C_center, n_pt_in,mu_in)

      NAI_pol_mult_erf_ao_extra += integral * ao_extra_coef_normalized_ordered_transp(j,j_ao) * ao_extra_coef_normalized_ordered_transp(i,i_ao)
    enddo
  enddo

end function NAI_pol_mult_erf_ao_extra

! ---

double precision function NAI_pol_mult_erf_ao_extra_mixed(i_ao, j_ao, mu_in, C_center)

  BEGIN_DOC
  !
  ! Computes the following integral :
  ! $\int_{-\infty}^{infty} dr \chi_i(r) \chi_j(r) \frac{\erf(\mu |r - R_C|)}{|r - R_C|}$.
  !
  !
  ! where $\chi_i(r)$ belongs to the extra basis and $\chi_j(r)$ to the regular basis 
  END_DOC

  implicit none
  integer,          intent(in)   :: i_ao, j_ao
  double precision, intent(in)   :: mu_in, C_center(3)

  integer                        :: i, j, num_A, num_B, power_A(3), power_B(3), n_pt_in
  double precision               :: A_center(3), B_center(3), integral, alpha, beta

  double precision               :: NAI_pol_mult_erf

  ! A = chi_i == extra basis 
  num_A         = ao_extra_nucl(i_ao)
  power_A(1:3)  = ao_extra_power(i_ao,1:3)
  A_center(1:3) = extra_nucl_coord(num_A,1:3)
  ! B = chi_j == regular basis 
  num_B         = ao_nucl(j_ao)
  power_B(1:3)  = ao_power(j_ao,1:3)
  B_center(1:3) = nucl_coord(num_B,1:3)

  n_pt_in = max(n_pt_max_integrals,n_pt_max_extra_basis_integrals)

  NAI_pol_mult_erf_ao_extra_mixed = 0.d0
  do i = 1, ao_extra_prim_num(i_ao)
    alpha = ao_extra_expo_ordered_transp(i,i_ao)
    do j = 1, ao_prim_num(j_ao)
      beta = ao_expo_ordered_transp(j,j_ao)

      integral = NAI_pol_mult_erf(A_center, B_center, power_A, power_B, alpha, beta, C_center, n_pt_in,mu_in)

      NAI_pol_mult_erf_ao_extra_mixed += integral * ao_coef_normalized_ordered_transp(j,j_ao) * ao_extra_coef_normalized_ordered_transp(i,i_ao)
    enddo
  enddo

end 

! ---

