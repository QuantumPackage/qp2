
! ---

subroutine give_all_erf_kl_ao(integrals_ao,mu_in,C_center)
  implicit none
  BEGIN_DOC
  ! Subroutine that returns all integrals over $r$ of type
  ! $\frac{ \erf(\mu * | r - R_C | ) }{ | r - R_C | }$
  END_DOC
  double precision, intent(in)   :: mu_in,C_center(3)
  double precision, intent(out)  :: integrals_ao(ao_num,ao_num)
  integer                        :: i,j,l,k,m
  double precision, allocatable :: integrals_ao_cart(:,:)
  allocate(integrals_ao_cart(ao_cart_num,ao_cart_num))
  call give_all_erf_kl_ao_cart(integrals_ao_cart,mu_in,C_center)
  call ao_cart_to_ao_basis(integrals_ao_cart, ao_cart_num, integrals_ao, ao_num)
end

! ---

double precision function NAI_pol_mult_erf_ao(i_ao, j_ao, mu_in, C_center)

  BEGIN_DOC
  !
  ! Computes the following integral :
  ! $\int_{-\infty}^{infty} dr \chi_i(r) \chi_j(r) \frac{\erf(\mu |r - R_C|)}{|r - R_C|}$.
  !
  END_DOC

  implicit none
  integer,          intent(in)   :: i_ao, j_ao
  double precision, intent(in)   :: mu_in, C_center(3)

  double precision, allocatable  :: integrals_ao(:,:),integrals_ao_cart(:,:)
  allocate(integrals_ao(ao_num,ao_num),integrals_ao_cart(ao_cart_num,ao_cart_num))
  call give_all_erf_kl_ao_cart(integrals_ao_cart,mu_in,C_center)
  call ao_cart_to_ao_basis(integrals_ao_cart, ao_cart_num, integrals_ao, ao_num)
  NAI_pol_mult_erf_ao = integrals_ao(i_ao, j_ao)
end function NAI_pol_mult_erf_ao

! ---
subroutine all_NAI_pol_mult_erf_ao_with1s(beta, B_center, mu_in, C_center)

  BEGIN_DOC
  !
  ! Computes ALL the following integral :
  ! $\int_{-\infty}^{infty} dr \chi_i(r) \chi_j(r) e^{-\beta (r - B_center)^2} \frac{\erf(\mu |r - R_C|)}{|r - R_C|}$.
  !
  END_DOC

  implicit none
  double precision, intent(in)   :: beta, B_center(3)
  double precision, intent(in)   :: mu_in, C_center(3)
  double precision, intent(out)  :: integrals_ao(ao_num,ao_num)

  double precision :: NAI_pol_mult_erf_ao_cart_with1s
  integer :: i,j
  double precision, allocatable  :: integrals_ao_cart(:,:)
  allocate(integrals_ao_cart(ao_cart_num,ao_cart_num))
  do i = 1, ao_cart_num
   do j = 1, ao_cart_num
    integrals_ao_car(j,i) = NAI_pol_mult_erf_ao_cart_with1s(i_ao, j_ao, beta, B_center, mu_in, C_center)
   enddo
  enddo
  call ao_cart_to_ao_basis(integrals_ao_cart, ao_cart_num, integrals_ao, ao_num)
end

double precision function NAI_pol_mult_erf_ao_with1s(i_ao, j_ao, beta, B_center, mu_in, C_center)

  BEGIN_DOC
  !
  ! Computes the following integral :
  ! $\int_{-\infty}^{infty} dr \chi_i(r) \chi_j(r) e^{-\beta (r - B_center)^2} \frac{\erf(\mu |r - R_C|)}{|r - R_C|}$.
  !
  END_DOC

  implicit none
  integer,          intent(in)   :: i_ao, j_ao
  double precision, intent(in)   :: beta, B_center(3)
  double precision, intent(in)   :: mu_in, C_center(3)
  double precision :: NAI_pol_mult_erf_ao_cart_with1s
  integer :: i,j
  double precision, allocatable  :: integrals_ao(:,:),integrals_ao_cart(:,:)
  allocate(integrals_ao(ao_num,ao_num),integrals_ao_cart(ao_cart_num,ao_cart_num))
  do i = 1, ao_cart_num
   do j = 1, ao_cart_num
    integrals_ao_car(j,i) = NAI_pol_mult_erf_ao_cart_with1s(i_ao, j_ao, beta, B_center, mu_in, C_center)
   enddo
  enddo
  call ao_cart_to_ao_basis(integrals_ao_cart, ao_cart_num, integrals_ao, ao_num)
  NAI_pol_mult_erf_ao_with1s = integrals_ao(i_ao, j_ao)

end function NAI_pol_mult_erf_ao_with1s

