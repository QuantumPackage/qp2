BEGIN_PROVIDER [ double precision, mo_r_coef_normalized, (ao_num,mo_num) ]
 implicit none
 integer :: i,j,k
 double precision :: norm
 do i = 1, mo_num
  norm = 0.d0
  do j = 1, ao_num
   do k = 1, ao_num
    norm += mo_r_coef(k,i) * mo_r_coef(j,i) * ao_overlap(k,j)
   enddo
  enddo
  norm = 1.d0/dsqrt(norm)
  do j = 1, ao_num
   mo_r_coef_normalized(j,i) = mo_r_coef(j,i) * norm
  enddo
 enddo
END_PROVIDER

BEGIN_PROVIDER [ double precision, tc_spin_dens_right_only, (ao_num, ao_num)]
 implicit none
 integer :: i,j,k
 tc_spin_dens_right_only = 0.d0
 do i = elec_beta_num+1, elec_alpha_num
  do j = 1, ao_num
   do k = 1, ao_num
    tc_spin_dens_right_only(k,j) += mo_r_coef_normalized(k,i) * mo_r_coef_normalized(j,i)
   enddo
  enddo
 enddo
END_PROVIDER
