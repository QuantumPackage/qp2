BEGIN_PROVIDER [ double precision, effective_ao_extra_dm, (ao_extra_num, ao_extra_num)]
 implicit none
 BEGIN_DOC
 ! effective density matrix : rho_pq x N_p x N_q x E_pq x (pi/gamma_pq)^3/2
 !
 ! where rho_pq is the usual density matrix 
 !
 ! N_p and N_q are the normalization factors associated with the two Gaussians 
 !
 ! E_pq is the prefactor resulting from the Gaussian product 
 !
 ! gamma_pq = gamma_p + gamm_q
 END_DOC
 integer :: i,j
 do i = 1, ao_extra_num
  do j = 1, ao_extra_num
   effective_ao_extra_dm(j,i) = ao_extra_one_e_dm(j,i,1) * ao_extra_coef_normalized(j,1) * ao_extra_coef_normalized(i,1) & 
       * inv_pi_gamma_pq_3_2_ao_extra(j,i) * E_pq_ao_extra(j,i)
  enddo
 enddo

END_PROVIDER 

 BEGIN_PROVIDER [ double precision, gamma_pq_ao_extra, (ao_extra_num,ao_extra_num)]
&BEGIN_PROVIDER [ double precision, inv_gamma_pq_ao_extra, (ao_extra_num,ao_extra_num)]
&BEGIN_PROVIDER [ double precision, inv_pi_gamma_pq_3_2_ao_extra, (ao_extra_num,ao_extra_num)]
&BEGIN_PROVIDER [ double precision, sqrt_gamma_pq_ao_extra, (ao_extra_num,ao_extra_num)]
 implicit none
 BEGIN_DOC
 ! gamma_pq_ao_extra = gamma_p + gamma_q  
 !
 ! inv_gamma_pq_ao_extra = 1/(gamma_p + gamma_q)
 !
 ! sqrt_gamma_pq_ao_extra = sqrt(gamma_p + gamma_q)
 !
 ! WARNING :: VALID ONLY IN THE CASE OF A PURELY 1S BASIS 
 END_DOC
 include 'constants.include.F'
 integer :: i,j
 do i = 1, ao_extra_num
  do j = 1, ao_extra_num
   gamma_pq_ao_extra(j,i) = ao_extra_expo(j,1) + ao_extra_expo(i,1)
   inv_gamma_pq_ao_extra(j,i) = 1.d0/gamma_pq_ao_extra(j,i)
   sqrt_gamma_pq_ao_extra(j,i) = dsqrt(gamma_pq_ao_extra(j,i))
   inv_pi_gamma_pq_3_2_ao_extra(j,i) = (pi * inv_gamma_pq_ao_extra(j,i))**(1.5d0)
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, E_pq_ao_extra, (ao_extra_num,ao_extra_num)]
 implicit none
 BEGIN_DOC
 ! E_pq_ao_extra = exp(-alpha_p alpha_q/gamma_pq (Q_p - Q_q)^2)
 END_DOC 
 integer :: i,j
 do i = 1, ao_extra_num
  do j = 1, ao_extra_num
   E_pq_ao_extra(j,i) = dexp(-ao_extra_center_1s_dist(j,i)**2 * ao_extra_expo(j,1)*ao_extra_expo(i,1)*inv_gamma_pq_ao_extra(j,i))
  enddo
 enddo

END_PROVIDER 

BEGIN_PROVIDER [ double precision, ao_extra_center_1s, (3,ao_extra_num)]
 implicit none
 BEGIN_DOC
 ! Original position of each extra AO 
 END_DOC
 integer :: i,i_nucl
 do i = 1, ao_extra_num
  i_nucl= ao_extra_nucl(i)
  ao_extra_center_1s(1:3,i) = extra_nucl_coord(i_nucl,1:3)
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, ao_extra_center_1s_dist, (ao_extra_num,ao_extra_num)]
 implicit none
 integer :: i,j
 do i = 1, ao_extra_num
  do j = 1, ao_extra_num
   ao_extra_center_1s_dist(j,i) = (ao_extra_center_1s(1,j) - ao_extra_center_1s(1,i))**2 &
                                + (ao_extra_center_1s(2,j) - ao_extra_center_1s(2,i))**2 &
                                + (ao_extra_center_1s(3,j) - ao_extra_center_1s(3,i))**2
    ao_extra_center_1s_dist(j,i)=dsqrt(ao_extra_center_1s_dist(j,i))
  enddo
 enddo
END_PROVIDER 
