 BEGIN_PROVIDER [ double precision, dressing_column_h, (N_det,N_states) ]
&BEGIN_PROVIDER [ double precision, dressing_column_s, (N_det,N_states) ]
 implicit none
 BEGIN_DOC
 ! \Delta_{state-specific}. \Psi
 END_DOC

 integer :: i,ii,k,j, l
 double precision :: f, tmp
 double precision, external :: u_dot_v
 logical, external :: detEq

 dressing_column_h(:,:) = 0.d0
 dressing_column_s(:,:) = 0.d0

 l = dressed_column_idx(1)
 do j = 1, n_det
   dressing_column_h(j,1) = 0.5d0*dmc_delta_h(j)
   dressing_column_h(l,1) -= 0.5d0 * psi_coef(j,1) * dmc_delta_h(j) /psi_coef(l,1)
 enddo
END_PROVIDER



