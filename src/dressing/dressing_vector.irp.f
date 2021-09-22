 BEGIN_PROVIDER [ double precision, dressing_column_h, (N_det,N_states) ]
&BEGIN_PROVIDER [ double precision, dressing_column_s, (N_det,N_states) ]
 implicit none
 BEGIN_DOC
 ! \Delta_{state-specific}. \Psi
 ! Diagonal element is divided by 2 because Delta = D + D^t
 END_DOC

 integer :: i,ii,k,j, l
 double precision :: f, tmp
 double precision, external :: u_dot_v
 logical, external :: detEq

 dressing_column_h(:,:) = 0.d0
 dressing_column_s(:,:) = 0.d0

 do k=1,N_states
   l = dressed_column_idx(k)
   do j = 1, n_det
     dressing_column_h(j,k) = delta_ij(k,j,1)
     dressing_column_s(j,k) = delta_ij(k,j,2)
     dressing_column_h(l,k) -= 0.5d0 * psi_coef(j,k) * delta_ij(k,j,1) /psi_coef(l,k)
     dressing_column_s(l,k) -= 0.5d0 * psi_coef(j,k) * delta_ij(k,j,2) /psi_coef(l,k)
   enddo
 enddo
END_PROVIDER



