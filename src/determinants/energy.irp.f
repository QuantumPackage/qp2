BEGIN_PROVIDER [ double precision, diagonal_H_matrix_on_psi_det, (N_det) ]
 implicit none
 BEGIN_DOC
 ! Diagonal of the Hamiltonian ordered as psi_det
 END_DOC
 double precision, external     :: diag_h_mat_elem
 integer                        :: i

 do i=1,N_det
   diagonal_H_matrix_on_psi_det(i) = diag_h_mat_elem(psi_det(1,1,i),N_int)
 enddo

END_PROVIDER

BEGIN_PROVIDER [ double precision, barycentric_electronic_energy, (N_states) ]
 implicit none
 BEGIN_DOC
 ! $E_n = \sum_i {c_i^{(n)}}^2 H_{ii}$
 END_DOC
 integer :: istate,i

 barycentric_electronic_energy(:) = 0.d0

 do istate=1,N_states
   do i=1,N_det
     barycentric_electronic_energy(istate) += psi_coef(i,istate)*psi_coef(i,istate)*diagonal_H_matrix_on_psi_det(i)
   enddo
 enddo

END_PROVIDER

