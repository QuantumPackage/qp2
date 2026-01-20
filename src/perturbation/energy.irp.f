BEGIN_PROVIDER [ double precision, energy_by_projection, (N_states)]
 implicit none
 BEGIN_DOC
 ! E = <Psi|H|0> / <Psi|0>
 END_DOC

 integer :: i, k, iref
 integer, allocatable :: degree(:), idx(:)
 double precision :: hij, tmp
 allocate(degree(N_det),idx(0:N_det))
 call i_h_j(psi_det(1,1,1),psi_det(1,1,1),N_int,hij)

 do i=1,N_states
   energy_by_projection(i) = 0.d0
   iref = maxloc(abs(psi_coef(:,i)),1)
   call get_excitation_degree_vector(psi_det,psi_det(1,1,iref),degree,N_int, N_det, idx)

   !$OMP PARALLEL PRIVATE(k,hij,tmp) 
   tmp = 0.d0
   !$OMP DO 
   do k=1,idx(0)
     call i_h_j(psi_det(1,1,iref),psi_det(1,1,idx(k)),N_int,hij)
     tmp = tmp + psi_coef(idx(k),i)*hij
   enddo
   !$OMP END DO
   !$OMP ATOMIC
   energy_by_projection(i) = energy_by_projection(i) + tmp
   !$OMP END PARALLEL

   energy_by_projection(i) = energy_by_projection(i) / psi_coef(iref,i) + nuclear_repulsion
 enddo
END_PROVIDER

