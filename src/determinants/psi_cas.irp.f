use bitmasks

 BEGIN_PROVIDER [ integer(bit_kind), psi_cas, (N_int,2,psi_det_size) ]
&BEGIN_PROVIDER [ double precision, psi_cas_coef,  (psi_det_size,n_states) ]
&BEGIN_PROVIDER [ integer, idx_cas, (psi_det_size) ]
&BEGIN_PROVIDER [ integer, N_det_cas ]
  implicit none
  BEGIN_DOC
  ! |CAS| wave function, defined from the application of the |CAS| bitmask on the
  ! determinants. idx_cas gives the indice of the |CAS| determinant in psi_det.
  END_DOC
  integer                        :: i, k, l
  logical                        :: good
  N_det_cas = 0
  do i=1,N_det
    do l = 1, N_states
     psi_cas_coef(i,l) = 0.d0
    enddo
    good = .True.
    do k=1,N_int
      good = good .and. (                                          &
          iand(not(act_bitmask(k,1)), psi_det(k,1,i)) ==         &
          iand(not(act_bitmask(k,1)), hf_bitmask(k,1)) ) .and. (  &
          iand(not(act_bitmask(k,2)), psi_det(k,2,i)) ==         &
          iand(not(act_bitmask(k,2)), hf_bitmask(k,2)) )
    enddo
    if (good) then
      exit
    endif
    if (good) then
      N_det_cas = N_det_cas+1
      do k=1,N_int
        psi_cas(k,1,N_det_cas) = psi_det(k,1,i)
        psi_cas(k,2,N_det_cas) = psi_det(k,2,i)
      enddo
      idx_cas(N_det_cas) = i
      do k=1,N_states
        psi_cas_coef(N_det_cas,k) = psi_coef(i,k)
      enddo
    endif
  enddo
  call write_int(6,N_det_cas, 'Number of determinants in the CAS')

END_PROVIDER


 BEGIN_PROVIDER [ integer(bit_kind), psi_cas_sorted_bit, (N_int,2,psi_det_size) ]
&BEGIN_PROVIDER [ double precision, psi_cas_coef_sorted_bit, (psi_det_size,N_states) ]
 implicit none
 BEGIN_DOC
 ! |CAS| determinants sorted to accelerate the search of a random determinant in the wave
 ! function.
 END_DOC
 call sort_dets_by_det_search_key(N_det_cas, psi_cas, psi_cas_coef, size(psi_cas_coef,1), &
     psi_cas_sorted_bit, psi_cas_coef_sorted_bit, N_states)

END_PROVIDER



 BEGIN_PROVIDER [ integer(bit_kind), psi_non_cas,  (N_int,2,psi_det_size) ]
&BEGIN_PROVIDER [ double precision, psi_non_cas_coef, (psi_det_size,n_states) ]
&BEGIN_PROVIDER [ integer, idx_non_cas,  (psi_det_size) ]
&BEGIN_PROVIDER [ integer, N_det_non_cas ]
 implicit none
 BEGIN_DOC
  ! Set of determinants which are not part of the |CAS|, defined from the application
  ! of the |CAS| bitmask on the determinants.
  ! idx_non_cas gives the indice of the determinant in psi_det.
 END_DOC
 integer                        :: i_non_cas,j,k
 integer                        :: degree
 logical                        :: in_cas
 i_non_cas =0
 do k=1,N_det
   in_cas = .False.
   do j=1,N_det_cas
     call get_excitation_degree(psi_cas(1,1,j), psi_det(1,1,k), degree, N_int)
     if (degree == 0) then
       in_cas = .True.
       exit
     endif
   enddo
   if (.not.in_cas) then
     double precision :: hij
     i_non_cas += 1
     do j=1,N_int
       psi_non_cas(j,1,i_non_cas) = psi_det(j,1,k)
       psi_non_cas(j,2,i_non_cas) = psi_det(j,2,k)
     enddo
     do j=1,N_states
       psi_non_cas_coef(i_non_cas,j) = psi_coef(k,j)
     enddo
     idx_non_cas(i_non_cas) = k
   endif
 enddo
 N_det_non_cas = i_non_cas
END_PROVIDER

 BEGIN_PROVIDER [ integer(bit_kind), psi_non_cas_sorted_bit, (N_int,2,psi_det_size) ]
&BEGIN_PROVIDER [ double precision, psi_non_cas_coef_sorted_bit, (psi_det_size,N_states) ]
 implicit none
 BEGIN_DOC
 ! |CAS| determinants sorted to accelerate the search of a random determinant in the wave
 ! function.
 END_DOC
 call sort_dets_by_det_search_key(N_det_cas, psi_non_cas, psi_non_cas_coef, size(psi_non_cas_coef,1), &
     psi_non_cas_sorted_bit, psi_non_cas_coef_sorted_bit, N_states)

END_PROVIDER


BEGIN_PROVIDER [double precision, H_matrix_cas, (N_det_cas,N_det_cas)]
 implicit none
 integer :: i,j
 double precision :: hij
  do i = 1, N_det_cas
   do j = 1, N_det_cas
    call i_H_j(psi_cas(1,1,i),psi_cas(1,1,j),N_int,hij)
    H_matrix_cas(i,j) = hij
   enddo
  enddo
END_PROVIDER

 BEGIN_PROVIDER [double precision, psi_coef_cas_diagonalized, (N_det_cas,N_states)]
&BEGIN_PROVIDER [double precision, psi_cas_energy_diagonalized, (N_states)]
 implicit none
 integer :: i,j
  double precision, allocatable  :: eigenvectors(:,:), eigenvalues(:)
  allocate (eigenvectors(size(H_matrix_cas,1),N_det_cas))
  allocate (eigenvalues(N_det_cas))
  call lapack_diag(eigenvalues,eigenvectors,                       &
      H_matrix_cas,size(H_matrix_cas,1),N_det_cas)
  do i = 1, N_states
   psi_cas_energy_diagonalized(i) = eigenvalues(i)
   do j = 1, N_det_cas
    psi_coef_cas_diagonalized(j,i) = eigenvectors(j,i)
   enddo
  enddo


 END_PROVIDER

 BEGIN_PROVIDER [double precision, psi_cas_energy, (N_states)]
 implicit none
 BEGIN_DOC
! Variational energy of $\Psi_{CAS}$, where $\Psi_{CAS} =  \sum_{I \in CAS} \I \rangle \langle I | \Psi \rangle$.
 END_DOC
 integer :: i,j,k
 double precision :: hij,norm,u_dot_v
  psi_cas_energy = 0.d0


  do k = 1, N_states
   norm = 0.d0
   do i = 1, N_det_cas
    norm += psi_cas_coef(i,k) * psi_cas_coef(i,k)
    do j = 1, N_det_cas
      psi_cas_energy(k) += psi_cas_coef(i,k) * psi_cas_coef(j,k) * H_matrix_cas(i,j)
    enddo
   enddo
   psi_cas_energy(k) = psi_cas_energy(k) /norm
  enddo

END_PROVIDER




