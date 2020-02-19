use bitmasks

 BEGIN_PROVIDER [ integer(bit_kind), psi_cas_complex, (N_int,2,psi_det_size) ]
&BEGIN_PROVIDER [ complex*16, psi_cas_coef_complex,  (psi_det_size,n_states) ]
&BEGIN_PROVIDER [ integer, idx_cas_complex, (psi_det_size) ]
&BEGIN_PROVIDER [ integer, N_det_cas_complex ]
  implicit none
  BEGIN_DOC
  ! |CAS| wave function, defined from the application of the |CAS| bitmask on the
  ! determinants. idx_cas gives the indice of the |CAS| determinant in psi_det.
  END_DOC
  integer                        :: i, k, l
  logical                        :: good
  n_det_cas_complex = 0
  do i=1,N_det
    do l = 1, N_states
     psi_cas_coef_complex(i,l) = (0.d0,0.d0)
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
      n_det_cas_complex = n_det_cas_complex+1
      do k=1,N_int
        psi_cas_complex(k,1,n_det_cas_complex) = psi_det(k,1,i)
        psi_cas_complex(k,2,n_det_cas_complex) = psi_det(k,2,i)
      enddo
      idx_cas(n_det_cas_complex) = i
      do k=1,N_states
        psi_cas_coef_complex(n_det_cas_complex,k) = psi_coef_complex(i,k)
      enddo
    endif
  enddo
  call write_int(6,n_det_cas_complex, 'Number of determinants in the CAS')

END_PROVIDER


 BEGIN_PROVIDER [ integer(bit_kind), psi_cas_sorted_bit_complex, (N_int,2,psi_det_size) ]
&BEGIN_PROVIDER [ complex*16, psi_cas_coef_sorted_bit_complex, (psi_det_size,N_states) ]
 implicit none
 BEGIN_DOC
 ! |CAS| determinants sorted to accelerate the search of a random determinant in the wave
 ! function.
 END_DOC
 call sort_dets_by_det_search_key_complex(n_det_cas_complex, psi_cas_complex, psi_cas_coef_complex, size(psi_cas_coef_complex,1), &
     psi_cas_sorted_bit_complex, psi_cas_coef_sorted_bit_complex, N_states)

END_PROVIDER



 BEGIN_PROVIDER [ integer(bit_kind), psi_non_cas_complex,  (N_int,2,psi_det_size) ]
&BEGIN_PROVIDER [ complex*16, psi_non_cas_coef,_complex (psi_det_size,n_states) ]
&BEGIN_PROVIDER [ integer, idx_non_cas_complex,  (psi_det_size) ]
&BEGIN_PROVIDER [ integer, N_det_non_cas_complex ]
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
   do j=1,N_det_cas_complex
     call get_excitation_degree(psi_cas_complex(1,1,j), psi_det(1,1,k), degree, N_int)
     if (degree == 0) then
       in_cas = .True.
       exit
     endif
   enddo
   if (.not.in_cas) then
     double precision :: hij
     i_non_cas += 1
     do j=1,N_int
       psi_non_cas_complex(j,1,i_non_cas) = psi_det(j,1,k)
       psi_non_cas_complex(j,2,i_non_cas) = psi_det(j,2,k)
     enddo
     do j=1,N_states
       psi_non_cas_coef_complex(i_non_cas,j) = psi_coef_complex(k,j)
     enddo
     idx_non_cas_complex(i_non_cas) = k
   endif
 enddo
 N_det_non_cas_complex = i_non_cas
END_PROVIDER

 BEGIN_PROVIDER [ integer(bit_kind), psi_non_cas_sorted_bit_complex, (N_int,2,psi_det_size) ]
&BEGIN_PROVIDER [ complex*16, psi_non_cas_coef_sorted_bit_complex, (psi_det_size,N_states) ]
 implicit none
 BEGIN_DOC
 ! |CAS| determinants sorted to accelerate the search of a random determinant in the wave
 ! function.
 END_DOC
 !TODO: should this be n_det_non_cas_complex?
 call sort_dets_by_det_search_key_complex(N_det_cas_complex, psi_non_cas_complex, psi_non_cas_coef_complex, size(psi_non_cas_coef_complex,1), &
     psi_non_cas_sorted_bit_complex, psi_non_cas_coef_sorted_bit_complex, N_states)

END_PROVIDER


BEGIN_PROVIDER [complex*16, H_matrix_cas_complex, (N_det_cas_complex,N_det_cas_complex)]
 implicit none
 integer :: i,j
 complex*16 :: hij
  do i = 1, N_det_cas_complex
   do j = 1, N_det_cas_complex
    call i_h_j_complex(psi_cas_complex(1,1,i),psi_cas_complex(1,1,j),N_int,hij)
    H_matrix_cas_complex(i,j) = hij
   enddo
  enddo
END_PROVIDER

 BEGIN_PROVIDER [complex*16, psi_coef_cas_diagonalized_complex, (N_det_cas_complex,N_states)]
&BEGIN_PROVIDER [double precision, psi_cas_energy_diagonalized_complex, (N_states)]
 implicit none
 integer :: i,j
  double precision, allocatable  :: eigenvalues(:)
  complex*16, allocatable  :: eigenvectors(:,:)
  allocate (eigenvectors(size(H_matrix_cas,1),N_det_cas))
  allocate (eigenvalues(N_det_cas))
  call lapack_diag_complex(eigenvalues,eigenvectors,                       &
      H_matrix_cas_complex,size(H_matrix_cas_complex,1),N_det_cas_complex)
  do i = 1, N_states
   psi_cas_energy_diagonalized_complex(i) = eigenvalues(i)
   do j = 1, N_det_cas_complex
    psi_coef_cas_diagonalized_complex(j,i) = eigenvectors(j,i)
   enddo
  enddo


 END_PROVIDER

