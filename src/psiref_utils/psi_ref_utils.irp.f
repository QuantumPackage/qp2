use bitmasks


 BEGIN_PROVIDER [ integer(bit_kind), psi_ref_sorted_bit, (N_int,2,psi_det_size) ]
&BEGIN_PROVIDER [ double precision, psi_ref_coef_sorted_bit, (psi_det_size,N_states) ]
 implicit none
 BEGIN_DOC
 ! Reference determinants sorted to accelerate the search of a random determinant in the wave
 ! function.
 END_DOC
 call sort_dets_by_det_search_key(N_det_ref, psi_ref, psi_ref_coef, &
     psi_ref_sorted_bit, psi_ref_coef_sorted_bit, N_states)

END_PROVIDER


BEGIN_PROVIDER [ double precision, psi_ref_coef_transp, (n_states,psi_det_size) ]
 implicit none
 BEGIN_DOC
! Transposed psi_ref_coef
 END_DOC
 integer :: i,j
 do j=1,N_det_ref
   do i=1, n_states
     psi_ref_coef_transp(i,j) = psi_ref_coef(j,i)
   enddo
 enddo
END_PROVIDER

BEGIN_PROVIDER [ double precision, psi_ref_coef_normalized,  (psi_det_size,n_states) ]
 implicit none
 BEGIN_DOC
! Normalized coefficients of the reference
 END_DOC
 integer :: i,j,k
 do k=1,N_states
   do j=1,N_det_ref
     psi_ref_coef_normalized(j,k) = psi_ref_coef(j,k)
   enddo
   call normalize(psi_ref_coef_normalized(1,k), N_det_ref)
 enddo

END_PROVIDER


BEGIN_PROVIDER [ double precision, psi_non_ref_coef_transp, (n_states,psi_det_size) ]
 implicit none
 BEGIN_DOC
! Transposed psi_non_ref_coef
 END_DOC
 integer :: i,j
 do j=1,N_det_non_ref
   do i=1, n_states
     psi_non_ref_coef_transp(i,j) = psi_non_ref_coef(j,i)
   enddo
 enddo
END_PROVIDER

 BEGIN_PROVIDER [ integer(bit_kind), psi_non_ref,  (N_int,2,psi_det_size) ]
&BEGIN_PROVIDER [ double precision, psi_non_ref_coef, (psi_det_size,n_states) ]
&BEGIN_PROVIDER [ integer, idx_non_ref,  (psi_det_size) ]
&BEGIN_PROVIDER [ integer, idx_non_ref_rev,  (psi_det_size) ]
&BEGIN_PROVIDER [ integer, N_det_non_ref ]
 implicit none
 BEGIN_DOC
  ! Set of determinants which are not part of the reference, defined from the application
  ! of the reference bitmask on the determinants.
  ! idx_non_ref gives the indice of the determinant in psi_det.
  ! idx_non_ref_rev gives the reverse.
 END_DOC
 integer                        :: i_non_ref,j,k
 integer                        :: degree
 logical                        :: in_ref
 i_non_ref =0
 idx_non_ref_rev = 0
 do k=1,N_det
   in_ref = .False.
   do j=1,N_det_ref
     call get_excitation_degree(psi_ref(1,1,j), psi_det(1,1,k), degree, N_int)
     if (degree == 0) then
       in_ref = .True.
       exit
     endif
   enddo
   if (.not.in_ref) then
     double precision :: hij
     i_non_ref += 1
     do j=1,N_int
       psi_non_ref(j,1,i_non_ref) = psi_det(j,1,k)
       psi_non_ref(j,2,i_non_ref) = psi_det(j,2,k)
     enddo
     do j=1,N_states
       psi_non_ref_coef(i_non_ref,j) = psi_coef(k,j)
     enddo
     idx_non_ref(i_non_ref) = k
     idx_non_ref_rev(k) = i_non_ref
   endif
 enddo
 N_det_non_ref = i_non_ref
 if (N_det_non_ref < 1) then
   print *,  'Warning : All determinants are in the reference'
 endif
END_PROVIDER

 BEGIN_PROVIDER [ integer(bit_kind), psi_non_ref_restart,  (N_int,2,psi_det_size) ]
&BEGIN_PROVIDER [ double precision, psi_non_ref_coef_restart, (psi_det_size,n_states) ]
 implicit none
 BEGIN_DOC
  ! Set of determinants which are not part of the reference, defined from the application
  ! of the reference bitmask on the determinants.
  ! idx_non_ref gives the indice of the determinant in psi_det.
  ! But this is with respect to the restart wave function.
 END_DOC
 integer                        :: i_non_ref,j,k
 integer                        :: degree
 logical                        :: in_ref
 integer, save                  :: ifirst = 0
 if(ifirst==0)then
  ifirst = 1
  i_non_ref =0
  do k=1,N_det
    in_ref = .False.
    do j=1,N_det_ref
      call get_excitation_degree(psi_ref(1,1,j), psi_det(1,1,k), degree, N_int)
      if (degree == 0) then
        in_ref = .True.
        exit
      endif
    enddo
    if (.not.in_ref) then
      double precision :: hij
      i_non_ref += 1
      do j=1,N_int
        psi_non_ref_restart(j,1,i_non_ref) = psi_det(j,1,k)
        psi_non_ref_restart(j,2,i_non_ref) = psi_det(j,2,k)
      enddo
      do j=1,N_states
        psi_non_ref_coef_restart(i_non_ref,j) = psi_coef(k,j)
      enddo
    endif
  enddo
 endif
END_PROVIDER



 BEGIN_PROVIDER [ integer(bit_kind), psi_non_ref_sorted_bit, (N_int,2,psi_det_size) ]
&BEGIN_PROVIDER [ double precision, psi_non_ref_coef_sorted_bit, (psi_det_size,N_states) ]
 implicit none
 BEGIN_DOC
 ! Reference determinants sorted to accelerate the search of a random determinant in the wave
 ! function.
 END_DOC
 call sort_dets_by_det_search_key(N_det_ref, psi_non_ref, psi_non_ref_coef, size(psi_non_ref_coef,1), &
     psi_non_ref_sorted_bit, psi_non_ref_coef_sorted_bit, N_states)

END_PROVIDER


BEGIN_PROVIDER [double precision, H_matrix_ref, (N_det_ref,N_det_ref)]
 implicit none
 integer :: i,j
 double precision :: hij
  do i = 1, N_det_ref
   do j = 1, N_det_ref
    call i_H_j(psi_ref(1,1,i),psi_ref(1,1,j),N_int,hij)
    H_matrix_ref(i,j) = hij
   enddo
  enddo
END_PROVIDER

 BEGIN_PROVIDER [double precision, psi_ref_coef_diagonalized, (N_det_ref,N_states)]
&BEGIN_PROVIDER [double precision, psi_ref_energy_diagonalized, (N_states)]
 implicit none
 integer :: i,j
  double precision, allocatable  :: eigenvectors(:,:), eigenvalues(:)
  allocate (eigenvectors(size(H_matrix_ref,1),N_det_ref))
  allocate (eigenvalues(N_det_ref))
  call lapack_diag(eigenvalues,eigenvectors,                       &
      H_matrix_ref,size(H_matrix_ref,1),N_det_ref)
  do i = 1, N_states
   psi_ref_energy_diagonalized(i) = eigenvalues(i)
   do j = 1, N_det_ref
    psi_ref_coef_diagonalized(j,i) = eigenvectors(j,i)
   enddo
  enddo
  deallocate (eigenvectors)
  deallocate (eigenvalues)


 END_PROVIDER

 BEGIN_PROVIDER [double precision, psi_ref_energy, (N_states)]
 implicit none
 integer :: i,j,k
 double precision :: hij,norm,u_dot_v
  psi_ref_energy = 0.d0


  do k = 1, N_states
   norm = 0.d0
   do i = 1, N_det_ref
    norm += psi_ref_coef(i,k) * psi_ref_coef(i,k)
    do j = 1, N_det_ref
      psi_ref_energy(k) += psi_ref_coef(i,k) * psi_ref_coef(j,k) * H_matrix_ref(i,j)
    enddo
   enddo
   psi_ref_energy(k) = psi_ref_energy(k) /norm
  enddo

END_PROVIDER


logical function is_in_psi_ref(key,Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
! True if the determinant ``det`` is in the wave function
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: key(Nint,2)
  integer, external              :: get_index_in_psi_ref_sorted_bit

  !DIR$ FORCEINLINE
  is_in_psi_ref = get_index_in_psi_ref_sorted_bit(key,Nint) > 0
end

integer function get_index_in_psi_ref_sorted_bit(key,Nint)
  use bitmasks
  BEGIN_DOC
! Returns the index of the determinant in the ``psi_ref_sorted_bit`` array
  END_DOC
  implicit none

  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: key(Nint,2)

  integer                        :: i, ibegin, iend, istep, l
  integer*8                      :: det_ref, det_search
  integer*8, external            :: det_search_key
  logical                        :: in_wavefunction

  in_wavefunction = .False.
  get_index_in_psi_ref_sorted_bit = 0
  ibegin = 1
  iend   = N_det+1

  !DIR$ FORCEINLINE
  det_ref = det_search_key(key,Nint)
  !DIR$ FORCEINLINE
  det_search = det_search_key(psi_ref_sorted_bit(1,1,1),Nint)

  istep = shiftr(iend-ibegin,1)
  i=ibegin+istep
  do while (istep > 0)
    !DIR$ FORCEINLINE
    det_search = det_search_key(psi_ref_sorted_bit(1,1,i),Nint)
    if ( det_search > det_ref ) then
      iend = i
    else if ( det_search == det_ref ) then
      exit
    else
      ibegin = i
    endif
    istep = shiftr(iend-ibegin,1)
    i = ibegin + istep
  end do

  !DIR$ FORCEINLINE
  do while (det_search_key(psi_ref_sorted_bit(1,1,i),Nint) == det_ref)
    i = i-1
    if (i == 0) then
      exit
    endif
  enddo
  i += 1

  if (i > N_det) then
    return
  endif

  !DIR$ FORCEINLINE
  do while (det_search_key(psi_ref_sorted_bit(1,1,i),Nint) == det_ref)
    if ( (key(1,1) /= psi_ref_sorted_bit(1,1,i)).or.                               &
          (key(1,2) /= psi_ref_sorted_bit(1,2,i)) ) then
      continue
    else
      in_wavefunction = .True.
      !DIR$ IVDEP
      !DIR$ LOOP COUNT MIN(3)
      do l=2,Nint
        if ( (key(l,1) /= psi_ref_sorted_bit(l,1,i)).or.                           &
              (key(l,2) /= psi_ref_sorted_bit(l,2,i)) ) then
          in_wavefunction = .False.
        endif
      enddo
      if (in_wavefunction) then
        get_index_in_psi_ref_sorted_bit = i
!        exit
        return
      endif
    endif
    i += 1
    if (i > N_det) then
!      exit
      return
    endif

  enddo

end

BEGIN_PROVIDER [double precision, ref_hamiltonian_matrix, (n_det_ref,n_det_ref)]
 BEGIN_DOC
 ! H matrix in the Reference space
 END_DOC
 implicit none
 integer :: i,j
 double precision :: hij
 do i = 1, N_det_ref
  do j = 1, N_det_ref
   call i_H_j(psi_ref(1,1,i),psi_ref(1,1,j),N_int,hij)
   ref_hamiltonian_matrix(i,j) = hij
  enddo
 enddo
END_PROVIDER


BEGIN_PROVIDER [ integer, idx_non_ref_from_sorted, (N_det) ]
  implicit none
  integer :: i,inpsisor

  idx_non_ref_from_sorted = 0

  do i=1,N_det
    inpsisor = psi_det_sorted_order(i)
    if(inpsisor <= 0) stop "idx_non_ref_from_sorted"
    idx_non_ref_from_sorted(inpsisor) = idx_non_ref_rev(i)
  end do
END_PROVIDER

