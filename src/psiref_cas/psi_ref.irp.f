use bitmasks

 BEGIN_PROVIDER [ integer(bit_kind), psi_ref, (N_int,2,psi_det_size) ]
&BEGIN_PROVIDER [ double precision, psi_ref_coef,  (psi_det_size,n_states) ]
&BEGIN_PROVIDER [ integer, idx_ref, (psi_det_size) ]
&BEGIN_PROVIDER [ integer, N_det_ref ]
  implicit none
  BEGIN_DOC
  ! CAS wave function, defined from the application of the CAS bitmask on the
  ! determinants. idx_cas gives the indice of the CAS determinant in psi_det.
  END_DOC
  integer :: i,j,k
  N_det_ref = N_det_cas
  do i=1,N_det_ref
    do k=1,N_int
      psi_ref(k,1,i) = psi_cas(k,1,i)
      psi_ref(k,2,i) = psi_cas(k,2,i)
    enddo
    idx_ref(i) = idx_cas(i)
  enddo
  do k=1,N_states
    do i=1,N_det_ref
      psi_ref_coef(i,k) = psi_cas_coef(i,k)
    enddo
  enddo

END_PROVIDER

BEGIN_PROVIDER [ double precision, psi_ref_coef_inv, (psi_det_size,n_states) ]
 implicit none
 BEGIN_DOC
 ! 1/psi_ref_coef
 END_DOC
 integer :: i, i_state
 do i_state=1,N_states
  do i=1,N_det_ref
    psi_ref_coef_inv(i,i_state) = 1.d0/psi_ref_coef(i,i_state)
  enddo
 enddo

END_PROVIDER


 BEGIN_PROVIDER [ integer(bit_kind), psi_ref_restart, (N_int,2,psi_det_size) ]
&BEGIN_PROVIDER [ double precision, psi_ref_coef_restart,  (psi_det_size,n_states) ]
  implicit none
  BEGIN_DOC
  ! Projection of the CAS wave function on the restart wave function.
  END_DOC
  integer :: i,j,k
  integer, save                  :: ifirst

  if(ifirst == 0)then
   ifirst = 1
   do i=1,N_det_ref
     do k=1,N_int
       psi_ref_restart(k,1,i) = psi_cas(k,1,i)
       psi_ref_restart(k,2,i) = psi_cas(k,2,i)
     enddo
   enddo
   do k=1,N_states
     do i=1,N_det_ref
       psi_ref_coef_restart(i,k) = psi_cas_coef(i,k)
     enddo
   enddo
  endif

END_PROVIDER

 BEGIN_PROVIDER [double precision, norm_psi_ref, (N_states)]
&BEGIN_PROVIDER [double precision, inv_norm_psi_ref, (N_states)]
  implicit none
  BEGIN_DOC
  ! norm_psi_ref(j) is the sum of psi_ref_coef(i,j)^2 over
  ! i = 1..N_det_ref, for each state j = 1..N_states, i.e. the
  ! squared norm of the reference coefficients of state j.
  ! inv_norm_psi_ref(j) is 1/sqrt(norm_psi_ref(j)).
  ! No in-tree caller; used by external plugins (author,
  ! 2026-10-31); supplies the normalization factor of the
  ! reference wave function coefficients.
  END_DOC
  integer :: i,j
  norm_psi_ref = 0.d0
  do j = 1, N_states
   do i = 1, N_det_ref
    norm_psi_ref(j) += psi_ref_coef(i,j) * psi_ref_coef(i,j)
   enddo
   inv_norm_psi_ref(j) = 1.d0/(dsqrt(norm_psi_Ref(j)))
   print *,  inv_norm_psi_ref(j)
  enddo

 END_PROVIDER

 BEGIN_PROVIDER [double precision, psi_ref_coef_interm_norm, (N_det_ref,N_states)]
  implicit none
  BEGIN_DOC
  ! Reference coefficients scaled per state j by inv_norm_psi_ref(j):
  ! psi_ref_coef_interm_norm(i,j) = inv_norm_psi_ref(j) *
  ! psi_ref_coef(i,j), for i = 1..N_det_ref, j = 1..N_states. With
  ! this scaling the coefficients of state j have squared sum 1.
  ! TODO(doc): no in-tree caller; the intended use of this
  ! intermediate-normalization provider is unknown to the author
  ! (author, 2026-10-31).
  END_DOC
  integer :: i,j
  do j = 1, N_states
   do i = 1, N_det_ref
    psi_ref_coef_interm_norm(i,j) = inv_norm_psi_ref(j) * psi_ref_coef(i,j)
   enddo
  enddo
 END_PROVIDER

 BEGIN_PROVIDER [double precision, psi_non_ref_coef_interm_norm, (N_det_non_ref,N_states)]
  implicit none
  BEGIN_DOC
  ! Non-reference coefficients scaled by the reference normalization
  ! factor: psi_non_ref_coef_interm_norm(i,j) is
  ! psi_non_ref_coef(i,j) * inv_norm_psi_ref(j), for i = 1..N_det_non_ref,
  ! j = 1..N_states. The scaling factor comes from the reference part.
  ! TODO(doc): no in-tree caller; the intended use of this
  ! intermediate-normalization provider is unknown to the author
  ! (author, 2026-10-31).
  END_DOC
  integer :: i,j
  do j = 1, N_states
   do i = 1, N_det_non_ref
    psi_non_ref_coef_interm_norm(i,j) = psi_non_ref_coef(i,j) * inv_norm_psi_ref(j)
   enddo
  enddo
 END_PROVIDER
