
! ---

 BEGIN_PROVIDER [ double precision, overlap_states,     (N_states,N_states) ]
&BEGIN_PROVIDER [ double precision, overlap_states_inv, (N_states,N_states) ]

  BEGIN_DOC
  !
  ! S_kl = ck.T x cl 
  !      = psi_coef(:,k).T x psi_coef(:,l)
  !
  END_DOC

  implicit none
  integer          :: i
  double precision :: o_tmp

  if(N_states == 1) then

    o_tmp = 0.d0
    do i = 1, N_det
      o_tmp = o_tmp + psi_coef(i,1) * psi_coef(i,1)
    enddo
    overlap_states    (1,1) = o_tmp
    overlap_states_inv(1,1) = 1.d0 / o_tmp

  else

    call dgemm( 'T', 'N', N_states, N_states, N_det, 1.d0                &
              , psi_coef, size(psi_coef, 1), psi_coef, size(psi_coef, 1) &
              , 0.d0, overlap_states, size(overlap_states, 1)            )

    call get_inverse(overlap_states, N_states, N_states, overlap_states_inv, N_states)

  endif

END_PROVIDER

! ---

