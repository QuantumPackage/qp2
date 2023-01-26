use bitmasks

BEGIN_PROVIDER [ double precision, threshold_selectors ]
 implicit none
 BEGIN_DOC
 ! Thresholds on selectors (fraction of the square of the norm)
 END_DOC
 threshold_selectors = dsqrt(threshold_generators)
END_PROVIDER

BEGIN_PROVIDER [ integer, N_det_selectors]
  implicit none
  BEGIN_DOC
  ! For Single reference wave functions, the number of selectors is 1 : the
  ! Hartree-Fock determinant
  END_DOC
  integer                        :: i
  double precision               :: norm, norm_max
  call write_time(6)
  N_det_selectors = N_det
  norm = 1.d0
  do i=1,N_det
    norm = norm - psi_average_norm_contrib_tc(i)
    if (norm - 1.d-10 < 1.d0 - threshold_selectors) then
      N_det_selectors = i
      exit
    endif
  enddo
  N_det_selectors = max(N_det_selectors,N_det_generators)
  call write_int(6,N_det_selectors,'Number of selectors')
END_PROVIDER

 BEGIN_PROVIDER [ integer(bit_kind), psi_selectors, (N_int,2,psi_selectors_size) ]
&BEGIN_PROVIDER [ double precision, psi_selectors_coef, (psi_selectors_size,N_states) ]
&BEGIN_PROVIDER [ double precision, psi_selectors_coef_tc, (psi_selectors_size,2,N_states) ]
  implicit none
  BEGIN_DOC
  ! Determinants on which we apply <i|H|psi> for perturbation.
  END_DOC
  integer                        :: i,k

  do i=1,N_det_selectors
    do k=1,N_int
      psi_selectors(k,1,i) = psi_det_sorted_tc(k,1,i)
      psi_selectors(k,2,i) = psi_det_sorted_tc(k,2,i)
    enddo
  enddo
  do k=1,N_states
    do i=1,N_det_selectors
      psi_selectors_coef(i,k) = psi_coef_sorted_tc_gen(i,k)
!      psi_selectors_coef_tc(i,1,k) = psi_r_coef_sorted_bi_ortho(i,k)
!      psi_selectors_coef_tc(i,2,k) = psi_l_coef_sorted_bi_ortho(i,k)
      psi_selectors_coef_tc(i,1,k) = 1.d0
      psi_selectors_coef_tc(i,2,k) = 1.d0
    enddo
  enddo

END_PROVIDER

 BEGIN_PROVIDER [ double precision, psi_selectors_coef_transp, (N_states,psi_selectors_size) ]
&BEGIN_PROVIDER [ double precision, psi_selectors_coef_transp_tc, (N_states,2,psi_selectors_size) ]
  implicit none
  BEGIN_DOC
  ! Transposed psi_selectors
  END_DOC
  integer                        :: i,k

  do i=1,N_det_selectors
    do k=1,N_states
      psi_selectors_coef_transp(k,i) = psi_selectors_coef(i,k)
      psi_selectors_coef_transp_tc(k,1,i) = psi_selectors_coef_tc(i,1,k)
      psi_selectors_coef_transp_tc(k,2,i) = psi_selectors_coef_tc(i,2,k)
    enddo
  enddo
END_PROVIDER

 BEGIN_PROVIDER [ double precision, psi_selectors_rcoef_bi_orth_transp, (N_states, psi_det_size) ]
&BEGIN_PROVIDER [ double precision, psi_selectors_lcoef_bi_orth_transp, (N_states, psi_det_size) ]

  implicit none
  integer :: i, k

  psi_selectors_rcoef_bi_orth_transp = 0.d0
  psi_selectors_lcoef_bi_orth_transp = 0.d0

  print*,'N_det,N_det_selectors',N_det,N_det_selectors
  do i = 1, N_det_selectors
    do k = 1, N_states
      psi_selectors_rcoef_bi_orth_transp(k,i) = psi_r_coef_sorted_bi_ortho(i,k)
      psi_selectors_lcoef_bi_orth_transp(k,i) = psi_l_coef_sorted_bi_ortho(i,k)
    enddo
  enddo

END_PROVIDER

BEGIN_PROVIDER [ integer, psi_selectors_size ]
 implicit none
 psi_selectors_size = psi_det_size
END_PROVIDER

