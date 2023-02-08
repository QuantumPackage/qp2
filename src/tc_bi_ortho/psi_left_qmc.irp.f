
! ---

BEGIN_PROVIDER  [ double precision, psi_bitcleft_bilinear_matrix_values, (N_det,N_states) ]

  BEGIN_DOC
  ! Sparse coefficient matrix if the wave function is expressed in a bilinear form :
  !  $D_\alpha^\dagger.C.D_\beta$
  !
  ! Rows are $\alpha$ determinants and columns are $\beta$.
  !
  ! Order refers to psi_det
  END_DOC

  use bitmasks

  implicit none
  integer :: k, l

  if(N_det .eq. 1) then

    do l = 1, N_states
      psi_bitcleft_bilinear_matrix_values(1,l) = 1.d0
    enddo

  else

    do l = 1, N_states
      do k = 1, N_det
        psi_bitcleft_bilinear_matrix_values(k,l) = psi_l_coef_bi_ortho(k,l)
      enddo
    enddo

    PROVIDE psi_bilinear_matrix_order
    do l = 1, N_states
      call dset_order(psi_bitcleft_bilinear_matrix_values(1,l), psi_bilinear_matrix_order, N_det)
    enddo

  endif

END_PROVIDER

! ---

