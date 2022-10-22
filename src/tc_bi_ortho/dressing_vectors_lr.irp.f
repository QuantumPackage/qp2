
! ---

subroutine get_delta_bitc_right(psidet, psicoef, ndet, Nint, delta)
  
  BEGIN_DOC
  !
  ! delta(I) = < I_left | H_TC - H | Psi_right >
  !
  END_DOC

  use bitmasks

  implicit none

  integer,           intent(in)  :: ndet, Nint
  double precision,  intent(in)  :: psicoef(ndet)
  integer(bit_kind), intent(in)  :: psidet(Nint,2,ndet)
  double precision, intent(out)  :: delta(ndet)

  integer                        :: i, j
  double precision               :: h_mono, h_twoe, h_tot
  double precision               :: htc_mono, htc_twoe, htc_three, htc_tot
  double precision               :: delta_mat

  print *, ' get_delta_bitc_right ...'

  i = 1
  j = 1
  call htilde_mu_mat_bi_ortho(psidet(1,1,i), psidet(1,1,j), Nint, htc_mono, htc_twoe, htc_three, htc_tot)
  call hmat_bi_ortho         (psidet(1,1,i), psidet(1,1,j), Nint, h_mono, h_twoe, h_tot)

  delta = 0.d0
 !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(dynamic,8)   &
 !$OMP SHARED(delta, ndet, psidet, psicoef, Nint)      &
 !$OMP PRIVATE(i, j, delta_mat, h_mono, h_twoe, h_tot, &
 !$OMP        htc_mono, htc_twoe, htc_three, htc_tot)
  do i = 1, ndet
    do j = 1, ndet

      ! < I | Htilde | J >
      call htilde_mu_mat_bi_ortho(psidet(1,1,i), psidet(1,1,j), Nint, htc_mono, htc_twoe, htc_three, htc_tot)
      ! < I | H | J >
      call hmat_bi_ortho(psidet(1,1,i), psidet(1,1,j), Nint, h_mono, h_twoe, h_tot)

      delta_mat = htc_tot - h_tot

      delta(i) = delta(i) + psicoef(j) * delta_mat
    enddo
  enddo
 !$OMP END PARALLEL DO

end subroutine get_delta_bitc_right

! ---

subroutine get_htc_bitc_right(psidet, psicoef, ndet, Nint, delta)
  
  BEGIN_DOC
  !
  ! delta(I) = < I_left | H_TC | Psi_right >
  !
  END_DOC

  use bitmasks

  implicit none

  integer,           intent(in)  :: ndet, Nint
  double precision,  intent(in)  :: psicoef(ndet)
  integer(bit_kind), intent(in)  :: psidet(Nint,2,ndet)
  double precision, intent(out)  :: delta(ndet)

  integer                        :: i, j
  double precision               :: htc_mono, htc_twoe, htc_three, htc_tot

  print *, ' get_htc_bitc_right ...'

  i = 1
  j = 1
  call htilde_mu_mat_bi_ortho(psidet(1,1,i), psidet(1,1,j), Nint, htc_mono, htc_twoe, htc_three, htc_tot)

  delta = 0.d0
 !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(dynamic,8)   &
 !$OMP SHARED(delta, ndet, psidet, psicoef, Nint)      &
 !$OMP PRIVATE(i, j, htc_mono, htc_twoe, htc_three, htc_tot)
  do i = 1, ndet
    do j = 1, ndet

      ! < I | Htilde | J >
      call htilde_mu_mat_bi_ortho(psidet(1,1,i), psidet(1,1,j), Nint, htc_mono, htc_twoe, htc_three, htc_tot)

      delta(i) = delta(i) + psicoef(j) * htc_tot
    enddo
  enddo
 !$OMP END PARALLEL DO

end subroutine get_htc_bitc_right

! ---

subroutine get_h_bitc_right(psidet, psicoef, ndet, Nint, delta)
  
  BEGIN_DOC
  !
  ! delta(I) = < I_left | H | Psi_right >
  !
  END_DOC

  use bitmasks

  implicit none

  integer,           intent(in)  :: ndet, Nint
  double precision,  intent(in)  :: psicoef(ndet)
  integer(bit_kind), intent(in)  :: psidet(Nint,2,ndet)
  double precision, intent(out)  :: delta(ndet)

  integer                        :: i, j
  double precision               :: h_mono, h_twoe, h_tot

  print *, ' get_h_bitc_right ...'

  i = 1
  j = 1
  call hmat_bi_ortho(psidet(1,1,i), psidet(1,1,j), Nint, h_mono, h_twoe, h_tot)

  !double precision :: norm
  !norm = 0.d0
  !do i = 1, ndet
  !  norm += psicoef(i) * psicoef(i)
  !enddo
  !print*, ' norm = ', norm
 
  call hmat_bi_ortho(psidet(1,1,i), psidet(1,1,j), Nint, h_mono, h_twoe, h_tot)

  delta = 0.d0
! !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(dynamic,8) &
! !$OMP SHARED(delta, ndet, psidet, psicoef, Nint)    &
! !$OMP PRIVATE(i, j, h_mono, h_twoe, h_tot)
  do i = 1, ndet
    do j = 1, ndet

      ! < I | H | J >
      call hmat_bi_ortho(psidet(1,1,i), psidet(1,1,j), Nint, h_mono, h_twoe, h_tot)

      delta(i) = delta(i) + psicoef(j) * h_tot
    enddo
  enddo
! !$OMP END PARALLEL DO

end subroutine get_h_bitc_right

! ---

