
! ---

program combine_lr_tcscf

  BEGIN_DOC
  ! TODO : Put the documentation of the program here
  END_DOC

  implicit none

  my_grid_becke  = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  bi_ortho = .True.
  touch bi_ortho

  call comb_orbitals()

end

! ---

subroutine comb_orbitals()

  implicit none
  integer                       :: i, m, n, nn, mm
  double precision              :: accu_d, accu_nd
  double precision, allocatable :: R(:,:), L(:,:), Rnew(:,:), tmp(:,:), S(:,:)

  n  = ao_num
  m  = mo_num
  nn = elec_alpha_num
  mm = m - nn

  allocate(L(n,m), R(n,m), Rnew(n,m), S(m,m))
  L = mo_l_coef
  R = mo_r_coef

  call check_weighted_biorthog(n, m, ao_overlap, L, R, accu_d, accu_nd, S, .true.)

  allocate(tmp(n,nn))
  do i = 1, nn 
    tmp(1:n,i) = R(1:n,i)
  enddo
  call impose_weighted_orthog_svd(n, nn, ao_overlap, tmp)
  do i = 1, nn
    Rnew(1:n,i) = tmp(1:n,i)
  enddo
  deallocate(tmp)

  allocate(tmp(n,mm))
  do i = 1, mm
    tmp(1:n,i) = L(1:n,i+nn)
  enddo
  call impose_weighted_orthog_svd(n, mm, ao_overlap, tmp)
  do i = 1, mm
    Rnew(1:n,i+nn) = tmp(1:n,i)
  enddo
  deallocate(tmp)

  call check_weighted_biorthog(n, m, ao_overlap, Rnew, Rnew, accu_d, accu_nd, S, .true.)

  mo_r_coef = Rnew
  call ezfio_set_bi_ortho_mos_mo_r_coef(mo_r_coef)

  deallocate(L, R, Rnew, S)

end subroutine comb_orbitals

! ---

