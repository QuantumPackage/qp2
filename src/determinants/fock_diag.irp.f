subroutine build_fock_tmp(fock_diag_tmp,det_ref,Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
! Build the diagonal of the Fock matrix corresponding to a generator
! determinant. $F_{00}$ is $\langle i|H|i \rangle = E_0$.
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: det_ref(Nint,2)
  double precision, intent(out)  :: fock_diag_tmp(2,mo_num+1)

  integer                        :: occ(Nint*bit_kind_size,2)
  integer                        :: ne(2), i, j, ii, jj
  double precision               :: E0

  ! Compute Fock matrix diagonal elements
  call bitstring_to_list_ab(det_ref,occ,Ne,Nint)

  fock_diag_tmp = 0.d0
  E0 = 0.d0

  if (Ne(1) /= elec_alpha_num) then
    print *,  'Error in build_fock_tmp (alpha)', Ne(1), Ne(2)
    call debug_det(det_ref,N_int)
    stop -1
  endif
  if (Ne(2) /= elec_beta_num) then
    print *, 'Error in build_fock_tmp (beta)', Ne(1), Ne(2)
    call debug_det(det_ref,N_int)
    stop -1
  endif

  ! Occupied MOs
  do ii=1,elec_alpha_num
    i = occ(ii,1)
    fock_diag_tmp(1,i) = fock_diag_tmp(1,i) + mo_bi_ortho_tc_one_e(i,i)
    E0 = E0 + mo_bi_ortho_tc_one_e(i,i)
    do jj=1,elec_alpha_num
      j = occ(jj,1)
      if (i==j) cycle
      fock_diag_tmp(1,i) = fock_diag_tmp(1,i) + mo_bi_ortho_tc_two_e_jj_anti(i,j)
      E0 = E0 + 0.5d0*mo_bi_ortho_tc_two_e_jj_anti(i,j)
    enddo
    do jj=1,elec_beta_num
      j = occ(jj,2)
      fock_diag_tmp(1,i) = fock_diag_tmp(1,i) + mo_bi_ortho_tc_two_e_jj(i,j)
      E0 = E0 + mo_bi_ortho_tc_two_e_jj(i,j)
    enddo
  enddo
  do ii=1,elec_beta_num
    i = occ(ii,2)
    fock_diag_tmp(2,i) = fock_diag_tmp(2,i) + mo_bi_ortho_tc_one_e(i,i)
    E0 = E0 + mo_bi_ortho_tc_one_e(i,i)
    do jj=1,elec_beta_num
      j = occ(jj,2)
      if (i==j) cycle
      fock_diag_tmp(2,i) = fock_diag_tmp(2,i) + mo_bi_ortho_tc_two_e_jj_anti(i,j)
      E0 = E0 + 0.5d0*mo_bi_ortho_tc_two_e_jj_anti(i,j)
    enddo
    do jj=1,elec_alpha_num
      j = occ(jj,1)
      fock_diag_tmp(2,i) = fock_diag_tmp(2,i) + mo_bi_ortho_tc_two_e_jj(i,j)
    enddo
  enddo

  ! Virtual MOs
  do i=1,mo_num
    if (fock_diag_tmp(1,i) /= 0.d0) cycle
    fock_diag_tmp(1,i) = fock_diag_tmp(1,i) + mo_bi_ortho_tc_one_e(i,i)
    do jj=1,elec_alpha_num
      j = occ(jj,1)
      fock_diag_tmp(1,i) = fock_diag_tmp(1,i) + mo_bi_ortho_tc_two_e_jj_anti(i,j)
    enddo
    do jj=1,elec_beta_num
      j = occ(jj,2)
      fock_diag_tmp(1,i) = fock_diag_tmp(1,i) + mo_bi_ortho_tc_two_e_jj(i,j)
    enddo
  enddo
  do i=1,mo_num
    if (fock_diag_tmp(2,i) /= 0.d0) cycle
    fock_diag_tmp(2,i) = fock_diag_tmp(2,i) + mo_bi_ortho_tc_one_e(i,i)
    do jj=1,elec_beta_num
      j = occ(jj,2)
      fock_diag_tmp(2,i) = fock_diag_tmp(2,i) + mo_bi_ortho_tc_two_e_jj_anti(i,j)
    enddo
    do jj=1,elec_alpha_num
      j = occ(jj,1)
      fock_diag_tmp(2,i) = fock_diag_tmp(2,i) + mo_bi_ortho_tc_two_e_jj(i,j)
    enddo
  enddo

  fock_diag_tmp(1,mo_num+1) = E0
  fock_diag_tmp(2,mo_num+1) = E0

end
