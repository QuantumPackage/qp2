program mp2
  call run
end

subroutine run
  implicit none
  double precision, allocatable  :: pt2(:), norm_pert(:)
  double precision               :: H_pert_diag, E_old
  integer                        :: N_st, iter
  PROVIDE all_mo_integrals  Fock_matrix_diag_mo H_apply_buffer_allocated
  N_st = N_states
  allocate (pt2(N_st), norm_pert(N_st))
  E_old = HF_energy

  pt2 = 0.d0
  call compute_mp2(pt2(1))
  print *,  'N_det    = ', N_det
  print *,  'N_states = ', N_states
  print *,  'MP2      = ', pt2
  print *,  'E        = ', E_old
  print *,  'E+MP2    = ', E_old+pt2(1)

  deallocate(pt2,norm_pert)
end


subroutine compute_mp2(emp2)
  implicit none
  double precision, intent(out) :: emp2
  integer :: i,j,a,b, k, sze
  double precision, allocatable :: ijab(:,:), ijba(:,:), ijab_anti, e(:,:), ei, ej, ea, eb
  double precision :: denom
  double precision :: aaaa, bbbb, abab, singles

  PROVIDE mo_two_e_integrals_in_map

  allocate(ijab(mo_num,mo_num), ijba(mo_num,mo_num), e(mo_num,3))
  do i=1,mo_num
    e(i,1) = Fock_matrix_mo_alpha(i,i)
    e(i,2) = Fock_matrix_mo_beta(i,i)
    e(i,3) = Fock_matrix_diag_mo(i)
  enddo

  singles = 0.d0
  aaaa = 0.d0
  bbbb = 0.d0
  abab = 0.d0


!  For OO and VV blocks, we take the common Fock matrix.
!  Otherwise, we use different alpha and beta Fock matrices

  do j=1,elec_beta_num
    if (list_core_reverse(j) > 0) cycle
    ej = e(j,3)

    do b=elec_beta_num+1,mo_num
      if (list_del_reverse(b) > 0) cycle
      if (b>elec_alpha_num) then
        eb = e(b,3)
      else
        eb = e(b,2)
      endif

      singles = singles + Fock_matrix_mo(b,j)*Fock_matrix_mo(b,j) / (ej - eb)

      call get_mo_two_e_integrals_i1j1(j,b,mo_num,ijab,mo_integrals_map)
      call get_mo_two_e_integrals_ij(b,j,mo_num,ijba,mo_integrals_map)

      do i=1,elec_alpha_num
        if (list_core_reverse(i) > 0) cycle
        if (i>elec_beta_num) then
          ei = e(i,1)
        else
          ei = e(i,3)
        endif

        do a=elec_alpha_num+1,mo_num
          if (list_del_reverse(a) > 0) cycle
          ea = e(a,3)
          denom = 1.d0 / (ei + ej - ea - eb)


          abab = abab + ijab(i,a)*ijab(i,a) * denom

        end do

      end do

      do i=j+1,elec_beta_num
        if (list_core_reverse(i) > 0) cycle
        ei = e(i,2)

        do a=b+1,mo_num
          if (list_del_reverse(a) > 0) cycle
          if (a>elec_alpha_num) then
            ea = e(a,3)
          else
            ea = e(a,2)
          endif
          denom = 1.d0 / (ei + ej - ea - eb)

          ijab_anti = ijab(i,a) - ijba(i,a)

          bbbb = bbbb + ijab_anti*ijab_anti * denom

        end do

      end do

    end do

  end do

  do j=1,elec_alpha_num
    if (list_core_reverse(j) > 0) cycle
    if (b>elec_beta_num) then
      ej = e(j,1)
    else
      ej = e(j,3)
    endif

    do b=elec_alpha_num+1,mo_num
      if (list_del_reverse(b) > 0) cycle
      eb = e(b,3)

      singles = singles + Fock_matrix_mo(b,j)*Fock_matrix_mo(b,j) / (ej - eb)

      call get_mo_two_e_integrals_i1j1(j,b,mo_num,ijab,mo_integrals_map)
      call get_mo_two_e_integrals_ij(b,j,mo_num,ijba,mo_integrals_map)

      do i=j+1,elec_alpha_num
        if (list_core_reverse(i) > 0) cycle
        if (i>elec_beta_num) then
          ei = e(i,1)
        else
          ei = e(i,3)
        endif

        do a=b+1,mo_num
          if (list_del_reverse(a) > 0) cycle
          ea = e(a,3)
          denom = 1.d0 / (ei + ej - ea - eb)

          ijab_anti = ijab(i,a) - ijba(i,a)

          aaaa = aaaa + ijab_anti*ijab_anti * denom

        end do

      end do


    end do

  end do

  print *, '       aaaa  correlation energy = ', aaaa, 'au'
  print *, '       abab  correlation energy = ', abab, 'au'
  print *, '       bbbb  correlation energy = ', bbbb, 'au'
  print *, '       non-Brillouin singles    = ', singles, 'au'

  emp2 = singles + aaaa + bbbb + abab
end
