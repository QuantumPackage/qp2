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
!
!  call H_apply_mp2(pt2, norm_pert, H_pert_diag,  N_st)
!  print *,  'N_det    = ', N_det
!  print *,  'N_states = ', N_states
!  print *,  'MP2      = ', pt2
!  print *,  'E        = ', E_old
!  print *,  'E+MP2    = ', E_old+pt2
!  print *, '----'
!
  deallocate(pt2,norm_pert)
end


subroutine compute_mp2(emp2)
  implicit none
  double precision, intent(out) :: emp2
  integer :: i,j,a,b, k, sze
  double precision, allocatable :: ijab(:,:), ijba(:,:), ijab_anti, e(:)
  double precision :: denom

  PROVIDE mo_two_e_integrals_in_map

  allocate(ijab(mo_num,mo_num), ijba(mo_num,mo_num), e(mo_num))
  do i=1,mo_num
    e(i) = Fock_matrix_diag_mo(i)
  enddo

  emp2 = 0.d0

  do j=1,elec_beta_num
    if (list_core_reverse(j) > 0) cycle

    do b=elec_beta_num+1,mo_num
      if (list_del_reverse(b) > 0) cycle

      emp2 = emp2 + Fock_matrix_mo(b,j)*Fock_matrix_mo(b,j) / (e(j) - e(b))

      call get_mo_two_e_integrals_i1j1(j,b,mo_num,ijab,mo_integrals_map)
      call get_mo_two_e_integrals_ij(b,j,mo_num,ijba,mo_integrals_map)

      do i=1,elec_alpha_num
        if (list_core_reverse(i) > 0) cycle

        do a=elec_alpha_num+1,mo_num
          if (list_del_reverse(a) > 0) cycle
          denom = 1.d0 / (e(i) + e(j) - e(a) - e(b))


          emp2 = emp2 + ijab(i,a)*ijab(i,a) * denom

        end do

      end do

      do i=j+1,elec_beta_num
        if (list_core_reverse(i) > 0) cycle

        do a=b+1,mo_num 
          if (list_del_reverse(a) > 0) cycle
          denom = 1.d0 / (e(i) + e(j) - e(a) - e(b))

          ijab_anti = ijab(i,a) - ijba(i,a)

          emp2 = emp2 + ijab_anti*ijab_anti * denom

        end do

      end do

    end do

  end do

  do j=1,elec_alpha_num
    if (list_core_reverse(j) > 0) cycle

    do b=elec_alpha_num+1,mo_num
      if (list_del_reverse(b) > 0) cycle

      emp2 = emp2 + Fock_matrix_mo(b,j)*Fock_matrix_mo(b,j) / (e(j) - e(b))

      call get_mo_two_e_integrals_i1j1(j,b,mo_num,ijab,mo_integrals_map)
      call get_mo_two_e_integrals_ij(b,j,mo_num,ijba,mo_integrals_map)

      do i=j+1,elec_alpha_num
        if (list_core_reverse(i) > 0) cycle

        do a=b+1,mo_num
          if (list_del_reverse(a) > 0) cycle
          denom = 1.d0 / (e(i) + e(j) - e(a) - e(b))

          ijab_anti = ijab(i,a) - ijba(i,a)

          emp2 = emp2 + ijab_anti*ijab_anti * denom

        end do

      end do


    end do

  end do
end
