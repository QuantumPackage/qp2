
 BEGIN_PROVIDER [double precision, TC_HF_energy        ]
&BEGIN_PROVIDER [double precision, TC_HF_one_e_energy  ]
&BEGIN_PROVIDER [double precision, TC_HF_two_e_energy  ]
&BEGIN_PROVIDER [double precision, TC_HF_three_e_energy]

  BEGIN_DOC
  ! TC Hartree-Fock energy containing the nuclear repulsion, and its one- and two-body components.
  END_DOC

  implicit none
  integer          :: i, j
  double precision :: t0, t1

  PROVIDE mo_l_coef mo_r_coef
  PROVIDE two_e_tc_integral_alpha two_e_tc_integral_beta

  TC_HF_energy = nuclear_repulsion
  TC_HF_one_e_energy = 0.d0
  TC_HF_two_e_energy = 0.d0

  do j = 1, ao_num
    do i = 1, ao_num
      TC_HF_two_e_energy += 0.5d0 * ( two_e_tc_integral_alpha(i,j) * TCSCF_density_matrix_ao_alpha(i,j) &
                                    + two_e_tc_integral_beta (i,j) * TCSCF_density_matrix_ao_beta (i,j) )
      TC_HF_one_e_energy += ao_one_e_integrals_tc_tot(i,j) &
                          * (TCSCF_density_matrix_ao_alpha(i,j) + TCSCF_density_matrix_ao_beta (i,j) )
    enddo
  enddo

  if((three_body_h_tc .eq. .False.) .and. (.not. noL_standard)) then
    TC_HF_three_e_energy = 0.d0
  else
    TC_HF_three_e_energy = noL_0e
  endif

  TC_HF_energy += TC_HF_one_e_energy + TC_HF_two_e_energy + TC_HF_three_e_energy

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, diag_three_elem_hf]

  BEGIN_DOC
  ! 
  ! < Phi_left | L | Phi_right >
  !
  !
  ! if three_body_h_tc == false and noL_standard == true ==> do a normal ordering
  !
  ! todo
  ! this should be equivalent to
  !     three_body_h_tc == true and noL_standard == false
  !
  ! if three_body_h_tc == false and noL_standard == false ==> this is equal to 0
  !
  END_DOC

  implicit none
  integer                       :: i, j, k, ipoint, mm
  double precision              :: contrib, weight, four_third, one_third, two_third, exchange_int_231
  double precision              :: integral_aaa, hthree, integral_aab, integral_abb, integral_bbb
  double precision, allocatable :: tmp(:)
  double precision, allocatable :: tmp_L(:,:), tmp_R(:,:)
  double precision, allocatable :: tmp_M(:,:), tmp_S(:), tmp_O(:), tmp_J(:,:)
  double precision, allocatable :: tmp_M_priv(:,:), tmp_S_priv(:), tmp_O_priv(:), tmp_J_priv(:,:)

  PROVIDE mo_l_coef mo_r_coef

  if(.not. three_body_h_tc) then

   if(noL_standard) then
      PROVIDE noL_0e
      diag_three_elem_hf = noL_0e
    else
      diag_three_elem_hf = 0.d0
    endif

  else

    PROVIDE int2_grad1_u12_bimo_t
    PROVIDE mos_l_in_r_array_transp
    PROVIDE mos_r_in_r_array_transp

    if(elec_alpha_num .eq. elec_beta_num) then
  
      allocate(tmp(elec_beta_num))
      allocate(tmp_L(n_points_final_grid,3), tmp_R(n_points_final_grid,3))
  
      !$OMP PARALLEL                                                 &
      !$OMP DEFAULT(NONE)                                            &
      !$OMP PRIVATE(j, i, ipoint, tmp_L, tmp_R)                      &
      !$OMP SHARED(elec_beta_num, n_points_final_grid,               & 
      !$OMP        mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
      !$OMP        int2_grad1_u12_bimo_t, tmp, final_weight_at_r_vector)
  
      !$OMP DO
      do j = 1, elec_beta_num
  
        tmp_L = 0.d0
        tmp_R = 0.d0
        do i = 1, elec_beta_num
          do ipoint = 1, n_points_final_grid
  
            tmp_L(ipoint,1) = tmp_L(ipoint,1) + int2_grad1_u12_bimo_t(ipoint,1,j,i) * mos_l_in_r_array_transp(ipoint,i)
            tmp_L(ipoint,2) = tmp_L(ipoint,2) + int2_grad1_u12_bimo_t(ipoint,2,j,i) * mos_l_in_r_array_transp(ipoint,i)
            tmp_L(ipoint,3) = tmp_L(ipoint,3) + int2_grad1_u12_bimo_t(ipoint,3,j,i) * mos_l_in_r_array_transp(ipoint,i)
  
            tmp_R(ipoint,1) = tmp_R(ipoint,1) + int2_grad1_u12_bimo_t(ipoint,1,i,j) * mos_r_in_r_array_transp(ipoint,i)
            tmp_R(ipoint,2) = tmp_R(ipoint,2) + int2_grad1_u12_bimo_t(ipoint,2,i,j) * mos_r_in_r_array_transp(ipoint,i)
            tmp_R(ipoint,3) = tmp_R(ipoint,3) + int2_grad1_u12_bimo_t(ipoint,3,i,j) * mos_r_in_r_array_transp(ipoint,i)
          enddo
        enddo
  
        tmp(j) = 0.d0
        do ipoint = 1, n_points_final_grid
          tmp(j) = tmp(j) + final_weight_at_r_vector(ipoint) * (tmp_L(ipoint,1)*tmp_R(ipoint,1) + tmp_L(ipoint,2)*tmp_R(ipoint,2) + tmp_L(ipoint,3)*tmp_R(ipoint,3))
        enddo
      enddo ! j
      !$OMP END DO
      !$OMP END PARALLEL
  
      diag_three_elem_hf = -2.d0 * sum(tmp)
  
      deallocate(tmp)
      deallocate(tmp_L, tmp_R)
  
      ! ---
  
      allocate(tmp_O(n_points_final_grid), tmp_J(n_points_final_grid,3))
      tmp_O = 0.d0
      tmp_J = 0.d0
  
      !$OMP PARALLEL                                                  &
      !$OMP DEFAULT(NONE)                                             &
      !$OMP PRIVATE(i, ipoint, tmp_O_priv, tmp_J_priv)                &
      !$OMP SHARED(elec_beta_num, n_points_final_grid,                & 
      !$OMP        mos_l_in_r_array_transp, mos_r_in_r_array_transp,  &
      !$OMP        int2_grad1_u12_bimo_t, tmp_O, tmp_J)
  
      allocate(tmp_O_priv(n_points_final_grid), tmp_J_priv(n_points_final_grid,3))
      tmp_O_priv = 0.d0
      tmp_J_priv = 0.d0
    
      !$OMP DO 
      do i = 1, elec_beta_num
        do ipoint = 1, n_points_final_grid
          tmp_O_priv(ipoint)   = tmp_O_priv(ipoint)   + mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,i)
          tmp_J_priv(ipoint,1) = tmp_J_priv(ipoint,1) + int2_grad1_u12_bimo_t(ipoint,1,i,i)
          tmp_J_priv(ipoint,2) = tmp_J_priv(ipoint,2) + int2_grad1_u12_bimo_t(ipoint,2,i,i)
          tmp_J_priv(ipoint,3) = tmp_J_priv(ipoint,3) + int2_grad1_u12_bimo_t(ipoint,3,i,i)
        enddo
      enddo
      !$OMP END DO NOWAIT
  
      !$OMP CRITICAL
      tmp_O = tmp_O + tmp_O_priv
      tmp_J = tmp_J + tmp_J_priv
      !$OMP END CRITICAL
  
      deallocate(tmp_O_priv, tmp_J_priv)
      !$OMP END PARALLEL
  
      allocate(tmp_M(n_points_final_grid,3), tmp_S(n_points_final_grid))
      tmp_M = 0.d0
      tmp_S = 0.d0
  
      !$OMP PARALLEL                                                 &
      !$OMP DEFAULT(NONE)                                            &
      !$OMP PRIVATE(i, j, ipoint, tmp_M_priv, tmp_S_priv)            &
      !$OMP SHARED(elec_beta_num, n_points_final_grid,               & 
      !$OMP        mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
      !$OMP        int2_grad1_u12_bimo_t, tmp_M, tmp_S)
  
      allocate(tmp_M_priv(n_points_final_grid,3), tmp_S_priv(n_points_final_grid))
      tmp_M_priv = 0.d0
      tmp_S_priv = 0.d0
    
      !$OMP DO COLLAPSE(2)
      do i = 1, elec_beta_num
        do j = 1, elec_beta_num
          do ipoint = 1, n_points_final_grid
  
            tmp_M_priv(ipoint,1) = tmp_M_priv(ipoint,1) + int2_grad1_u12_bimo_t(ipoint,1,j,i) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,j)
            tmp_M_priv(ipoint,2) = tmp_M_priv(ipoint,2) + int2_grad1_u12_bimo_t(ipoint,2,j,i) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,j)
            tmp_M_priv(ipoint,3) = tmp_M_priv(ipoint,3) + int2_grad1_u12_bimo_t(ipoint,3,j,i) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,j)
  
            tmp_S_priv(ipoint)   = tmp_S_priv(ipoint)   + int2_grad1_u12_bimo_t(ipoint,1,i,j) * int2_grad1_u12_bimo_t(ipoint,1,j,i) &
                                                        + int2_grad1_u12_bimo_t(ipoint,2,i,j) * int2_grad1_u12_bimo_t(ipoint,2,j,i) &
                                                        + int2_grad1_u12_bimo_t(ipoint,3,i,j) * int2_grad1_u12_bimo_t(ipoint,3,j,i)
          enddo
        enddo
      enddo
      !$OMP END DO NOWAIT
  
      !$OMP CRITICAL
      tmp_M = tmp_M + tmp_M_priv
      tmp_S = tmp_S + tmp_S_priv
      !$OMP END CRITICAL
  
      deallocate(tmp_M_priv, tmp_S_priv)
      !$OMP END PARALLEL
  
      allocate(tmp(n_points_final_grid))
  
      do ipoint = 1, n_points_final_grid
  
        tmp_S(ipoint) = 2.d0 * (tmp_J(ipoint,1)*tmp_J(ipoint,1) + tmp_J(ipoint,2)*tmp_J(ipoint,2) + tmp_J(ipoint,3)*tmp_J(ipoint,3)) - tmp_S(ipoint)
  
        tmp(ipoint) = final_weight_at_r_vector(ipoint) * ( tmp_O(ipoint) * tmp_S(ipoint)              &
                                                         - 2.d0 * ( tmp_J(ipoint,1) * tmp_M(ipoint,1) &
                                                                  + tmp_J(ipoint,2) * tmp_M(ipoint,2) &
                                                                  + tmp_J(ipoint,3) * tmp_M(ipoint,3)))
      enddo
  
      diag_three_elem_hf = diag_three_elem_hf -2.d0 * (sum(tmp))
  
      deallocate(tmp)
  
    else ! elec_alpha_num .neq. elec_beta_num
  
      allocate(tmp(elec_alpha_num))
      allocate(tmp_L(n_points_final_grid,3), tmp_R(n_points_final_grid,3))
  
      !$OMP PARALLEL                                                   &
      !$OMP DEFAULT(NONE)                                              &
      !$OMP PRIVATE(j, i, ipoint, tmp_L, tmp_R)                        &
      !$OMP SHARED(elec_beta_num, elec_alpha_num, n_points_final_grid, & 
      !$OMP        mos_l_in_r_array_transp, mos_r_in_r_array_transp,   &
      !$OMP        int2_grad1_u12_bimo_t, tmp, final_weight_at_r_vector)
  
      !$OMP DO
      do j = 1, elec_beta_num
  
        tmp_L = 0.d0
        tmp_R = 0.d0
        do i = elec_beta_num+1, elec_alpha_num
          do ipoint = 1, n_points_final_grid
  
            tmp_L(ipoint,1) = tmp_L(ipoint,1) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,1,j,i) * mos_l_in_r_array_transp(ipoint,i)
            tmp_L(ipoint,2) = tmp_L(ipoint,2) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,2,j,i) * mos_l_in_r_array_transp(ipoint,i)
            tmp_L(ipoint,3) = tmp_L(ipoint,3) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,3,j,i) * mos_l_in_r_array_transp(ipoint,i)
  
            tmp_R(ipoint,1) = tmp_R(ipoint,1) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,1,i,j) * mos_r_in_r_array_transp(ipoint,i)
            tmp_R(ipoint,2) = tmp_R(ipoint,2) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,2,i,j) * mos_r_in_r_array_transp(ipoint,i)
            tmp_R(ipoint,3) = tmp_R(ipoint,3) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,3,i,j) * mos_r_in_r_array_transp(ipoint,i)
          enddo
        enddo
  
        tmp(j) = 0.d0
        do ipoint = 1, n_points_final_grid
          tmp(j) = tmp(j) + final_weight_at_r_vector(ipoint) * (tmp_L(ipoint,1)*tmp_R(ipoint,1) + tmp_L(ipoint,2)*tmp_R(ipoint,2) + tmp_L(ipoint,3)*tmp_R(ipoint,3))
        enddo
  
        do i = 1, elec_beta_num
          do ipoint = 1, n_points_final_grid
  
            tmp_L(ipoint,1) = tmp_L(ipoint,1) + int2_grad1_u12_bimo_t(ipoint,1,j,i) * mos_l_in_r_array_transp(ipoint,i)
            tmp_L(ipoint,2) = tmp_L(ipoint,2) + int2_grad1_u12_bimo_t(ipoint,2,j,i) * mos_l_in_r_array_transp(ipoint,i)
            tmp_L(ipoint,3) = tmp_L(ipoint,3) + int2_grad1_u12_bimo_t(ipoint,3,j,i) * mos_l_in_r_array_transp(ipoint,i)
  
            tmp_R(ipoint,1) = tmp_R(ipoint,1) + int2_grad1_u12_bimo_t(ipoint,1,i,j) * mos_r_in_r_array_transp(ipoint,i)
            tmp_R(ipoint,2) = tmp_R(ipoint,2) + int2_grad1_u12_bimo_t(ipoint,2,i,j) * mos_r_in_r_array_transp(ipoint,i)
            tmp_R(ipoint,3) = tmp_R(ipoint,3) + int2_grad1_u12_bimo_t(ipoint,3,i,j) * mos_r_in_r_array_transp(ipoint,i)
          enddo
        enddo
  
        do ipoint = 1, n_points_final_grid
          tmp(j) = tmp(j) + final_weight_at_r_vector(ipoint) * (tmp_L(ipoint,1)*tmp_R(ipoint,1) + tmp_L(ipoint,2)*tmp_R(ipoint,2) + tmp_L(ipoint,3)*tmp_R(ipoint,3))
        enddo
      enddo ! j
      !$OMP END DO
      !$OMP END PARALLEL
  
      ! ---
  
      !$OMP PARALLEL                                                   &
      !$OMP DEFAULT(NONE)                                              &
      !$OMP PRIVATE(j, i, ipoint, tmp_L, tmp_R)                        &
      !$OMP SHARED(elec_beta_num, elec_alpha_num, n_points_final_grid, & 
      !$OMP        mos_l_in_r_array_transp, mos_r_in_r_array_transp,   &
      !$OMP        int2_grad1_u12_bimo_t, tmp, final_weight_at_r_vector)
  
      !$OMP DO
      do j = elec_beta_num+1, elec_alpha_num
  
        tmp_L = 0.d0
        tmp_R = 0.d0
        do i = 1, elec_alpha_num
          do ipoint = 1, n_points_final_grid
            tmp_L(ipoint,1) = tmp_L(ipoint,1) + int2_grad1_u12_bimo_t(ipoint,1,j,i) * mos_l_in_r_array_transp(ipoint,i)
            tmp_L(ipoint,2) = tmp_L(ipoint,2) + int2_grad1_u12_bimo_t(ipoint,2,j,i) * mos_l_in_r_array_transp(ipoint,i)
            tmp_L(ipoint,3) = tmp_L(ipoint,3) + int2_grad1_u12_bimo_t(ipoint,3,j,i) * mos_l_in_r_array_transp(ipoint,i)
  
            tmp_R(ipoint,1) = tmp_R(ipoint,1) + int2_grad1_u12_bimo_t(ipoint,1,i,j) * mos_r_in_r_array_transp(ipoint,i)
            tmp_R(ipoint,2) = tmp_R(ipoint,2) + int2_grad1_u12_bimo_t(ipoint,2,i,j) * mos_r_in_r_array_transp(ipoint,i)
            tmp_R(ipoint,3) = tmp_R(ipoint,3) + int2_grad1_u12_bimo_t(ipoint,3,i,j) * mos_r_in_r_array_transp(ipoint,i)
          enddo
        enddo
  
        tmp(j) = 0.d0
        do ipoint = 1, n_points_final_grid
          tmp(j) = tmp(j) + 0.5d0 * final_weight_at_r_vector(ipoint) * (tmp_L(ipoint,1)*tmp_R(ipoint,1) + tmp_L(ipoint,2)*tmp_R(ipoint,2) + tmp_L(ipoint,3)*tmp_R(ipoint,3))
        enddo
      enddo ! j
      !$OMP END DO
      !$OMP END PARALLEL
  
      diag_three_elem_hf = -2.d0 * sum(tmp)
  
      deallocate(tmp)
      deallocate(tmp_L, tmp_R)
  
      ! ---
  
      allocate(tmp_O(n_points_final_grid), tmp_J(n_points_final_grid,3))
      tmp_O = 0.d0
      tmp_J = 0.d0
  
      !$OMP PARALLEL                                                   &
      !$OMP DEFAULT(NONE)                                              &
      !$OMP PRIVATE(i, ipoint, tmp_O_priv, tmp_J_priv)                 &
      !$OMP SHARED(elec_beta_num, elec_alpha_num, n_points_final_grid, & 
      !$OMP        mos_l_in_r_array_transp, mos_r_in_r_array_transp,   &
      !$OMP        int2_grad1_u12_bimo_t, tmp_O, tmp_J)
  
      allocate(tmp_O_priv(n_points_final_grid), tmp_J_priv(n_points_final_grid,3))
      tmp_O_priv = 0.d0
      tmp_J_priv = 0.d0
    
      !$OMP DO 
      do i = 1, elec_beta_num
        do ipoint = 1, n_points_final_grid
          tmp_O_priv(ipoint)   = tmp_O_priv(ipoint)   + mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,i)
          tmp_J_priv(ipoint,1) = tmp_J_priv(ipoint,1) + int2_grad1_u12_bimo_t(ipoint,1,i,i)
          tmp_J_priv(ipoint,2) = tmp_J_priv(ipoint,2) + int2_grad1_u12_bimo_t(ipoint,2,i,i)
          tmp_J_priv(ipoint,3) = tmp_J_priv(ipoint,3) + int2_grad1_u12_bimo_t(ipoint,3,i,i)
        enddo
      enddo
      !$OMP END DO NOWAIT
  
      !$OMP DO 
      do i = elec_beta_num+1, elec_alpha_num
        do ipoint = 1, n_points_final_grid
          tmp_O_priv(ipoint)   = tmp_O_priv(ipoint)   + 0.5d0 * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,i)
          tmp_J_priv(ipoint,1) = tmp_J_priv(ipoint,1) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,1,i,i)
          tmp_J_priv(ipoint,2) = tmp_J_priv(ipoint,2) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,2,i,i)
          tmp_J_priv(ipoint,3) = tmp_J_priv(ipoint,3) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,3,i,i)
        enddo
      enddo
      !$OMP END DO NOWAIT
  
      !$OMP CRITICAL
      tmp_O = tmp_O + tmp_O_priv
      tmp_J = tmp_J + tmp_J_priv
      !$OMP END CRITICAL
  
      deallocate(tmp_O_priv, tmp_J_priv)
      !$OMP END PARALLEL
  
      ! ---
  
      allocate(tmp_M(n_points_final_grid,3), tmp_S(n_points_final_grid))
      tmp_M = 0.d0
      tmp_S = 0.d0
  
      !$OMP PARALLEL                                                   &
      !$OMP DEFAULT(NONE)                                              &
      !$OMP PRIVATE(i, j, ipoint, tmp_M_priv, tmp_S_priv)              &
      !$OMP SHARED(elec_beta_num, elec_alpha_num, n_points_final_grid, & 
      !$OMP        mos_l_in_r_array_transp, mos_r_in_r_array_transp,   &
      !$OMP        int2_grad1_u12_bimo_t, tmp_M, tmp_S)
  
      allocate(tmp_M_priv(n_points_final_grid,3), tmp_S_priv(n_points_final_grid))
      tmp_M_priv = 0.d0
      tmp_S_priv = 0.d0
    
      !$OMP DO COLLAPSE(2)
      do i = 1, elec_beta_num
        do j = 1, elec_beta_num
          do ipoint = 1, n_points_final_grid
  
            tmp_M_priv(ipoint,1) = tmp_M_priv(ipoint,1) + int2_grad1_u12_bimo_t(ipoint,1,j,i) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,j)
            tmp_M_priv(ipoint,2) = tmp_M_priv(ipoint,2) + int2_grad1_u12_bimo_t(ipoint,2,j,i) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,j)
            tmp_M_priv(ipoint,3) = tmp_M_priv(ipoint,3) + int2_grad1_u12_bimo_t(ipoint,3,j,i) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,j)
  
            tmp_S_priv(ipoint)   = tmp_S_priv(ipoint)   + int2_grad1_u12_bimo_t(ipoint,1,i,j) * int2_grad1_u12_bimo_t(ipoint,1,j,i) &
                                                        + int2_grad1_u12_bimo_t(ipoint,2,i,j) * int2_grad1_u12_bimo_t(ipoint,2,j,i) &
                                                        + int2_grad1_u12_bimo_t(ipoint,3,i,j) * int2_grad1_u12_bimo_t(ipoint,3,j,i)
          enddo
        enddo
      enddo
      !$OMP END DO NOWAIT
  
      !$OMP DO COLLAPSE(2)
      do i = elec_beta_num+1, elec_alpha_num
        do j = 1, elec_beta_num
          do ipoint = 1, n_points_final_grid
  
            tmp_M_priv(ipoint,1) = tmp_M_priv(ipoint,1) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,1,j,i) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,j)
            tmp_M_priv(ipoint,2) = tmp_M_priv(ipoint,2) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,2,j,i) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,j)
            tmp_M_priv(ipoint,3) = tmp_M_priv(ipoint,3) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,3,j,i) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,j)
  
            tmp_M_priv(ipoint,1) = tmp_M_priv(ipoint,1) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,1,i,j) * mos_l_in_r_array_transp(ipoint,j) * mos_r_in_r_array_transp(ipoint,i)
            tmp_M_priv(ipoint,2) = tmp_M_priv(ipoint,2) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,2,i,j) * mos_l_in_r_array_transp(ipoint,j) * mos_r_in_r_array_transp(ipoint,i)
            tmp_M_priv(ipoint,3) = tmp_M_priv(ipoint,3) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,3,i,j) * mos_l_in_r_array_transp(ipoint,j) * mos_r_in_r_array_transp(ipoint,i)
  
            tmp_S_priv(ipoint)   = tmp_S_priv(ipoint)   + int2_grad1_u12_bimo_t(ipoint,1,i,j) * int2_grad1_u12_bimo_t(ipoint,1,j,i) &
                                                        + int2_grad1_u12_bimo_t(ipoint,2,i,j) * int2_grad1_u12_bimo_t(ipoint,2,j,i) &
                                                        + int2_grad1_u12_bimo_t(ipoint,3,i,j) * int2_grad1_u12_bimo_t(ipoint,3,j,i)
          enddo
        enddo
      enddo
      !$OMP END DO NOWAIT
  
      !$OMP DO COLLAPSE(2)
      do i = elec_beta_num+1, elec_alpha_num
        do j = elec_beta_num+1, elec_alpha_num
          do ipoint = 1, n_points_final_grid
  
            tmp_M_priv(ipoint,1) = tmp_M_priv(ipoint,1) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,1,j,i) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,j)
            tmp_M_priv(ipoint,2) = tmp_M_priv(ipoint,2) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,2,j,i) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,j)
            tmp_M_priv(ipoint,3) = tmp_M_priv(ipoint,3) + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,3,j,i) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,j)
  
            tmp_S_priv(ipoint)   = tmp_S_priv(ipoint)   + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,1,i,j) * int2_grad1_u12_bimo_t(ipoint,1,j,i) &
                                                        + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,2,i,j) * int2_grad1_u12_bimo_t(ipoint,2,j,i) &
                                                        + 0.5d0 * int2_grad1_u12_bimo_t(ipoint,3,i,j) * int2_grad1_u12_bimo_t(ipoint,3,j,i)
          enddo
        enddo
      enddo
      !$OMP END DO NOWAIT
  
      !$OMP CRITICAL
      tmp_M = tmp_M + tmp_M_priv
      tmp_S = tmp_S + tmp_S_priv
      !$OMP END CRITICAL
  
      deallocate(tmp_M_priv, tmp_S_priv)
      !$OMP END PARALLEL
  
      allocate(tmp(n_points_final_grid))
  
      do ipoint = 1, n_points_final_grid
  
        tmp_S(ipoint) = 2.d0 * (tmp_J(ipoint,1)*tmp_J(ipoint,1) + tmp_J(ipoint,2)*tmp_J(ipoint,2) + tmp_J(ipoint,3)*tmp_J(ipoint,3)) - tmp_S(ipoint)
  
        tmp(ipoint) = final_weight_at_r_vector(ipoint) * ( tmp_O(ipoint) * tmp_S(ipoint)              &
                                                         - 2.d0 * ( tmp_J(ipoint,1) * tmp_M(ipoint,1) &
                                                                  + tmp_J(ipoint,2) * tmp_M(ipoint,2) &
                                                                  + tmp_J(ipoint,3) * tmp_M(ipoint,3)))
      enddo
  
      diag_three_elem_hf = diag_three_elem_hf - 2.d0 * (sum(tmp))
  
      deallocate(tmp)
  
    endif ! alpha/beta condition

  endif ! three_body_h_tc

END_PROVIDER 

! ---

