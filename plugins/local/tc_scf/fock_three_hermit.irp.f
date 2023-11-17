
! ---

BEGIN_PROVIDER [ double precision, fock_3_mat, (mo_num, mo_num)] 

  implicit none
  integer :: i,j
  double precision :: contrib

  fock_3_mat = 0.d0
  if(.not.bi_ortho .and. three_body_h_tc) then

    call give_fock_ia_three_e_total(1, 1, contrib)
    !!  !$OMP PARALLEL                  &
    !!  !$OMP DEFAULT (NONE)            &
    !!  !$OMP PRIVATE (i,j,m,integral) & 
    !!  !$OMP SHARED (mo_num,three_body_3_index)
    !!  !$OMP DO SCHEDULE (guided) COLLAPSE(3)
    do i = 1, mo_num
      do j = 1, mo_num
        call give_fock_ia_three_e_total(j,i,contrib)
        fock_3_mat(j,i) = -contrib
       enddo
     enddo
    !else if(bi_ortho.and.three_body_h_tc) then
    !!  !$OMP END DO
    !!  !$OMP END PARALLEL
    !!  do i = 1, mo_num
    !!   do j = 1, i-1
    !!    mat_three(j,i) = mat_three(i,j)
    !!   enddo
    !!  enddo
  endif

END_PROVIDER 


subroutine give_fock_ia_three_e_total(i,a,contrib)
 implicit none
 BEGIN_DOC
! contrib is the TOTAL (same spins / opposite spins) contribution from the three body term to the Fock operator 
!
 END_DOC
 integer, intent(in) :: i,a
 double precision, intent(out) :: contrib
 double precision :: int_1, int_2, int_3
 double precision :: mos_i, mos_a, w_ia
 double precision :: mos_ia, weight

 integer :: mm, ipoint,k,l

 int_1 = 0.d0
 int_2 = 0.d0
 int_3 = 0.d0
 do mm = 1, 3
  do ipoint = 1, n_points_final_grid
   weight = final_weight_at_r_vector(ipoint)                                                                          
   mos_i  = mos_in_r_array_transp(ipoint,i) 
   mos_a  = mos_in_r_array_transp(ipoint,a) 
   mos_ia = mos_a * mos_i
   w_ia   = x_W_ij_erf_rk(ipoint,mm,i,a) 
     
   int_1  += weight * fock_3_w_kk_sum(ipoint,mm) * (4.d0 * fock_3_rho_beta(ipoint) * w_ia               & 
                                                  + 2.0d0 * mos_ia * fock_3_w_kk_sum(ipoint,mm)         & 
                                                  - 2.0d0 * fock_3_w_ki_mos_k(ipoint,mm,i) * mos_a      & 
                                                  - 2.0d0 * fock_3_w_ki_mos_k(ipoint,mm,a) * mos_i      )
   int_2  += weight * (-1.d0) * ( 2.0d0 * fock_3_w_kl_mo_k_mo_l(ipoint,mm) * w_ia                     & 
                                + 2.0d0 * fock_3_rho_beta(ipoint) * fock_3_w_ki_wk_a(ipoint,mm,i,a)   & 
                                + 1.0d0 * mos_ia * fock_3_trace_w_tilde(ipoint,mm)                    )

   int_3  += weight *   1.d0  * (fock_3_w_kl_wla_phi_k(ipoint,mm,i) * mos_a + fock_3_w_kl_wla_phi_k(ipoint,mm,a) * mos_i & 
                                +fock_3_w_ki_mos_k(ipoint,mm,i)     * fock_3_w_ki_mos_k(ipoint,mm,a)                     )
  enddo
 enddo
 contrib = int_1 + int_2 + int_3

end

! ---

BEGIN_PROVIDER [double precision, diag_three_elem_hf]

  implicit none
  integer                       :: i, j, k, ipoint, mm
  double precision              :: contrib, weight, four_third, one_third, two_third, exchange_int_231
  double precision              :: integral_aaa, hthree, integral_aab, integral_abb, integral_bbb
  double precision, allocatable :: tmp(:)
  double precision, allocatable :: tmp_L(:,:), tmp_R(:,:)
  double precision, allocatable :: tmp_M(:,:), tmp_S(:), tmp_O(:), tmp_J(:,:)
  double precision, allocatable :: tmp_M_priv(:,:), tmp_S_priv(:), tmp_O_priv(:), tmp_J_priv(:,:)

  PROVIDE mo_l_coef mo_r_coef

  !print *, ' providing diag_three_elem_hf'

  if(.not. three_body_h_tc) then

   if(noL_standard) then
      PROVIDE noL_0e
      diag_three_elem_hf = noL_0e
    else
      diag_three_elem_hf = 0.d0
    endif

  else

    if(.not. bi_ortho) then

      ! ---

      one_third  = 1.d0/3.d0
      two_third  = 2.d0/3.d0
      four_third = 4.d0/3.d0
      diag_three_elem_hf = 0.d0
      do i = 1, elec_beta_num
        do j = 1, elec_beta_num
          do k = 1, elec_beta_num
            call give_integrals_3_body(k, j, i, j, i, k, exchange_int_231)   
            diag_three_elem_hf += two_third * exchange_int_231
          enddo
        enddo
      enddo
      do mm = 1, 3
        do ipoint = 1, n_points_final_grid
          weight  = final_weight_at_r_vector(ipoint)                                                                          
          contrib = 3.d0 * fock_3_w_kk_sum(ipoint,mm) * fock_3_rho_beta(ipoint) * fock_3_w_kk_sum(ipoint,mm) & 
                  - 2.d0 * fock_3_w_kl_mo_k_mo_l(ipoint,mm) * fock_3_w_kk_sum(ipoint,mm)                     & 
                  - 1.d0 * fock_3_rho_beta(ipoint) * fock_3_w_kl_w_kl(ipoint,mm)
          contrib *= four_third
          contrib += -two_third  * fock_3_rho_beta(ipoint)    * fock_3_w_kl_w_kl(ipoint,mm) & 
                     -four_third * fock_3_w_kk_sum(ipoint,mm) * fock_3_w_kl_mo_k_mo_l(ipoint,mm)
          diag_three_elem_hf += weight * contrib
       enddo
      enddo

      diag_three_elem_hf = - diag_three_elem_hf

      ! ---

    else

      ! ------------
      ! SLOW VERSION
      ! ------------

      !call give_aaa_contrib(integral_aaa)
      !call give_aab_contrib(integral_aab)
      !call give_abb_contrib(integral_abb)
      !call give_bbb_contrib(integral_bbb)
      !diag_three_elem_hf = integral_aaa + integral_aab + integral_abb + integral_bbb

      ! ------------
      ! ------------

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
    
      else
    
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
    
      endif


    endif

  endif

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, fock_3_mat_a_op_sh, (mo_num, mo_num)]
 implicit none 
 integer :: h,p,i,j
 double precision :: direct_int, exch_int, exchange_int_231, exchange_int_312
 double precision :: exchange_int_23, exchange_int_12, exchange_int_13 

 fock_3_mat_a_op_sh = 0.d0
 do h = 1, mo_num
  do p = 1, mo_num
   !F_a^{ab}(h,p) 
   do i = 1, elec_beta_num ! beta 
    do j = elec_beta_num+1, elec_alpha_num ! alpha
     call  give_integrals_3_body(h,j,i,p,j,i,direct_int)    ! <hji|pji>
     call  give_integrals_3_body(h,j,i,j,p,i,exch_int)   
     fock_3_mat_a_op_sh(h,p) -= direct_int - exch_int
    enddo
   enddo
   !F_a^{aa}(h,p)
   do i = 1, elec_beta_num ! alpha 
    do j = elec_beta_num+1, elec_alpha_num ! alpha
       call  give_integrals_3_body(h,j,i,p,j,i,direct_int) 
       call  give_integrals_3_body(h,j,i,i,p,j,exchange_int_231)
       call  give_integrals_3_body(h,j,i,j,i,p,exchange_int_312) 
       call  give_integrals_3_body(h,j,i,p,i,j,exchange_int_23) 
       call  give_integrals_3_body(h,j,i,i,j,p,exchange_int_12)
       call  give_integrals_3_body(h,j,i,j,p,i,exchange_int_13)  
       fock_3_mat_a_op_sh(h,p) -= ( direct_int + exchange_int_231 + exchange_int_312 & 
              -  exchange_int_23 & ! i <-> j
              -  exchange_int_12 & ! p <-> j
              -  exchange_int_13  )! p <-> i
    enddo 
   enddo
  enddo
 enddo
! symmetrized 
! do p = 1, elec_beta_num
!  do h = elec_alpha_num +1, mo_num
!   fock_3_mat_a_op_sh(h,p) = fock_3_mat_a_op_sh(p,h)
!  enddo
! enddo
 
! do h = elec_beta_num+1, elec_alpha_num
!  do p = elec_alpha_num +1, mo_num
!   !F_a^{bb}(h,p) 
!   do i = 1, elec_beta_num
!    do j = i+1, elec_beta_num
!     call  give_integrals_3_body(h,j,i,p,j,i,direct_int)   
!     call  give_integrals_3_body(h,j,i,p,i,j,exch_int)   
!     fock_3_mat_a_op_sh(h,p) -= direct_int - exch_int
!    enddo
!   enddo
!  enddo
! enddo

END_PROVIDER 

BEGIN_PROVIDER [ double precision, fock_3_mat_b_op_sh, (mo_num, mo_num)]
 implicit none 
 integer :: h,p,i,j
 double precision :: direct_int, exch_int
 fock_3_mat_b_op_sh = 0.d0
 do h = 1, elec_beta_num
  do p = elec_alpha_num +1, mo_num
   !F_b^{aa}(h,p) 
   do i = 1, elec_beta_num
    do j = elec_beta_num+1, elec_alpha_num
     call  give_integrals_3_body(h,j,i,p,j,i,direct_int)   
     call  give_integrals_3_body(h,j,i,p,i,j,exch_int)   
     fock_3_mat_b_op_sh(h,p) += direct_int - exch_int
    enddo
   enddo

   !F_b^{ab}(h,p) 
   do i = elec_beta_num+1, elec_beta_num
    do j = 1, elec_beta_num
     call  give_integrals_3_body(h,j,i,p,j,i,direct_int)   
     call  give_integrals_3_body(h,j,i,j,p,i,exch_int)   
     fock_3_mat_b_op_sh(h,p) += direct_int - exch_int
    enddo
   enddo
 
  enddo
 enddo

END_PROVIDER 


BEGIN_PROVIDER [ double precision, fock_3_w_kk_sum, (n_points_final_grid,3)]
 implicit none
 integer :: mm, ipoint,k
 double precision :: w_kk
 fock_3_w_kk_sum = 0.d0
 do k = 1, elec_beta_num
  do mm = 1, 3
   do ipoint = 1, n_points_final_grid
    w_kk   = x_W_ij_erf_rk(ipoint,mm,k,k) 
    fock_3_w_kk_sum(ipoint,mm) += w_kk
   enddo
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, fock_3_w_ki_mos_k, (n_points_final_grid,3,mo_num)]
 implicit none
 integer :: mm, ipoint,k,i
 double precision :: w_ki, mo_k
 fock_3_w_ki_mos_k = 0.d0
 do i = 1, mo_num
  do k = 1, elec_beta_num
   do mm = 1, 3
    do ipoint = 1, n_points_final_grid
     w_ki   = x_W_ij_erf_rk(ipoint,mm,k,i) 
     mo_k = mos_in_r_array(k,ipoint)
     fock_3_w_ki_mos_k(ipoint,mm,i) += w_ki * mo_k
    enddo
   enddo
  enddo
 enddo

END_PROVIDER 

BEGIN_PROVIDER [ double precision, fock_3_w_kl_w_kl, (n_points_final_grid,3)]
 implicit none
 integer :: k,j,ipoint,mm
 double precision :: w_kj
 fock_3_w_kl_w_kl = 0.d0
 do j = 1, elec_beta_num
  do k = 1, elec_beta_num
   do mm = 1, 3
    do ipoint = 1, n_points_final_grid
     w_kj   = x_W_ij_erf_rk(ipoint,mm,k,j) 
     fock_3_w_kl_w_kl(ipoint,mm) += w_kj * w_kj
    enddo
   enddo
  enddo
 enddo


END_PROVIDER 

BEGIN_PROVIDER [ double precision, fock_3_rho_beta, (n_points_final_grid)]
 implicit none
 integer :: ipoint,k
 fock_3_rho_beta = 0.d0
 do ipoint = 1, n_points_final_grid
  do k = 1, elec_beta_num
   fock_3_rho_beta(ipoint) += mos_in_r_array(k,ipoint) * mos_in_r_array(k,ipoint)
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, fock_3_w_kl_mo_k_mo_l, (n_points_final_grid,3)]
 implicit none
 integer :: ipoint,k,l,mm
 double precision :: mos_k, mos_l, w_kl
 fock_3_w_kl_mo_k_mo_l = 0.d0
 do k = 1, elec_beta_num
  do l = 1, elec_beta_num
   do mm = 1, 3
    do ipoint = 1, n_points_final_grid
     mos_k  = mos_in_r_array_transp(ipoint,k) 
     mos_l  = mos_in_r_array_transp(ipoint,l) 
     w_kl   = x_W_ij_erf_rk(ipoint,mm,l,k)
     fock_3_w_kl_mo_k_mo_l(ipoint,mm) += w_kl * mos_k * mos_l 
    enddo
   enddo
  enddo
 enddo

END_PROVIDER 

BEGIN_PROVIDER [ double precision, fock_3_w_ki_wk_a, (n_points_final_grid,3,mo_num, mo_num)]
 implicit none
 integer :: ipoint,i,a,k,mm
 double precision :: w_ki,w_ka
 fock_3_w_ki_wk_a = 0.d0
 do i = 1, mo_num
  do a = 1, mo_num
   do mm = 1, 3
    do ipoint = 1, n_points_final_grid
     do k = 1, elec_beta_num
      w_ki   = x_W_ij_erf_rk(ipoint,mm,k,i)
      w_ka   = x_W_ij_erf_rk(ipoint,mm,k,a)
      fock_3_w_ki_wk_a(ipoint,mm,a,i) += w_ki * w_ka
     enddo
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, fock_3_trace_w_tilde, (n_points_final_grid,3)]
 implicit none
 integer :: ipoint,k,mm
 fock_3_trace_w_tilde = 0.d0
 do k = 1, elec_beta_num
   do mm = 1, 3
    do ipoint = 1, n_points_final_grid
     fock_3_trace_w_tilde(ipoint,mm) += fock_3_w_ki_wk_a(ipoint,mm,k,k)
    enddo
   enddo
 enddo

END_PROVIDER 

BEGIN_PROVIDER [ double precision, fock_3_w_kl_wla_phi_k, (n_points_final_grid,3,mo_num)]
 implicit none
 integer :: ipoint,a,k,mm,l
 double precision :: w_kl,w_la, mo_k
 fock_3_w_kl_wla_phi_k = 0.d0
 do a = 1, mo_num
  do k = 1, elec_beta_num 
   do l = 1, elec_beta_num
    do mm = 1, 3
     do ipoint = 1, n_points_final_grid
      w_kl   = x_W_ij_erf_rk(ipoint,mm,l,k)
      w_la   = x_W_ij_erf_rk(ipoint,mm,l,a)
      mo_k  = mos_in_r_array_transp(ipoint,k) 
      fock_3_w_kl_wla_phi_k(ipoint,mm,a) += w_kl * w_la * mo_k
     enddo
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 





