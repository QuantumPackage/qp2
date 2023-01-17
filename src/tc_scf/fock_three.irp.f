BEGIN_PROVIDER [ double precision, fock_3_mat, (mo_num, mo_num)] 
 implicit none
  integer :: i,j
  double precision :: contrib
  fock_3_mat = 0.d0
  if(.not.bi_ortho.and.three_body_h_tc)then
   call give_fock_ia_three_e_total(1,1,contrib)
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
  else if(bi_ortho.and.three_body_h_tc)then
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
  integer          :: i, j, k, ipoint, mm
  double precision :: contrib, weight, four_third, one_third, two_third, exchange_int_231
  double precision :: integral_aaa, hthree, integral_aab, integral_abb, integral_bbb

  PROVIDE mo_l_coef mo_r_coef

  !print *, ' providing diag_three_elem_hf'

  if(.not. three_body_h_tc) then

    diag_three_elem_hf = 0.d0

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
            call give_integrals_3_body(k, j, i, j, i, k,exchange_int_231)   
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

      provide mo_l_coef mo_r_coef
      call give_aaa_contrib(integral_aaa)
      call give_aab_contrib(integral_aab)
      call give_abb_contrib(integral_abb)
      call give_bbb_contrib(integral_bbb)
      diag_three_elem_hf = integral_aaa + integral_aab + integral_abb + integral_bbb
!      print*,'integral_aaa + integral_aab + integral_abb + integral_bbb'
!      print*,integral_aaa , integral_aab , integral_abb , integral_bbb

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
       direct_int = three_body_4_index(j,i,h,p)                    
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
