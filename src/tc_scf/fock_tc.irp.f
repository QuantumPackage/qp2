
! ---

 BEGIN_PROVIDER [ double precision, two_e_tc_non_hermit_integral_alpha, (ao_num, ao_num)]
&BEGIN_PROVIDER [ double precision, two_e_tc_non_hermit_integral_beta , (ao_num, ao_num)]
 BEGIN_DOC
! two_e_tc_non_hermit_integral_alpha(k,i) = <k| F^tc_alpha |i> 
!
! where F^tc is the two-body part of the TC Fock matrix and k,i are AO basis functions
 END_DOC
  implicit none
  integer          :: i, j, k, l
  double precision :: density, density_a, density_b

  two_e_tc_non_hermit_integral_alpha = 0.d0
  two_e_tc_non_hermit_integral_beta  = 0.d0

  !! TODO :: parallelization properly done
  do i = 1, ao_num
    do k = 1, ao_num
!!$OMP PARALLEL                  &
!!$OMP DEFAULT (NONE)            &
!!$OMP PRIVATE (j,l,density_a,density_b,density) & 
!!$OMP SHARED (i,k,ao_num,SCF_density_matrix_ao_alpha,SCF_density_matrix_ao_beta,ao_non_hermit_term_chemist) & 
!!$OMP SHARED (two_e_tc_non_hermit_integral_alpha,two_e_tc_non_hermit_integral_beta)
!!$OMP DO SCHEDULE (dynamic)
      do j = 1, ao_num
        do l = 1, ao_num

          density_a = TCSCF_density_matrix_ao_alpha(l,j)
          density_b = TCSCF_density_matrix_ao_beta (l,j)
          density   = density_a + density_b                      

          !                                         rho(l,j)   *      < k l| T | i j>
          two_e_tc_non_hermit_integral_alpha(k,i) += density   * ao_two_e_tc_tot(l,j,k,i)
          !                                         rho(l,j)   *      < k l| T | i j>
          two_e_tc_non_hermit_integral_beta (k,i) += density   * ao_two_e_tc_tot(l,j,k,i)
          !                                         rho_a(l,j) *      < l k| T | i j>
          two_e_tc_non_hermit_integral_alpha(k,i) -= density_a * ao_two_e_tc_tot(k,j,l,i)
          !                                         rho_b(l,j) *      < l k| T | i j>
          two_e_tc_non_hermit_integral_beta (k,i) -= density_b * ao_two_e_tc_tot(k,j,l,i)

        enddo
      enddo
!!$OMP END DO
!!$OMP END PARALLEL
    enddo
  enddo

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, Fock_matrix_tc_ao_alpha, (ao_num, ao_num)]
  implicit none
  BEGIN_DOC
 ! Total alpha TC Fock matrix : h_c + Two-e^TC terms on the AO basis
  END_DOC
  Fock_matrix_tc_ao_alpha =  ao_one_e_integrals_tc_tot &
                          + two_e_tc_non_hermit_integral_alpha 

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, Fock_matrix_tc_ao_beta, (ao_num, ao_num)]

  BEGIN_DOC
 ! Total beta TC Fock matrix : h_c + Two-e^TC terms on the AO basis
  END_DOC
  implicit none

  Fock_matrix_tc_ao_beta = ao_one_e_integrals_tc_tot &
                         + two_e_tc_non_hermit_integral_beta 

END_PROVIDER 
! ---

BEGIN_PROVIDER [ double precision, Fock_matrix_tc_ao_tot, (ao_num, ao_num) ]
  implicit none
  BEGIN_DOC
 ! Total alpha+beta TC Fock matrix : h_c + Two-e^TC terms on the AO basis
  END_DOC
  Fock_matrix_tc_ao_tot = 0.5d0 * (Fock_matrix_tc_ao_alpha + Fock_matrix_tc_ao_beta)
END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, Fock_matrix_tc_mo_alpha, (mo_num, mo_num) ]
  implicit none
  BEGIN_DOC
 ! Total alpha TC Fock matrix : h_c + Two-e^TC terms on the MO basis
  END_DOC
  if(bi_ortho)then
   call ao_to_mo_bi_ortho( Fock_matrix_tc_ao_alpha, size(Fock_matrix_tc_ao_alpha, 1) &
                         , Fock_matrix_tc_mo_alpha, size(Fock_matrix_tc_mo_alpha, 1) )
  else
   call ao_to_mo(  Fock_matrix_tc_ao_alpha, size(Fock_matrix_tc_ao_alpha, 1) &
                 , Fock_matrix_tc_mo_alpha, size(Fock_matrix_tc_mo_alpha, 1) )
  endif
END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, Fock_matrix_tc_mo_beta, (mo_num,mo_num) ]
  implicit none
  BEGIN_DOC
 ! Total beta  TC Fock matrix : h_c + Two-e^TC terms on the MO basis
  END_DOC
  if(bi_ortho)then
   call ao_to_mo_bi_ortho( Fock_matrix_tc_ao_beta, size(Fock_matrix_tc_ao_beta, 1) &
                         , Fock_matrix_tc_mo_beta, size(Fock_matrix_tc_mo_beta, 1) )
  else
   call ao_to_mo(  Fock_matrix_tc_ao_beta, size(Fock_matrix_tc_ao_beta, 1) &
                 , Fock_matrix_tc_mo_beta, size(Fock_matrix_tc_mo_beta, 1) )
  endif
END_PROVIDER


BEGIN_PROVIDER [ double precision, Fock_matrix_tc_mo_tot, (mo_num, mo_num)]
  implicit none
  BEGIN_DOC
 ! Total alpha+beta  TC Fock matrix : h_c + Two-e^TC terms on the MO basis
  END_DOC
  Fock_matrix_tc_mo_tot = 0.5d0 * (Fock_matrix_tc_mo_alpha + Fock_matrix_tc_mo_beta)
  if(three_body_h_tc) then
    Fock_matrix_tc_mo_tot += fock_3_mat
  endif
  !call restore_symmetry(mo_num, mo_num, Fock_matrix_tc_mo_tot, mo_num, 1.d-10)
END_PROVIDER 

! ---

 BEGIN_PROVIDER [ double precision, grad_non_hermit_left]
&BEGIN_PROVIDER [ double precision, grad_non_hermit_right]
&BEGIN_PROVIDER [ double precision, grad_non_hermit]
 implicit none
  integer :: i, k
  grad_non_hermit_left = 0.d0
  grad_non_hermit_right = 0.d0
  do i = 1, elec_beta_num ! doc --> SOMO
    do k = elec_beta_num+1, elec_alpha_num
      grad_non_hermit_left+= dabs(Fock_matrix_tc_mo_tot(k,i))
      grad_non_hermit_right+= dabs(Fock_matrix_tc_mo_tot(i,k))
    enddo
  enddo
  do i = 1, elec_beta_num ! doc --> virt 
    do k = elec_alpha_num+1, mo_num
      grad_non_hermit_left+= dabs(Fock_matrix_tc_mo_tot(k,i))
      grad_non_hermit_right+= dabs(Fock_matrix_tc_mo_tot(i,k))
    enddo
  enddo
  do i = elec_beta_num+1, elec_alpha_num ! SOMO --> virt 
    do k = elec_alpha_num+1, mo_num
      grad_non_hermit_left+= dabs(Fock_matrix_tc_mo_tot(k,i))
      grad_non_hermit_right+= dabs(Fock_matrix_tc_mo_tot(i,k))
    enddo
  enddo
 grad_non_hermit = grad_non_hermit_left + grad_non_hermit_right
END_PROVIDER 
