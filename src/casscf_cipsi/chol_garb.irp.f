
!!!!! FUNCTIONS THAT WORK BUT WHICH ARE USELESS AS THE ARRAYS CAN ALWAYS BE STORED
!double precision function bielecCI_chol(i_a, j_a, k_a, i_mo)
!  BEGIN_DOC
!  ! function that computes (i_a j_a |k_a j_mo) with Cholesky decomposition 
!  ! 
!  ! where i_a, j_a, k_a are in [1:n_act_orb] !!! ONLY ON ACTIVE 
!  END_DOC
! implicit none
! integer, intent(in) :: i_a, j_a, k_a, i_mo
! integer :: ii_a, jj_a, kk_a
! double precision :: mo_two_e_integral
! ii_a = list_act(i_a)
! jj_a = list_act(j_a)
! kk_a = list_act(k_a)
! bielecCI_chol = mo_two_e_integral(ii_a,kk_a,jj_a,i_mo)
!end

!double precision function bielecCI_no_chol(i_ca, j_ca, k_ca, i_mo)
!  BEGIN_DOC
!  ! function that computes (i_ca j_ca |k_ca j_mo) with Cholesky decomposition on the NO basis for active orbitals 
!  ! 
!  ! where i_ca, j_ca, k_ca are in [1:n_core_inact_act_orb]
!  END_DOC
! implicit none 
! integer, intent(in) :: i_ca, j_ca, k_ca, i_mo
! integer :: ii_ca, jj_ca, kk_ca
! double precision :: bielec_no_basis_chol
! ii_ca = list_act(i_ca)
! jj_ca = list_act(j_ca)
! kk_ca = list_act(k_ca)
! bielecCI_no_chol = bielec_no_basis_chol(ii_ca, jj_ca, kk_ca, i_mo)
!
!end
