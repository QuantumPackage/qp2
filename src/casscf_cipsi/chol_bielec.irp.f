
BEGIN_PROVIDER [double precision, cholesky_no_1_idx_transp, (cholesky_mo_num, n_act_orb, mo_num)]
 BEGIN_DOC
 ! Cholesky vectors with ONE orbital on the active natural orbital basis 
 END_DOC
 implicit none
 integer :: i_chol,i_act,i_mo,jj_act
 double precision, allocatable :: chol_tmp(:,:)
 allocate(chol_tmp(cholesky_mo_num,n_act_orb))
 cholesky_no_1_idx_transp = 0.D0
 do i_mo = 1, mo_num
  ! Get all the integrals corresponding to the "i_mo"
  do i_act = 1, n_act_orb
   jj_act = list_act(i_act)
   do i_chol = 1, cholesky_mo_num
    chol_tmp(i_chol, i_act) = cholesky_mo_transp(i_chol, jj_act, i_mo)
   enddo
  enddo
!  ! Do the matrix product 
!  do i_act = 1, n_act_orb
!   do jj_act = 1, n_act_orb
!    do i_chol = 1, cholesky_mo_num
!     cholesky_no_1_idx_transp(i_chol, i_act, i_mo) += chol_tmp(i_chol, jj_act) * natorbsCI(jj_act,i_act) 
!    enddo
!   enddo
!  enddo
  call dgemm('N','N',cholesky_mo_num,n_act_orb,n_act_orb,1.d0,  &
        chol_tmp, size(chol_tmp,1),                              &
        natorbsCI, size(natorbsCI,1),                                              &
        0.d0,                                                      &
        cholesky_no_1_idx_transp(1,1,i_mo), size(cholesky_no_1_idx_transp,1))
 enddo

END_PROVIDER 

BEGIN_PROVIDER [double precision, cholesky_no_2_idx_transp_old, (cholesky_mo_num, n_act_orb, n_act_orb)]
 BEGIN_DOC
 ! Cholesky vectors with TWO orbital on the active natural orbital basis 
 END_DOC
 implicit none
 integer :: i_chol,i_act,j_act,jj_act,jjj_act
 double precision, allocatable :: chol_tmp(:,:)
 allocate(chol_tmp(cholesky_mo_num,n_act_orb))
 cholesky_no_2_idx_transp_old = 0.D0
 do jj_act = 1, n_act_orb
  jjj_act = list_act(jj_act)
  do j_act = 1, n_act_orb
   do i_act = 1, n_act_orb
    do i_chol = 1, cholesky_mo_num
     cholesky_no_2_idx_transp_old(i_chol, i_act, j_act) += cholesky_no_1_idx_transp(i_chol, i_act,jjj_act) * natorbsCI(jj_act,j_act) 
    enddo
   enddo
  enddo
 enddo

END_PROVIDER 


BEGIN_PROVIDER [double precision, cholesky_no_2_idx_transp, (cholesky_mo_num, n_act_orb, n_act_orb)]
 BEGIN_DOC
 ! Cholesky vectors with TWO orbital on the active natural orbital basis 
 END_DOC
 implicit none
 integer :: i_chol,i_act,j_act,jj_act
 double precision, allocatable :: chol_tmp(:,:),chol_tmp_bis(:,:)
 allocate(chol_tmp(cholesky_mo_num,n_act_orb),chol_tmp_bis(cholesky_mo_num,n_act_orb))
 cholesky_no_2_idx_transp = 0.D0
 do i_act = 1, n_act_orb
  ! Get all the integrals corresponding to the "j_act"
  do j_act = 1, n_act_orb
   jj_act = list_act(j_act)
   do i_chol = 1, cholesky_mo_num
    chol_tmp(i_chol, j_act) = cholesky_no_1_idx_transp(i_chol, i_act, jj_act)
   enddo
  enddo
  call dgemm('N','N',cholesky_mo_num,n_act_orb,n_act_orb,1.d0,  &
        chol_tmp, size(chol_tmp,1),                              &
        natorbsCI, size(natorbsCI,1),                                              &
        0.d0,                                                      &
        cholesky_no_2_idx_transp(1,1,i_act), size(cholesky_no_2_idx_transp,1))
 enddo

END_PROVIDER 

BEGIN_PROVIDER [ double precision, cholesky_no_total_transp, (cholesky_mo_num, mo_num, mo_num)]
 implicit none
 BEGIN_DOC
 ! Cholesky vectors defined on all basis including the NO basis
 END_DOC
 integer :: i_chol, i_act, ii_act, j_act, jj_act, i_core_inact, j_core_inact, ii_core_inact, jj_core_inact
 integer :: i_virt, ii_virt, j_virt, jj_virt
 ! Block when two orbitals belong to the core/inact 
 do j_core_inact = 1, n_core_inact_orb
  jj_core_inact = list_core_inact(j_core_inact)
  do i_core_inact = 1, n_core_inact_orb
   ii_core_inact = list_core_inact(i_core_inact)
   do i_chol = 1, cholesky_mo_num
    cholesky_no_total_transp(i_chol, ii_core_inact, jj_core_inact) = cholesky_mo_transp(i_chol,ii_core_inact,jj_core_inact)
   enddo
  enddo
 enddo

 ! Block when one orbitals belongs to the core/inact and one belongs to the active
 do j_core_inact = 1, n_core_inact_orb
  jj_core_inact = list_core_inact(j_core_inact)
  do i_act = 1, n_act_orb
   ii_act = list_act(i_act)
   do i_chol = 1, cholesky_mo_num 
    cholesky_no_total_transp(i_chol,ii_act,j_core_inact) = cholesky_no_1_idx_transp(i_chol,i_act,jj_core_inact)
   enddo
  enddo
 enddo
 do j_core_inact = 1, n_core_inact_orb
  jj_core_inact = list_core_inact(j_core_inact)
  do i_act = 1, n_act_orb
   ii_act = list_act(i_act)
   do i_chol = 1, cholesky_mo_num 
    cholesky_no_total_transp(i_chol,j_core_inact,ii_act) = cholesky_no_1_idx_transp(i_chol,i_act,jj_core_inact)
   enddo
  enddo
 enddo
 
 ! Block when two orbitals belong to the active 
 do j_act = 1, n_act_orb
  jj_act = list_act(j_act)
  do i_act = 1, n_act_orb
   ii_act = list_act(i_act)
   do i_chol = 1, cholesky_mo_num 
    cholesky_no_total_transp(i_chol,ii_act,jj_act) = cholesky_no_2_idx_transp(i_chol,i_act,j_act)
   enddo
  enddo
 enddo

 ! Block when two orbitals belong to the virtuals 
 do i_virt = 1, n_virt_orb
  ii_virt = list_virt(i_virt)
  do j_virt = 1, n_virt_orb
   jj_virt = list_virt(j_virt)
   do i_chol = 1, cholesky_mo_num
    cholesky_no_total_transp(i_chol,jj_virt,ii_virt) = cholesky_mo_transp(i_chol,jj_virt,ii_virt)
   enddo
  enddo
 enddo
 
 ! Block when one orbital is in active and the other in the virtuals 
 do i_virt = 1, n_virt_orb
  ii_virt = list_virt(i_virt)
  do i_act = 1, n_act_orb
   ii_act = list_act(i_act)
   do i_chol = 1, cholesky_mo_num
    cholesky_no_total_transp(i_chol,ii_act,ii_virt) = cholesky_no_1_idx_transp(i_chol, i_act,ii_virt)
   enddo
  enddo
 enddo
 do i_virt = 1, n_virt_orb
  ii_virt = list_virt(i_virt)
  do i_act = 1, n_act_orb
   ii_act = list_act(i_act)
   do i_chol = 1, cholesky_mo_num
    cholesky_no_total_transp(i_chol,ii_virt,ii_act) = cholesky_no_1_idx_transp(i_chol, i_act,ii_virt)
   enddo
  enddo
 enddo
 ! Block when one orbital is in the virtual and one in the core-inact 
 do i_virt = 1, n_virt_orb
  ii_virt = list_virt(i_virt)
  do i_core_inact = 1, n_core_inact_orb
   ii_core_inact = list_core_inact(i_core_inact)
   do i_chol = 1, cholesky_mo_num
    cholesky_no_total_transp(i_chol, ii_core_inact, ii_virt) = cholesky_mo_transp(i_chol, ii_core_inact, ii_virt)
   enddo
  enddo
 enddo
 do i_core_inact = 1, n_core_inact_orb
  ii_core_inact = list_core_inact(i_core_inact)
  do i_virt = 1, n_virt_orb
   ii_virt = list_virt(i_virt)
   do i_chol = 1, cholesky_mo_num
    cholesky_no_total_transp(i_chol, ii_virt, ii_core_inact) = cholesky_mo_transp(i_chol, ii_virt, ii_core_inact)
   enddo
  enddo
 enddo
END_PROVIDER 


double precision function bielec_no_basis_chol(i_1,j_1,i_2,j_2)
 implicit none
 integer, intent(in) :: i_1,j_1,i_2,j_2
  BEGIN_DOC
  ! integral (i_1 j_1|i_2 j_2) in the mixed basis of both MOs and natural MOs
  ! 
  END_DOC
  integer :: i_chol 
 bielec_no_basis_chol = 0.d0
 do i_chol = 1, cholesky_mo_num
  bielec_no_basis_chol += cholesky_no_total_transp(i_chol,i_1, j_1) * cholesky_no_total_transp(i_chol,i_2,j_2)
 enddo
end

double precision function bielec_PQxx_no_chol(i_mo, j_mo, i_ca, j_ca)
 implicit none
 BEGIN_DOC
 ! function that computes (i_mo j_mo| i_ca j_ca) with Cholesky decomposition 
 ! 
 ! indices are unshifted orbital numbers
 END_DOC
 integer, intent(in) :: i_ca, j_ca, i_mo, j_mo
 integer :: ii_ca, jj_ca
 double precision :: bielec_no_basis_chol
 ii_ca = list_core_inact_act(i_ca)
 jj_ca = list_core_inact_act(j_ca)
 bielec_PQxx_no_chol = bielec_no_basis_chol(i_mo,j_mo,ii_ca,jj_ca)

end

double precision function bielec_PxxQ_no_chol(i_mo, j_ca, i_ca, j_mo)
 implicit none 
  BEGIN_DOC
  ! function that computes (i_mo j_ca |i_ca j_mo) with Cholesky decomposition 
  ! 
  ! indices are unshifted orbital numbers
  END_DOC
 integer, intent(in) :: i_ca, j_ca, i_mo, j_mo
 integer :: ii_ca, jj_ca
 double precision :: bielec_no_basis_chol
 ii_ca = list_core_inact_act(i_ca)
 jj_ca = list_core_inact_act(j_ca)
 bielec_PxxQ_no_chol = bielec_no_basis_chol(i_mo, jj_ca, ii_ca, j_mo)

end

double precision function bielecCI_no_chol(i_ca, j_ca, k_ca, i_mo)
 implicit none 
 integer, intent(in) :: i_ca, j_ca, k_ca, i_mo
 integer :: ii_ca, jj_ca, kk_ca
 double precision :: bielec_no_basis_chol
 ii_ca = list_act(i_ca)
 jj_ca = list_act(j_ca)
 kk_ca = list_act(k_ca)
 bielecCI_no_chol = bielec_no_basis_chol(ii_ca, jj_ca, kk_ca, i_mo)

end
