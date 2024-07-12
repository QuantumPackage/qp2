
BEGIN_PROVIDER [double precision, cholesky_no_1_idx_transp, (cholesky_mo_num, n_act_orb, mo_num)]
 BEGIN_DOC
 ! Cholesky vectors with ONE orbital on the active natural orbital basis 
 END_DOC
 implicit none
 integer :: i_chol,i_act,i_mo,jj_act
 double precision, allocatable :: chol_tmp(:,:)
 double precision :: wall0,wall1
 call wall_time(wall0)
 print*,'Providing cholesky_no_1_idx_transp'
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
  call dgemm('N','N',cholesky_mo_num,n_act_orb,n_act_orb,1.d0,  &
        chol_tmp, size(chol_tmp,1),                              &
        natorbsCI, size(natorbsCI,1),                                              &
        0.d0,                                                      &
        cholesky_no_1_idx_transp(1,1,i_mo), size(cholesky_no_1_idx_transp,1))
 enddo
 call wall_time(wall1)
 print*,'Time to provide cholesky_no_1_idx_transp = ', wall1 - wall0

END_PROVIDER 

BEGIN_PROVIDER [double precision, cholesky_no_2_idx_transp, (cholesky_mo_num, n_act_orb, n_act_orb)]
 BEGIN_DOC
 ! Cholesky vectors with TWO orbital on the active natural orbital basis 
 END_DOC
 implicit none
 integer :: i_chol,i_act,j_act,jj_act
 double precision, allocatable :: chol_tmp(:,:),chol_tmp_bis(:,:)
 allocate(chol_tmp(cholesky_mo_num,n_act_orb),chol_tmp_bis(cholesky_mo_num,n_act_orb))
 double precision :: wall0,wall1
 call wall_time(wall0)
 print*,'Providing cholesky_no_2_idx_transp'
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
 call wall_time(wall1)
 print*,'Time to provide  cholesky_no_2_idx_transp = ', wall1 - wall0

END_PROVIDER 

BEGIN_PROVIDER [ double precision, cholesky_no_total_transp, (cholesky_mo_num, mo_num, mo_num)]
 implicit none
 BEGIN_DOC
 ! Cholesky vectors defined on all basis including the NO basis
 END_DOC
 integer :: i_chol, i_act, ii_act, j_act, jj_act, i_core_inact, j_core_inact, ii_core_inact, jj_core_inact
 integer :: i_virt, ii_virt, j_virt, jj_virt
 double precision :: wall0,wall1
 call wall_time(wall0)
 print*,'Providing cholesky_no_total_transp '
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

 call wall_time(wall1)
 print*,'Time to provide cholesky_no_total_transp = ', wall1 - wall0
END_PROVIDER 


double precision function bielec_no_basis(i_1,j_1,i_2,j_2)
 implicit none
 integer, intent(in) :: i_1,j_1,i_2,j_2
  BEGIN_DOC
  ! integral (i_1 j_1|i_2 j_2) in the mixed basis of both MOs and natural MOs
  ! 
  END_DOC
  integer :: i 
 bielec_no_basis = 0.d0
 do i = 1, cholesky_mo_num
  bielec_no_basis += cholesky_no_total_transp(i,i_1, j_1) * cholesky_no_total_transp(i,i_2,j_2)
 enddo
end

double precision function bielec_PQxx_no(i_mo, j_mo, i_ca, j_ca)
 implicit none
 BEGIN_DOC
 ! function that computes (i_mo j_mo| i_ca j_ca) with Cholesky decomposition  on the NO basis for active orbitals 
 ! 
 ! where i_ca, j_ca are in [1:n_core_inact_act_orb]
 END_DOC
 integer, intent(in) :: i_ca, j_ca, i_mo, j_mo
 integer :: ii_ca, jj_ca
 double precision :: bielec_no_basis
 ii_ca = list_core_inact_act(i_ca)
 jj_ca = list_core_inact_act(j_ca)
 bielec_PQxx_no = bielec_no_basis(i_mo,j_mo,ii_ca,jj_ca)
end

double precision function bielec_PxxQ_no(i_mo, j_ca, i_ca, j_mo)
 implicit none 
  BEGIN_DOC
  ! function that computes (i_mo j_ca |i_ca j_mo) with Cholesky decomposition  on the NO basis for active orbitals 
  ! 
  ! where i_ca, j_ca are in [1:n_core_inact_act_orb]
  END_DOC
 integer, intent(in) :: i_ca, j_ca, i_mo, j_mo
 integer :: ii_ca, jj_ca
 double precision :: bielec_no_basis
 ii_ca = list_core_inact_act(i_ca)
 jj_ca = list_core_inact_act(j_ca)
 bielec_PxxQ_no = bielec_no_basis(i_mo, jj_ca, ii_ca, j_mo)

end


double precision function bielec_PQxx(i_mo, j_mo, i_ca, j_ca)
  BEGIN_DOC
  ! function that computes (i_mo j_mo |i_ca j_ca) with Cholesky decomposition 
  ! 
  ! indices are unshifted orbital numbers
  ! 
  ! where i_ca, j_ca are in [1:n_core_inact_act_orb]
  END_DOC
 implicit none 
 integer, intent(in) :: i_ca, j_ca, j_mo, i_mo
 double precision :: mo_two_e_integral
 integer :: ii_ca, jj_ca
 ii_ca = list_core_inact_act(i_ca)
 jj_ca = list_core_inact_act(j_ca)
 bielec_PQxx = mo_two_e_integral(i_mo,ii_ca,j_mo,jj_ca)
end

double precision function bielec_PxxQ(i_mo, i_ca, j_ca, j_mo)
  BEGIN_DOC
  ! function that computes (i_mo j_mo |i_ca j_ca) with Cholesky decomposition 
  ! 
  ! where i_ca, j_ca are in [1:n_core_inact_act_orb]
  END_DOC
 implicit none
 integer, intent(in) :: i_ca, j_ca, j_mo, i_mo
 double precision :: mo_two_e_integral
 integer :: ii_ca, jj_ca
 ii_ca = list_core_inact_act(i_ca)
 jj_ca = list_core_inact_act(j_ca)
 bielec_PxxQ = mo_two_e_integral(i_mo,jj_ca,ii_ca,j_mo)
end

