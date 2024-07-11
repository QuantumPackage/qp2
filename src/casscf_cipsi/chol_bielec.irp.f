
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


BEGIN_PROVIDER [double precision, cholesky_no_2_idx_transp, (cholesky_mo_num, n_act_orb, n_act_orb)]
 BEGIN_DOC
 ! Cholesky vectors with TWO orbital on the active natural orbital basis 
 END_DOC
 implicit none
 integer :: i_chol,i_act,j_act,jj_act
 double precision, allocatable :: chol_tmp(:,:)
 allocate(chol_tmp(cholesky_mo_num,n_act_orb))
 cholesky_no_2_idx_transp = 0.D0
 do j_act = 1, n_act_orb
  do i_act = 1, n_act_orb
   do jj_act = 1, n_act_orb
    do i_chol = 1, cholesky_mo_num
     cholesky_no_2_idx_transp(i_chol, i_act, j_act) += cholesky_no_1_idx_transp(i_chol, i_act,jj_act) * natorbsCI(jj_act,i_act) 
    enddo
   enddo
  enddo
 enddo

END_PROVIDER 


BEGIN_PROVIDER [double precision, cholesky_no_2_idx_transp_dgemm, (cholesky_mo_num, n_act_orb, n_act_orb)]
 BEGIN_DOC
 ! Cholesky vectors with TWO orbital on the active natural orbital basis 
 END_DOC
 implicit none
 integer :: i_chol,i_act,j_act,jj_act
 double precision, allocatable :: chol_tmp(:,:)
 allocate(chol_tmp(cholesky_mo_num,n_act_orb))
 cholesky_no_2_idx_transp_dgemm = 0.D0
 do j_act = 1, n_act_orb
  ! Get all the integrals corresponding to the "j_act"
  do i_act = 1, n_act_orb
   jj_act = list_act(i_act)
   do i_chol = 1, cholesky_mo_num
    chol_tmp(i_chol, i_act) = cholesky_no_1_idx_transp(i_chol, j_act, jj_act)
   enddo
  enddo
!  ! Do the matrix product 
!  do i_act = 1, n_act_orb
!   do jj_act = 1, n_act_orb
!    do i_chol = 1, cholesky_mo_num
!     cholesky_no_1_idx_transp(i_chol, i_act, j_act) += chol_tmp(i_chol, jj_act) * natorbsCI(jj_act,i_act) 
!    enddo
!   enddo
!  enddo
  call dgemm('N','N',cholesky_mo_num,n_act_orb,n_act_orb,1.d0,  &
        chol_tmp, size(chol_tmp,1),                              &
        natorbsCI, size(natorbsCI,1),                                              &
        0.d0,                                                      &
        cholesky_no_2_idx_transp_dgemm(1,1,j_act), size(cholesky_no_2_idx_transp_dgemm,1))
 enddo

END_PROVIDER 


