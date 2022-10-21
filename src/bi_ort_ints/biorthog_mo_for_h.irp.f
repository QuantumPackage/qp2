
! ---

double precision function bi_ortho_mo_coul_ints(l, k, j, i)

  BEGIN_DOC
  !
  ! < mo^L_k mo^L_l | 1/r12 | mo^R_i mo^R_j >
  !
  END_DOC

  implicit none
  integer, intent(in) :: i, j, k, l
  integer             :: m, n, p, q

  bi_ortho_mo_coul_ints = 0.d0
  do m = 1, ao_num
    do p = 1, ao_num
      do n = 1, ao_num
        do q = 1, ao_num
          !                                   p1h1p2h2   l1                  l2              r1               r2
          bi_ortho_mo_coul_ints += ao_two_e_coul(n,q,m,p) * mo_l_coef(m,l) * mo_l_coef(n,k) * mo_r_coef(p,j) * mo_r_coef(q,i)
        enddo
      enddo
    enddo
  enddo

end function bi_ortho_mo_coul_ints

! ---

! TODO :: transform into DEGEMM

BEGIN_PROVIDER [double precision, mo_bi_ortho_coul_e_chemist, (mo_num, mo_num, mo_num, mo_num)]

  BEGIN_DOC
  !
  ! mo_bi_ortho_coul_e_chemist(k,i,l,j) = < k l | 1/r12 | i j > where i,j are right MOs and k,l are left MOs
  !
  END_DOC

  implicit none
  integer                       :: i, j, k, l, m, n, p, q
  double precision, allocatable :: mo_tmp_1(:,:,:,:), mo_tmp_2(:,:,:,:)

  allocate(mo_tmp_1(mo_num,ao_num,ao_num,ao_num))
  mo_tmp_1 = 0.d0

  do m = 1, ao_num
    do p = 1, ao_num
      do n = 1, ao_num
        do q = 1, ao_num
          do k = 1, mo_num
            !       (k n|p m)    = sum_q c_qk * (q n|p m)
            mo_tmp_1(k,n,p,m) += mo_l_coef_transp(k,q) * ao_two_e_coul(q,n,p,m)
          enddo
        enddo
      enddo
    enddo
  enddo

  allocate(mo_tmp_2(mo_num,mo_num,ao_num,ao_num))
  mo_tmp_2 = 0.d0

  do m = 1, ao_num
    do p = 1, ao_num
      do n = 1, ao_num
        do i = 1, mo_num
          do k = 1, mo_num
            !       (k i|p m) = sum_n c_ni * (k n|p m)
            mo_tmp_2(k,i,p,m) += mo_r_coef_transp(i,n) * mo_tmp_1(k,n,p,m)
          enddo
        enddo
      enddo
    enddo
  enddo
  deallocate(mo_tmp_1)

  allocate(mo_tmp_1(mo_num,mo_num,mo_num,ao_num))
  mo_tmp_1 = 0.d0
  do m = 1, ao_num
    do p = 1, ao_num
      do l = 1, mo_num
        do i = 1, mo_num
          do k = 1, mo_num
            mo_tmp_1(k,i,l,m) += mo_l_coef_transp(l,p) * mo_tmp_2(k,i,p,m)
          enddo
        enddo
      enddo
    enddo
  enddo
  deallocate(mo_tmp_2)

  mo_bi_ortho_coul_e_chemist = 0.d0 
  do m = 1, ao_num
    do j = 1, mo_num
      do l = 1, mo_num
        do i = 1, mo_num
          do k = 1, mo_num
            mo_bi_ortho_coul_e_chemist(k,i,l,j) += mo_r_coef_transp(j,m) * mo_tmp_1(k,i,l,m)
          enddo
        enddo
      enddo
    enddo
  enddo
  deallocate(mo_tmp_1)

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, mo_bi_ortho_coul_e, (mo_num, mo_num, mo_num, mo_num)]

  BEGIN_DOC
  !
  ! mo_bi_ortho_coul_e(k,l,i,j) = < k l | 1/r12 | i j > where i,j are right MOs and k,l are left MOs
  !
  END_DOC

  implicit none
  integer :: i, j, k, l

  do j = 1, mo_num
    do i = 1, mo_num
      do l = 1, mo_num
        do k = 1, mo_num
           !    < k l | V12 | i j >                  (k i|l j)
           mo_bi_ortho_coul_e(k,l,i,j) = mo_bi_ortho_coul_e_chemist(k,i,l,j)
        enddo
      enddo
    enddo
  enddo

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, mo_bi_ortho_one_e, (mo_num, mo_num)]

  BEGIN_DOC 
  !
  ! mo_bi_ortho_one_e(k,i) = < MO^L_k | h_c | MO^R_i >
  !
  END_DOC

  implicit none

  call ao_to_mo_bi_ortho(ao_one_e_integrals, ao_num, mo_bi_ortho_one_e , mo_num)

END_PROVIDER 

! ---

