
! ---

double precision function bi_ortho_mo_ints(l, k, j, i)

  BEGIN_DOC
  !
  ! <mo^L_k mo^L_l | V^TC(r_12) | mo^R_i mo^R_j>
  !
  ! WARNING :: very naive, super slow, only used to DEBUG.
  !
  END_DOC

  implicit none
  integer, intent(in) :: i, j, k, l
  integer             :: m, n, p, q

  bi_ortho_mo_ints = 0.d0
  do p = 1, ao_num
    do m = 1, ao_num
      do q = 1, ao_num
        do n = 1, ao_num
          !                                   p1h1p2h2   l1                  l2              r1               r2
          bi_ortho_mo_ints += ao_two_e_tc_tot(n,q,m,p) * mo_l_coef(m,l) * mo_l_coef(n,k) * mo_r_coef(p,j) * mo_r_coef(q,i)
        enddo
      enddo
    enddo
  enddo

end

! ---

BEGIN_PROVIDER [double precision, mo_bi_ortho_tc_two_e_chemist, (mo_num, mo_num, mo_num, mo_num)]

  BEGIN_DOC
  !
  ! mo_bi_ortho_tc_two_e_chemist(k,i,l,j) = <k l|V(r_12)|i j> where i,j are right MOs and k,l are left MOs
  !
  END_DOC

  implicit none
  integer                       :: i, j, k, l, m, n, p, q, s, r
  double precision              :: t1, t2, tt1, tt2
  double precision, allocatable :: a1(:,:,:,:), a2(:,:,:,:)
  double precision, allocatable :: a_jkp(:,:,:), a_kpq(:,:,:), ao_two_e_tc_tot_tmp(:,:,:)

  print *, ' PROVIDING mo_bi_ortho_tc_two_e_chemist ...'
  call wall_time(t1)
  call print_memory_usage()

  PROVIDE mo_r_coef mo_l_coef

  if(ao_to_mo_tc_n3) then

    print*, ' memory scale of TC ao -> mo: O(N3) '

    if(.not.read_tc_integ) then
       stop 'read_tc_integ needs to be set to true'
    endif

    allocate(a_jkp(ao_num,ao_num,mo_num))
    allocate(a_kpq(ao_num,mo_num,mo_num))
    allocate(ao_two_e_tc_tot_tmp(ao_num,ao_num,ao_num))

    open(unit=11, form="unformatted", file=trim(ezfio_filename)//'/work/ao_two_e_tc_tot', action="read")

    call wall_time(tt1)

    mo_bi_ortho_tc_two_e_chemist(:,:,:,:) = 0.d0
    do l = 1, ao_num
      read(11) ao_two_e_tc_tot_tmp(:,:,:)

      do s = 1, mo_num

        call dgemm( 'T', 'N', ao_num*ao_num, mo_num, ao_num, 1.d0              &
                  , ao_two_e_tc_tot_tmp(1,1,1), ao_num, mo_l_coef(1,1), ao_num &
                  , 0.d0, a_jkp(1,1,1), ao_num*ao_num)

        call dgemm( 'T', 'N', ao_num*mo_num, mo_num, ao_num, 1.d0 &
                  , a_jkp(1,1,1), ao_num, mo_r_coef(1,1), ao_num  &
                  , 0.d0, a_kpq(1,1,1), ao_num*mo_num)

        call dgemm( 'T', 'N', mo_num*mo_num, mo_num, ao_num, mo_r_coef(l,s) &
                  , a_kpq(1,1,1), ao_num, mo_l_coef(1,1), ao_num            &
                  , 1.d0, mo_bi_ortho_tc_two_e_chemist(1,1,1,s), mo_num*mo_num)

      enddo ! s

      if(l == 2) then
        call wall_time(tt2)
        print*, ' 1 / mo_num done in (min)', (tt2-tt1)/60.d0
        print*, ' estimated time required (min)', dble(mo_num-1)*(tt2-tt1)/60.d0
      elseif(l == 11) then
        call wall_time(tt2)
        print*, ' 10 / mo_num done in (min)', (tt2-tt1)/60.d0
        print*, ' estimated time required (min)', dble(mo_num-10)*(tt2-tt1)/(60.d0*10.d0)
      elseif(l == 101) then
        call wall_time(tt2)
        print*, ' 100 / mo_num done in (min)', (tt2-tt1)/60.d0
        print*, ' estimated time required (min)', dble(mo_num-100)*(tt2-tt1)/(60.d0*100.d0)
      endif
    enddo ! l

    close(11)

    deallocate(a_jkp, a_kpq, ao_two_e_tc_tot_tmp)

  else

    print*, ' memory scale of TC ao -> mo: O(N4) '

    allocate(a2(ao_num,ao_num,ao_num,mo_num))
  
    call dgemm( 'T', 'N', ao_num*ao_num*ao_num, mo_num, ao_num, 1.d0     &
              , ao_two_e_tc_tot(1,1,1,1), ao_num, mo_l_coef(1,1), ao_num &
              , 0.d0, a2(1,1,1,1), ao_num*ao_num*ao_num)
  
    FREE ao_two_e_tc_tot

    allocate(a1(ao_num,ao_num,mo_num,mo_num))
  
    call dgemm( 'T', 'N', ao_num*ao_num*mo_num, mo_num, ao_num, 1.d0 &
              , a2(1,1,1,1), ao_num, mo_r_coef(1,1), ao_num          &
              , 0.d0, a1(1,1,1,1), ao_num*ao_num*mo_num)
  
    deallocate(a2)
    allocate(a2(ao_num,mo_num,mo_num,mo_num))
  
    call dgemm( 'T', 'N', ao_num*mo_num*mo_num, mo_num, ao_num, 1.d0 &
              , a1(1,1,1,1), ao_num, mo_l_coef(1,1), ao_num          &
              , 0.d0, a2(1,1,1,1), ao_num*mo_num*mo_num)
  
    deallocate(a1)
  
    call dgemm( 'T', 'N', mo_num*mo_num*mo_num, mo_num, ao_num, 1.d0 &
              , a2(1,1,1,1), ao_num, mo_r_coef(1,1), ao_num          &
              , 0.d0, mo_bi_ortho_tc_two_e_chemist(1,1,1,1), mo_num*mo_num*mo_num)
  
    deallocate(a2)
  
  endif

  !allocate(a1(mo_num,ao_num,ao_num,ao_num))
  !a1 = 0.d0

  !do m = 1, ao_num
  !  do p = 1, ao_num
  !    do n = 1, ao_num
  !      do q = 1, ao_num
  !        do k = 1, mo_num
  !          !       (k n|p m)    = sum_q c_qk * (q n|p m)
  !          a1(k,n,p,m) += mo_l_coef_transp(k,q) * ao_two_e_tc_tot(q,n,p,m)
  !        enddo
  !      enddo
  !    enddo
  !  enddo
  !enddo

  !allocate(a2(mo_num,mo_num,ao_num,ao_num))
  !a2 = 0.d0

  !do m = 1, ao_num
  !  do p = 1, ao_num
  !    do n = 1, ao_num
  !      do i = 1, mo_num
  !        do k = 1, mo_num
  !          !       (k i|p m) = sum_n c_ni * (k n|p m)
  !          a2(k,i,p,m) += mo_r_coef_transp(i,n) * a1(k,n,p,m)
  !        enddo
  !      enddo
  !    enddo
  !  enddo
  !enddo
  !deallocate(a1)

  !allocate(a1(mo_num,mo_num,mo_num,ao_num))
  !a1 = 0.d0
  !do m = 1, ao_num
  !  do p = 1, ao_num
  !    do l = 1, mo_num
  !      do i = 1, mo_num
  !        do k = 1, mo_num
  !          a1(k,i,l,m) += mo_l_coef_transp(l,p) * a2(k,i,p,m)
  !        enddo
  !      enddo
  !    enddo
  !  enddo
  !enddo
  !deallocate(a2)

  !mo_bi_ortho_tc_two_e_chemist = 0.d0 
  !do m = 1, ao_num
  !  do j = 1, mo_num
  !    do l = 1, mo_num
  !      do i = 1, mo_num
  !        do k = 1, mo_num
  !          mo_bi_ortho_tc_two_e_chemist(k,i,l,j) += mo_r_coef_transp(j,m) * a1(k,i,l,m)
  !        enddo
  !      enddo
  !    enddo
  !  enddo
  !enddo
  !deallocate(a1)

  call wall_time(t2)
  print *, ' WALL TIME for PROVIDING mo_bi_ortho_tc_two_e_chemist (min)', (t2-t1)/60.d0
  call print_memory_usage()

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, mo_bi_ortho_tc_two_e, (mo_num, mo_num, mo_num, mo_num)]

  BEGIN_DOC
  !
  ! mo_bi_ortho_tc_two_e(k,l,i,j) = <k l| V(r_12) |i j> where i,j are right MOs and k,l are left MOs
  !
  ! the potential V(r_12) contains ALL TWO-E CONTRIBUTION OF THE TC-HAMILTONIAN
  !
  END_DOC

  implicit none
  integer :: i, j, k, l

  PROVIDE mo_bi_ortho_tc_two_e_chemist

  do j = 1, mo_num
    do i = 1, mo_num
      do l = 1, mo_num
        do k = 1, mo_num
           !              < k l | V12 | i j >                          (k i|l j)
           mo_bi_ortho_tc_two_e(k,l,i,j) = mo_bi_ortho_tc_two_e_chemist(k,i,l,j)
        enddo
      enddo
    enddo
  enddo

  FREE mo_bi_ortho_tc_two_e_chemist

  if(noL_standard) then
    PROVIDE noL_2e
    ! x 2 because of the Slater-Condon rules convention
    mo_bi_ortho_tc_two_e = mo_bi_ortho_tc_two_e + 2.d0 * noL_2e
    FREE noL_2e
  endif

END_PROVIDER 

! ---

 BEGIN_PROVIDER [ double precision, mo_bi_ortho_tc_two_e_jj,          (mo_num,mo_num)]
&BEGIN_PROVIDER [ double precision, mo_bi_ortho_tc_two_e_jj_exchange, (mo_num,mo_num)]
&BEGIN_PROVIDER [ double precision, mo_bi_ortho_tc_two_e_jj_anti,     (mo_num,mo_num)]

  BEGIN_DOC
  !
  ! mo_bi_ortho_tc_two_e_jj         (i,j) = J_ij = <ji|W-K|ji>
  ! mo_bi_ortho_tc_two_e_jj_exchange(i,j) = K_ij = <ij|W-K|ji>
  ! mo_bi_ortho_tc_two_e_jj_anti    (i,j) = J_ij - K_ij
  !
  END_DOC

  implicit none
  integer :: i, j

  mo_bi_ortho_tc_two_e_jj          = 0.d0
  mo_bi_ortho_tc_two_e_jj_exchange = 0.d0

  do i = 1, mo_num
    do j = 1, mo_num
      mo_bi_ortho_tc_two_e_jj         (i,j) = mo_bi_ortho_tc_two_e(j,i,j,i)
      mo_bi_ortho_tc_two_e_jj_exchange(i,j) = mo_bi_ortho_tc_two_e(i,j,j,i)
      mo_bi_ortho_tc_two_e_jj_anti    (i,j) = mo_bi_ortho_tc_two_e_jj(i,j) - mo_bi_ortho_tc_two_e_jj_exchange(i,j)
    enddo
  enddo

END_PROVIDER

! ---

 BEGIN_PROVIDER [double precision, tc_2e_3idx_coulomb_integrals , (mo_num,mo_num,mo_num)]
&BEGIN_PROVIDER [double precision, tc_2e_3idx_exchange_integrals, (mo_num,mo_num,mo_num)]

  BEGIN_DOC
  ! tc_2e_3idx_coulomb_integrals (j,k,i) = <jk|ji> 
  ! tc_2e_3idx_exchange_integrals(j,k,i) = <kj|ji> 
  END_DOC

  implicit none
  integer :: i, j, k

  do i = 1, mo_num
    do k = 1, mo_num
      do j = 1, mo_num
        tc_2e_3idx_coulomb_integrals(j, k,i) = mo_bi_ortho_tc_two_e(j ,k ,j ,i ) 
        tc_2e_3idx_exchange_integrals(j,k,i) = mo_bi_ortho_tc_two_e(k ,j ,j ,i ) 
      enddo
    enddo
  enddo

END_PROVIDER

! ---

