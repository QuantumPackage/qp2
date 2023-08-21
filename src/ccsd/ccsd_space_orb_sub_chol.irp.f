subroutine ccsd_energy_space_chol(nO,nV,tau,t1,energy)

  implicit none

  integer, intent(in)           :: nO, nV
  double precision, intent(in)  :: tau(nO,nO,nV,nV)
  double precision, intent(in)  :: t1(nO,nV)
  double precision, intent(out) :: energy

  ! internal
  integer :: i,j,a,b
  double precision :: e

  energy = 0d0
  !$omp parallel &
  !$omp shared(nO,nV,energy,tau,t1,&
  !$omp cc_space_f_vo,cc_space_w_oovv) &
  !$omp private(i,j,a,b,e) &
  !$omp default(none)
  e = 0d0
  !$omp do
  do a = 1, nV
    do i = 1, nO
      e = e + 2d0 * cc_space_f_vo(a,i) * t1(i,a)
    enddo
  enddo
  !$omp end do nowait
  !$omp do
  do b = 1, nV
    do a = 1, nV
      do j = 1, nO
        do i = 1, nO
          e = e + tau(i,j,a,b) * cc_space_w_oovv(i,j,a,b)
       enddo
      enddo
    enddo
  enddo
  !$omp end do nowait
  !$omp critical
  energy = energy + e
  !$omp end critical
  !$omp end parallel

end

! Tau

subroutine update_tau_space_chol(nO,nV,t1,t2,tau)

  implicit none

  ! in
  integer, intent(in)           :: nO, nV
  double precision, intent(in)  :: t1(nO,nV), t2(nO,nO,nV,nV)

  ! out
  double precision, intent(out) :: tau(nO,nO,nV,nV)

  ! internal
  integer                       :: i,j,a,b

  !$OMP PARALLEL &
  !$OMP SHARED(nO,nV,tau,t2,t1) &
  !$OMP PRIVATE(i,j,a,b) &
  !$OMP DEFAULT(NONE)
  !$OMP DO
  do b = 1, nV
    do a = 1, nV
      do j = 1, nO
        do i = 1, nO
          tau(i,j,a,b) = t2(i,j,a,b) + t1(i,a) * t1(j,b)
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

end

! R1

subroutine compute_r1_space_chol(nO,nV,t1,t2,tau,H_oo,H_vv,H_vo,r1,max_r1)

  implicit none

  ! in
  integer, intent(in)           :: nO, nV
  double precision, intent(in)  :: t1(nO,nV), t2(nO,nO,nV,nV), tau(nO,nO,nV,nV)
  double precision, intent(in)  :: H_oo(nO,nO), H_vv(nV,nV), H_vo(nV,nO)

  ! out
  double precision, intent(out) :: r1(nO,nV), max_r1

  ! internal
  integer                       :: u,i,j,beta,a,b

  !$omp parallel &
  !$omp shared(nO,nV,r1,cc_space_f_ov) &
  !$omp private(u,beta) &
  !$omp default(none)
  !$omp do
  do beta = 1, nV
    do u = 1, nO
      r1(u,beta) = cc_space_f_ov(u,beta)
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  double precision, allocatable :: X_oo(:,:)
  allocate(X_oo(nO,nO))
  call dgemm('N','N', nO, nO, nV, &
             -2d0, t1    , size(t1,1), &
                   cc_space_f_vo, size(cc_space_f_vo,1), &
              0d0, X_oo  , size(X_oo,1))

  call dgemm('T','N', nO, nV, nO, &
             1d0, X_oo, size(X_oo,2), &
                  t1  , size(t1,1), &
             1d0, r1  , size(r1,1))
  deallocate(X_oo)

  call dgemm('N','N', nO, nV, nV, &
             1d0, t1  , size(t1,1), &
                  H_vv, size(H_vv,1), &
             1d0, r1  , size(r1,1))

  call dgemm('N','N', nO, nV, nO, &
             -1d0, H_oo, size(H_oo,1), &
                   t1  , size(t1,1), &
              1d0, r1, size(r1,1))

  double precision, allocatable :: X_voov(:,:,:,:)
  allocate(X_voov(nV, nO, nO, nV))

  !$omp parallel &
  !$omp shared(nO,nV,X_voov,t2,t1) &
  !$omp private(u,beta,i,a) &
  !$omp default(none)
  !$omp do
  do beta = 1, nV
    do u = 1, nO
      do i = 1, nO
        do a = 1, nV
          X_voov(a,i,u,beta) = 2d0 * t2(i,u,a,beta) - t2(u,i,a,beta) + t1(u,a) * t1(i,beta)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  call dgemv('T', nV*nO, nO*nV, &
             1d0, X_voov, size(X_voov,1) * size(X_voov,2), &
                  H_vo  , 1, &
             1d0, r1    , 1)

  deallocate(X_voov)

  double precision, allocatable :: X_ovov(:,:,:,:)
  allocate(X_ovov(nO, nV, nO, nV))

  !$omp parallel &
  !$omp shared(nO,nV,cc_space_v_ovov,cc_space_v_voov,X_ovov) &
  !$omp private(u,beta,i,a) &
  !$omp default(none)
  !$omp do
  do beta = 1, nV
    do u = 1, nO
      do a = 1, nv
        do i = 1, nO
          X_ovov(i,a,u,beta) = 2d0 * cc_space_v_voov(a,u,i,beta) - cc_space_v_ovov(u,a,i,beta)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  call dgemv('T', nO*nV, nO*nV, &
             1d0, X_ovov, size(X_ovov,1) * size(X_ovov,2), &
                  t1     , 1, &
             1d0, r1     , 1)

  deallocate(X_ovov)

  integer :: iblock, block_size, nVmax
  double precision, allocatable :: W_vvov(:,:,:,:), W_vvov_tmp(:,:,:,:), T_vvoo(:,:,:,:)
  block_size = 16
  allocate(W_vvov(nV,nV,nO,block_size), W_vvov_tmp(nV,nO,nV,block_size), T_vvoo(nV,nV,nO,nO))

  !$omp parallel &
  !$omp private(u,i,b,a) &
  !$omp default(shared)
  !$omp do
  do u = 1, nO
    do i = 1, nO
      do b = 1, nV
        do a = 1, nV
          T_vvoo(a,b,i,u) = tau(i,u,a,b)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  do iblock = 1, nV, block_size
    nVmax = min(block_size,nV-iblock+1)

    call dgemm('T','N', nV*nO, nV*nVmax, cholesky_mo_num, 1.d0, &
      cc_space_v_vo_chol            , cholesky_mo_num, &
      cc_space_v_vv_chol(1,1,iblock), cholesky_mo_num, &
      0.d0, W_vvov_tmp, nV*nO)

    !$omp parallel &
    !$omp private(b,i,a,beta) &
    !$omp default(shared)
    do beta = 1,  nVmax
      do i = 1, nO
        !$omp do
        do b = 1, nV
          do a = 1, nV
            W_vvov(a,b,i,beta) = 2d0 * W_vvov_tmp(a,i,b,beta) - W_vvov_tmp(b,i,a,beta)
          enddo
        enddo
        !$omp end do nowait
      enddo
    enddo
    !$omp barrier
    !$omp end parallel

    call dgemm('T','N',nO,nVmax,nO*nV*nV, &
             1d0, T_vvoo, nV*nV*nO, &
                  W_vvov, nO*nV*nV, &
             1d0, r1(1,iblock), nO)
  enddo

  deallocate(W_vvov,T_vvoo)


  double precision, allocatable :: W_oovo(:,:,:,:)
  allocate(W_oovo(nO,nO,nV,nO))

  !$omp parallel &
  !$omp shared(nO,nV,cc_space_v_oovo,W_oovo) &
  !$omp private(u,a,i,j) &
  !$omp default(none)
  do u = 1, nO
    !$omp do
    do a = 1, nV
      do j = 1, nO
        do i = 1, nO
!          W_oovo(i,j,a,u) = 2d0 * cc_space_v_vooo(a,u,i,j) - cc_space_v_vooo(a,u,j,i)
          W_oovo(i,j,a,u) = 2d0 * cc_space_v_oovo(i,j,a,u) - cc_space_v_oovo(j,i,a,u)
        enddo
      enddo
    enddo
    !$omp end do nowait
  enddo
  !$omp end parallel

  call dgemm('T','N', nO, nV, nO*nO*nV, &
             -1d0, W_oovo, size(W_oovo,1) * size(W_oovo,2) * size(W_oovo,3), &
                   tau   , size(tau,1) * size(tau,2) * size(tau,3), &
              1d0, r1    , size(r1,1))

  deallocate(W_oovo)

  max_r1 = 0d0
  do a = 1, nV
    do i = 1, nO
      max_r1 = max(dabs(r1(i,a)), max_r1)
    enddo
  enddo

  ! Change the sign for consistency with the code in spin orbitals
  !$omp parallel &
  !$omp shared(nO,nV,r1) &
  !$omp private(a,i) &
  !$omp default(none)
  !$omp do
  do a = 1, nV
    do i = 1, nO
      r1(i,a) = -r1(i,a)
    enddo
  enddo
  !$omp end do
  !$omp end parallel

end

! H_oo

subroutine compute_H_oo_chol(nO,nV,tau_x,H_oo)

  implicit none

  integer, intent(in)           :: nO,nV
  double precision, intent(in)  :: tau_x(nO, nO, nV, nV)
  double precision, intent(out) :: H_oo(nO, nO)

  integer :: a,b,i,j,u,k

  double precision, allocatable :: tau_kau(:,:,:), tmp_vov(:,:,:)

  allocate(tau_kau(cholesky_mo_num,nV,nO))
  !$omp parallel &
  !$omp default(shared) &
  !$omp private(i,u,j,k,a,b,tmp_vov)
  allocate(tmp_vov(nV,nO,nV) )
  !$omp do
  do u = 1, nO
    do b=1,nV
      do j=1,nO
        do a=1,nV
          tmp_vov(a,j,b) = tau_x(u,j,a,b)
        enddo
      enddo
    enddo
    call dgemm('N','T',cholesky_mo_num,nV,nO*nV,1.d0, &
      cc_space_v_ov_chol, cholesky_mo_num, tmp_vov, nV, &
      0.d0, tau_kau(1,1,u), cholesky_mo_num)
  enddo
  !$omp end do nowait
  deallocate(tmp_vov)
  !$omp do
  do i = 1, nO
    do u = 1, nO
      H_oo(u,i) = cc_space_f_oo(u,i)
    enddo
  enddo
  !$omp end do nowait
  !$omp barrier
  !$omp end  parallel
  call dgemm('T', 'N', nO, nO, cholesky_mo_num*nV, 1.d0, &
    tau_kau, cholesky_mo_num*nV,  cc_space_v_vo_chol, cholesky_mo_num*nV, &
    1.d0, H_oo, nO)

end

! H_vv

subroutine compute_H_vv_chol(nO,nV,tau_x,H_vv)

  implicit none

  integer, intent(in)           :: nO,nV
  double precision, intent(in)  :: tau_x(nO, nO, nV, nV)
  double precision, intent(out) :: H_vv(nV, nV)

  integer :: a,b,i,j,u,k, beta

  double precision, allocatable :: tau_kia(:,:,:), tmp_oov(:,:,:)

  allocate(tau_kia(cholesky_mo_num,nO,nV))
  !$omp parallel &
  !$omp default(shared) &
  !$omp private(i,beta,j,k,a,b,tmp_oov)
  allocate(tmp_oov(nO,nO,nV) )
  !$omp do
  do a = 1, nV
    do b=1,nV
      do j=1,nO
        do i=1,nO
          tmp_oov(i,j,b) = tau_x(i,j,a,b)
        enddo
      enddo
    enddo
    call dgemm('N','T',cholesky_mo_num,nO,nO*nV,1.d0, &
      cc_space_v_ov_chol, cholesky_mo_num, tmp_oov, nO, &
      0.d0, tau_kia(1,1,a), cholesky_mo_num)
  enddo
  !$omp end do nowait
  deallocate(tmp_oov)

  !$omp do
  do beta = 1, nV
    do a = 1, nV
      H_vv(a,beta) = cc_space_f_vv(a,beta)
    enddo
  enddo
  !$omp end do nowait
  !$omp barrier
  !$omp end  parallel
  call dgemm('T', 'N', nV, nV, cholesky_mo_num*nO, -1.d0, &
    tau_kia, cholesky_mo_num*nO,  cc_space_v_ov_chol, cholesky_mo_num*nO, &
    1.d0, H_vv, nV)

end

! H_vo
subroutine compute_H_vo_chol(nO,nV,t1,H_vo)

  implicit none

  integer, intent(in)           :: nO,nV
  double precision, intent(in)  :: t1(nO, nV)
  double precision, intent(out) :: H_vo(nV, nO)

  integer :: a,b,i,j,u,k

  double precision, allocatable :: tmp_k(:), tmp(:,:,:), tmp2(:,:,:)
  do i=1,nO
    do a=1,nV
      H_vo(a,i) = cc_space_f_vo(a,i)
    enddo
  enddo

  allocate(tmp_k(cholesky_mo_num))
  call dgemm('N', 'N', cholesky_mo_num, 1, nO*nV, 2.d0, &
     cc_space_v_ov_chol, cholesky_mo_num, &
     t1, nO*nV, 0.d0, tmp_k, cholesky_mo_num)

  call dgemm('T','N',nV*nO,1,cholesky_mo_num,1.d0, &
      cc_space_v_vo_chol, cholesky_mo_num, tmp_k, cholesky_mo_num, 1.d0, &
      H_vo, nV*nO)
  deallocate(tmp_k)

  allocate(tmp(cholesky_mo_num,nO,nO))
  allocate(tmp2(cholesky_mo_num,nO,nO))

  call dgemm('N','T', cholesky_mo_num*nO, nO, nV, 1.d0, &
    cc_space_v_ov_chol, cholesky_mo_num*nO, t1, nO, 0.d0, tmp, cholesky_mo_num*nO)

  do i=1,nO
    do j=1,nO
      do k=1,cholesky_mo_num
        tmp2(k,j,i) = tmp(k,i,j)
      enddo
    enddo
  enddo
  deallocate(tmp)

  call dgemm('T','N', nV, nO, cholesky_mo_num*nO, -1.d0, &
    cc_space_v_ov_chol, cholesky_mo_num*nO, tmp2, cholesky_mo_num*nO, &
    1.d0, H_vo, nV)

end


! R2

subroutine compute_r2_space_chol(nO,nV,t1,t2,tau,H_oo,H_vv,H_vo,r2,max_r2)

  implicit none

  ! in
  integer, intent(in)           :: nO, nV
  double precision, intent(in)  :: t1(nO,nV), t2(nO,nO,nV,nV), tau(nO,nO,nV,nV)
  double precision, intent(in)  :: H_oo(nO,nO), H_vv(nV,nV), H_vo(nV,nO)

  ! out
  double precision, intent(out) :: r2(nO,nO,nV,nV), max_r2

  ! internal
  integer                       :: u,v,i,j,beta,gam,a,b
  double precision              :: max_r2_local

  call set_multiple_levels_omp(.False.)

  !$omp parallel &
  !$omp shared(nO,nV,r2,cc_space_v_oovv) &
  !$omp private(u,v,gam,beta) &
  !$omp default(none)
  !$omp do
  do gam = 1, nV
    do beta = 1, nV
      do v = 1, nO
        do u = 1, nO
         r2(u,v,beta,gam) = cc_space_v_oovv(u,v,beta,gam)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  double precision, allocatable :: A1(:,:,:,:)
  allocate(A1(nO,nO,nO,nO))
  call compute_A1_chol(nO,nV,t1,t2,tau,A1)
  call dgemm('N','N',nO*nO,nV*nV,nO*nO, &
             1d0, A1, size(A1,1) * size(A1,2), &
                  tau, size(tau,1) * size(tau,2), &
             1d0, r2, size(r2,1) * size(r2,2))

  deallocate(A1)
  integer :: block_size, iblock, k
  block_size = 16
  double precision, dimension(:,:,:), allocatable :: B1, tmp_cc, tmpB1
  double precision, dimension(:,:), allocatable :: tmp_cc2

  allocate(tmp_cc(cholesky_mo_num,nV,nV))
  call dgemm('N','N', cholesky_mo_num*nV, nV, nO, 1.d0, &
      cc_space_v_vo_chol, cholesky_mo_num*nV, t1, nO, 0.d0, tmp_cc, cholesky_mo_num*nV)

  call set_multiple_levels_omp(.False.)

  !$OMP PARALLEL PRIVATE(gam, iblock, B1, tmpB1, tmp_cc2, beta, b, a)
  allocate(B1(nV,nV,block_size), tmpB1(nV,block_size,nV), tmp_cc2(cholesky_mo_num,nV))
  !$OMP DO
  do gam = 1, nV

    do a=1,nV
      do k=1,cholesky_mo_num
        tmp_cc2(k,a) = cc_space_v_vv_chol(k,a,gam) - tmp_cc(k,a,gam)
      enddo
    enddo

    do iblock = 1, nV, block_size

        call dgemm('T', 'N', nV*min(block_size, nV-iblock+1), nV, cholesky_mo_num, &
                -1.d0, tmp_cc(1,1,iblock), cholesky_mo_num, &
                cc_space_v_vv_chol(1,1,gam), cholesky_mo_num, &
                0.d0, tmpB1, nV*block_size)

        call dgemm('T','N', nV*min(block_size, nV-iblock+1), nV, cholesky_mo_num, &
                1.d0, cc_space_v_vv_chol(1,1,iblock), cholesky_mo_num, &
                tmp_cc2, cholesky_mo_num, &
                1.d0, tmpB1, nV*block_size)

        do beta = iblock, min(nV, iblock+block_size-1)
          do b = 1, nV
            do a = 1, nV
              B1(a,b,beta-iblock+1) = tmpB1(a,beta-iblock+1,b)
            enddo
          enddo
        enddo

        call dgemm('N','N',nO*nO,min(block_size, nV-iblock+1),nV*nV, &
              1d0, tau, size(tau,1) * size(tau,2), &
                   B1 , size(B1 ,1) * size(B1 ,2), &
              1d0, r2(1,1,iblock,gam),  size(r2 ,1) * size(r2 ,2))
      enddo

  enddo
  !$OMP ENDDO

  deallocate(B1, tmpB1, tmp_cc2)
  !$OMP END PARALLEL

  deallocate(tmp_cc)


  double precision, allocatable :: X_oovv(:,:,:,:)
  allocate(X_oovv(nO,nO,nV,nV))
  !$omp parallel &
  !$omp shared(nO,nV,t2,X_oovv) &
  !$omp private(u,v,gam,a) &
  !$omp default(none)
  !$omp do
  do a = 1, nV
    do gam = 1, nV
      do v = 1, nO
        do u = 1, nO
          X_oovv(u,v,gam,a) = t2(u,v,gam,a)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  double precision, allocatable :: g_vir(:,:)
  allocate(g_vir(nV,nV))
  call compute_g_vir_chol(nO,nV,t1,t2,H_vv,g_vir)

  double precision, allocatable :: Y_oovv(:,:,:,:)
  allocate(Y_oovv(nO,nO,nV,nV))

  call dgemm('N','N',nO*nO*nV,nV,nV, &
             1d0, X_oovv, size(X_oovv,1) * size(X_oovv,2) * size(X_oovv,3), &
                  g_vir, size(g_vir,1), &
             0d0, Y_oovv, size(Y_oovv,1) * size(Y_oovv,2) * size(Y_oovv,3))
  deallocate(g_vir)
  deallocate(X_oovv)

  !$omp parallel &
  !$omp shared(nO,nV,r2,Y_oovv) &
  !$omp private(u,v,gam,beta) &
  !$omp default(none)
  !$omp do
  do gam = 1, nV
    do beta = 1, nV
      do v = 1, nO
        do u = 1, nO
          r2(u,v,beta,gam) = r2(u,v,beta,gam) + Y_oovv(u,v,beta,gam) + Y_oovv(v,u,gam,beta)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel
  deallocate(Y_oovv)

  double precision, allocatable :: g_occ(:,:)
  allocate(g_occ(nO,nO))
  call compute_g_occ_chol(nO,nV,t1,t2,H_oo,g_occ)

  allocate(X_oovv(nO,nO,nV,nV))
  call dgemm('N','N',nO,nO*nV*nV,nO, &
             1d0, g_occ , size(g_occ,1), &
                  t2    , size(t2,1), &
             0d0, X_oovv, size(X_oovv,1))
  deallocate(g_occ)

  !$omp parallel &
  !$omp shared(nO,nV,r2,X_oovv) &
  !$omp private(u,v,gam,beta) &
  !$omp default(none)
  !$omp do
  do gam = 1, nV
    do beta = 1, nV
      do v = 1, nO
        do u = 1, nO
          r2(u,v,beta,gam) = r2(u,v,beta,gam) - X_oovv(u,v,beta,gam) - X_oovv(v,u,gam,beta)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  deallocate(X_oovv)

  double precision, allocatable :: X_vovv(:,:,:,:)

  allocate(X_vovv(nV,nO,nV,block_size))
  allocate(Y_oovv(nO,nO,nV,nV))

  do iblock = 1, nV, block_size
    do gam = iblock, min(nV, iblock+block_size-1)
      call dgemm('T','N',nV, nO*nV, cholesky_mo_num, 1.d0, &
        cc_space_v_vv_chol(1,1,gam), cholesky_mo_num, cc_space_v_ov_chol, &
        cholesky_mo_num, 0.d0, X_vovv(1,1,1,gam-iblock+1), nV)

    enddo
    call dgemm('N','N',nO,nO*nV*min(block_size, nV-iblock+1),nV, &
             1d0, t1    , size(t1,1), &
                  X_vovv, size(X_vovv,1), &
             0d0, Y_oovv(1,1,1,iblock), size(Y_oovv,1))

  enddo
  deallocate(X_vovv)

  !$omp parallel &
  !$omp shared(nO,nV,r2,Y_oovv) &
  !$omp private(u,v,gam,beta) &
  !$omp default(none)
  !$omp do
  do gam = 1, nV
    do beta = 1, nV
      do v = 1, nO
        do u = 1, nO
          r2(u,v,beta,gam) = r2(u,v,beta,gam) + Y_oovv(v,u,beta,gam) + Y_oovv(u,v,gam,beta)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel
  deallocate(Y_oovv)

  double precision, allocatable :: X_ovvo(:,:,:,:)
  double precision, allocatable :: tcc(:,:,:), tcc2(:,:,:)
  allocate(tcc2(cholesky_mo_num,nV,nO), X_ovvo(nO,nV,nV,nO))
  allocate(tcc(cholesky_mo_num,nO,nV))

  call dgemm('N','T', cholesky_mo_num*nV, nO, nV, 1.d0, &
     cc_space_v_vv_chol, cholesky_mo_num*nV, t1, nO, &
     0.d0, tcc2, cholesky_mo_num*nV)

  call dgemm('N','N', cholesky_mo_num*nO, nV, nO, 1.d0, &
     cc_space_v_oo_chol, cholesky_mo_num*nO, t1, nO, &
     0.d0, tcc, cholesky_mo_num*nO)

  call dgemm('T','N', nO*nV, nV*nO, cholesky_mo_num, 1.d0, &
              tcc, cholesky_mo_num, tcc2, cholesky_mo_num, 0.d0, &
              X_ovvo, nO*nV)

  deallocate(tcc, tcc2)

  !$omp parallel &
  !$omp shared(nO,nV,r2,X_ovvo) &
  !$omp private(u,v,gam,beta) &
  !$omp default(none)
  !$omp do
  do gam = 1, nV
    do beta = 1, nV
      do v = 1, nO
        do u = 1, nO
          r2(u,v,beta,gam) = r2(u,v,beta,gam) - X_ovvo(u,beta,gam,v)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp do
  do beta = 1, nV
    do gam = 1, nV
      do v = 1, nO
        do u = 1, nO
          r2(v,u,gam,beta) = r2(v,u,gam,beta) - X_ovvo(u,beta,gam,v)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  deallocate(X_ovvo)
  !-----

  allocate(X_oovv(nO,nO,nV,nV))

  call dgemm('N','N',nO*nO*nV,nV,nO, &
             1d0, cc_space_v_oovo, size(cc_space_v_oovo,1) * size(cc_space_v_oovo,2) * size(cc_space_v_oovo,3), &
                  t1 , size(t1,1), &
             0d0, X_oovv, size(X_oovv,1) * size(X_oovv,2) * size(X_oovv,3))

  !$omp parallel &
  !$omp shared(nO,nV,r2,X_oovv) &
  !$omp private(u,v,gam,beta) &
  !$omp default(none)
  !$omp do
  do gam = 1, nV
    do beta = 1, nV
      do v = 1, nO
        do u = 1, nO
           r2(u,v,beta,gam) = r2(u,v,beta,gam) - X_oovv(u,v,beta,gam) - X_oovv(v,u,gam,beta)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel
  deallocate(X_oovv)

  double precision, allocatable :: X_vovo(:,:,:,:), Y_oovo(:,:,:,:)
  allocate(X_vovo(nV,nO,nV,nO))

  !$omp parallel &
  !$omp shared(nO,nV,X_vovo,cc_space_v_ovvo) &
  !$omp private(a,v,gam,i) &
  !$omp default(none)
  do i = 1, nO
    !$omp do
    do gam = 1, nV
      do v = 1, nO
        do a = 1, nV
          X_vovo(a,v,gam,i) = cc_space_v_ovvo(v,a,gam,i)
        enddo
      enddo
    enddo
    !$omp end do nowait
  enddo
  !$omp end parallel

  allocate(Y_oovo(nO,nO,nV,nO))
  call dgemm('N','N',nO,nO*nV*nO,nV, &
             1d0, t1, size(t1,1), &
                  X_vovo, size(X_vovo,1), &
             0d0, Y_oovo, size(Y_oovo,1))

  deallocate(X_vovo)
  allocate(X_oovv(nO,nO,nV,nV))
  call dgemm('N','N',nO*nO*nV, nV, nO, &
             1d0, Y_oovo, size(Y_oovo,1) * size(Y_oovo,2) * size(Y_oovo,3), &
                  t1    , size(t1,1), &
             0d0, X_oovv, size(X_oovv,1) * size(X_oovv,2) * size(X_oovv,3))
  deallocate(Y_oovo)

  !$omp parallel &
  !$omp shared(nO,nV,r2,X_oovv) &
  !$omp private(u,v,gam,beta) &
  !$omp default(none)
  !$omp do
  do gam = 1, nV
    do beta = 1, nV
      do v = 1, nO
        do u = 1, nO
           r2(u,v,beta,gam) = r2(u,v,beta,gam) - X_oovv(u,v,gam,beta) - X_oovv(v,u,beta,gam)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel
  deallocate(X_oovv)


  double precision, allocatable :: J1(:,:,:,:)
  allocate(J1(nO,nV,nV,nO))
  call compute_J1_chol(nO,nV,t1,t2,cc_space_v_ovvo,cc_space_v_ovoo, &
       cc_space_v_vvoo,J1)

  double precision, allocatable :: K1(:,:,:,:)
  allocate(K1(nO,nV,nO,nV))
  call compute_K1_chol(nO,nV,t1,t2,cc_space_v_ovoo,cc_space_v_vvoo, &
       cc_space_v_ovov,K1)

  allocate(X_ovvo(nO,nV,nV,nO))
  !$omp parallel &
  !$omp private(u,v,gam,beta,i,a) &
  !$omp default(shared)
  do i = 1, nO
    !$omp do
    do a = 1, nV
      do beta = 1, nV
        do u = 1, nO
          X_ovvo(u,beta,a,i) = (J1(u,a,beta,i) - 0.5d0 * K1(u,a,i,beta))
        enddo
      enddo
    enddo
    !$omp end do nowait
  enddo
  !$omp end parallel
  deallocate(J1)

  double precision, allocatable :: Y_voov(:,:,:,:)
  allocate(Y_voov(nV,nO,nO,nV))

  !$omp parallel &
  !$omp private(u,v,gam,beta,i,a) &
  !$omp default(shared)
  !$omp do
  do gam = 1, nV
    do v = 1, nO
      do i = 1, nO
        do a = 1, nV
          Y_voov(a,i,v,gam) = 2d0 * t2(i,v,a,gam) - t2(i,v,gam,a)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  double precision, allocatable :: Z_ovov(:,:,:,:)
  allocate(Z_ovov(nO,nV,nO,nV))

  call dgemm('N','N', nO*nV,nO*nV,nV*nO, &
             1d0, X_ovvo, size(X_ovvo,1) * size(X_ovvo,2), &
                  Y_voov, size(Y_voov,1) * size(Y_voov,2), &
             0d0, Z_ovov, size(Z_ovov,1) * size(Z_ovov,2))

  deallocate(X_ovvo,Y_voov)

  !$omp parallel &
  !$omp shared(nO,nV,r2,Z_ovov) &
  !$omp private(u,v,gam,beta) &
  !$omp default(none)
  !$omp do
  do gam = 1, nV
    do beta = 1, nV
      do v = 1, nO
        do u = 1, nO
           r2(u,v,beta,gam) = r2(u,v,beta,gam) + Z_ovov(u,beta,v,gam) + Z_ovov(v,gam,u,beta)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  deallocate(Z_ovov)

  double precision, allocatable :: Y_ovov(:,:,:,:), X_ovov(:,:,:,:)
  allocate(X_ovov(nO,nV,nO,nV))
  allocate(Y_ovov(nO,nV,nO,nV))

  !$omp parallel &
  !$omp shared(nO,nV,r2,K1,X_ovov,Y_ovov,t2) &
  !$omp private(u,a,i,beta,gam) &
  !$omp default(none)
  !$omp do
  do beta = 1, nV
    do u = 1, nO
      do a = 1, nV
        do i = 1, nO
          X_ovov(i,a,u,beta) = 0.5d0 * K1(u,a,i,beta)
        enddo
      enddo
    enddo
  enddo
  !$omp end do nowait

  !$omp do
  do gam = 1, nV
    do v = 1, nO
      do a = 1, nV
        do i = 1, nO
          Y_ovov(i,a,v,gam) = t2(i,v,gam,a)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  allocate(Z_ovov(nO,nV,nO,nV))
  call dgemm('T','N',nO*nV,nO*nV,nO*nV, &
             1d0, X_ovov, size(X_ovov,1) * size(X_ovov,2), &
                  Y_ovov, size(Y_ovov,1) * size(Y_ovov,2), &
             0d0, Z_ovov, size(Y_ovov,1) * size(Y_ovov,2))
  deallocate(X_ovov, Y_ovov)

  !$omp parallel &
  !$omp shared(nO,nV,r2,Z_ovov) &
  !$omp private(u,v,gam,beta) &
  !$omp default(none)
  !$omp do
  do gam = 1, nV
    do beta = 1, nV
      do v = 1, nO
        do u = 1, nO
           r2(u,v,beta,gam) = r2(u,v,beta,gam) - Z_ovov(u,beta,v,gam) - Z_ovov(v,gam,u,beta)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel
  deallocate(Z_ovov)

  allocate(X_ovov(nO,nV,nO,nV),Y_ovov(nO,nV,nO,nV))
  !$omp parallel &
  !$omp shared(nO,nV,K1,X_ovov,Y_ovov,t2) &
  !$omp private(u,v,gam,beta,i,a) &
  !$omp default(none)
  !$omp do
  do a = 1, nV
    do i = 1, nO
      do gam = 1, nV
        do u = 1, nO
          X_ovov(u,gam,i,a) = K1(u,a,i,gam)
        enddo
      enddo
    enddo
  enddo
  !$omp end do nowait

  !$omp do
  do beta = 1, nV
    do v = 1, nO
      do a = 1, nV
        do i = 1, nO
          Y_ovov(i,a,v,beta) = t2(i,v,beta,a)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  deallocate(K1)

  allocate(Z_ovov(nO,nV,nO,nV))
  call dgemm('N','N',nO*nV,nO*nV,nO*nV, &
             1d0, X_ovov, size(X_ovov,1) * size(X_ovov,2), &
                  Y_ovov, size(Y_ovov,1) * size(Y_ovov,2), &
             0d0, Z_ovov, size(Y_ovov,1) * size(Y_ovov,2))

  deallocate(X_ovov,Y_ovov)

  !$omp parallel &
  !$omp shared(nO,nV,r2,Z_ovov) &
  !$omp private(u,v,gam,beta) &
  !$omp default(none)
  !$omp do
  do gam = 1, nV
    do beta = 1, nV
      do v = 1, nO
        do u = 1, nO
           r2(u,v,beta,gam) = r2(u,v,beta,gam) - Z_ovov(u,gam,v,beta) - Z_ovov(v,beta,u,gam)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  deallocate(Z_ovov)

  ! Change the sign for consistency with the code in spin orbitals

  max_r2 = 0d0
  !$omp parallel &
  !$omp shared(nO,nV,r2,max_r2) &
  !$omp private(i,j,a,b,max_r2_local) &
  !$omp default(none)
  max_r2_local = 0.d0
  !$omp do
  do b = 1, nV
    do a = 1, nV
      do j = 1, nO
        do i = 1, nO
          r2(i,j,a,b) = -r2(i,j,a,b)
          max_r2_local = max(r2(i,j,a,b), max_r2_local)
        enddo
      enddo
    enddo
  enddo
  !$omp end do nowait
  !$omp critical
  max_r2 = max(max_r2, max_r2_local)
  !$omp end critical
  !$omp end parallel

end

! A1

subroutine compute_A1_chol(nO,nV,t1,t2,tau,A1)

  implicit none

  integer, intent(in)           :: nO,nV
  double precision, intent(in)  :: t1(nO, nV)
  double precision, intent(in)  :: t2(nO, nO, nV, nV)
  double precision, intent(in)  :: tau(nO, nO, nV, nV)
  double precision, intent(out) :: A1(nO, nO, nO, nO)

  integer :: a,tmp_a,b,k,l,c,d,tmp_c,tmp_d,i,j,u,v, beta

  double precision, allocatable :: Y_oooo(:,:,:,:)
  allocate(Y_oooo(nO,nO,nO,nO))

  ! A1(u,v,i,j) = cc_space_v_oooo(u,v,i,j)
  ! A1(u,v,i,j) += cc_space_v_ovoo(u,a,i,j) * t1(v,a) &

  call dgemm('N','N', nO, nO*nO*nO, nV, &
             1d0, t1    , size(t1,1), &
                  cc_space_v_vooo, size(cc_space_v_vooo,1), &
             0d0, Y_oooo, size(Y_oooo,1))

  !$omp parallel &
  !$omp private(u,v,i,j) &
  !$omp default(shared)
  !$omp do collapse(2)
  do j = 1, nO
    do i = 1, nO
      do v = 1, nO
        do u = 1, nO
          A1(u,v,i,j) = cc_space_v_oooo(u,v,i,j) + Y_oooo(v,u,j,i) + Y_oooo(u,v,i,j)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  deallocate(Y_oooo)

  ! A1(u,v,i,j) += cc_space_v_vvoo(a,b,i,j) * tau(u,v,a,b)
  call dgemm('N','N', nO*nO, nO*nO, nV*nV, &
             1d0, tau     , size(tau,1) * size(tau,2), &
                  cc_space_v_vvoo, size(cc_space_v_vvoo,1) * size(cc_space_v_vvoo,2), &
             1d0, A1      , size(A1,1) * size(A1,2))

end

! g_occ

subroutine compute_g_occ_chol(nO,nV,t1,t2,H_oo,g_occ)

  implicit none

  integer, intent(in)           :: nO,nV
  double precision, intent(in)  :: t1(nO, nV), H_oo(nO, nO)
  double precision, intent(in)  :: t2(nO, nO, nV, nV)
  double precision, intent(out) :: g_occ(nO, nO)

  g_occ = H_oo

  call dgemm('N','N',nO,nO,nV, &
             1d0, t1, size(t1,1), &
                  cc_space_f_vo, size(cc_space_f_vo,1), &
             1d0, g_occ, size(g_occ,1))

  double precision, allocatable :: X(:)
  allocate(X(cholesky_mo_num))
  call dgemv('N',cholesky_mo_num,nO*nV,2.d0, &
    cc_space_v_ov_chol, cholesky_mo_num, &
    t1, 1, 0.d0, X, 1)

  call dgemv('T',cholesky_mo_num,nO*nO,1.d0, &
    cc_space_v_oo_chol, cholesky_mo_num, &
    X, 1, 1.d0, g_occ, 1)
  deallocate(X)

  call dgemv('T',nO*nV,nO*nO,-1.d0, &
    cc_space_v_ovoo, nO*nV, &
    t1, 1, 1.d0, g_occ, 1)

end

! g_vir

subroutine compute_g_vir_chol(nO,nV,t1,t2,H_vv,g_vir)

  implicit none

  integer, intent(in)           :: nO,nV
  double precision, intent(in)  :: t1(nO, nV), H_vv(nV, nV)
  double precision, intent(in)  :: t2(nO, nO, nV, nV)
  double precision, intent(out) :: g_vir(nV, nV)

  integer :: a,tmp_a,b,k,l,c,d,tmp_c,tmp_d,i,j,u,v, beta, gam

  call dgemm('N','N',nV,nV,nO, &
             -1d0, cc_space_f_vo , size(cc_space_f_vo,1), &
                   t1   , size(t1,1), &
              0d0, g_vir, size(g_vir,1))

  double precision, allocatable :: tmp_k(:), tmp_vo(:,:,:), tmp_vo2(:,:,:)
  allocate(tmp_k(cholesky_mo_num))
  call dgemm('N','N', cholesky_mo_num, 1, nO*nV, 1.d0, &
    cc_space_v_ov_chol, cholesky_mo_num, t1, nO*nV, 0.d0, tmp_k, cholesky_mo_num)

  call dgemm('T','N', nV*nV, 1, cholesky_mo_num, 2.d0, &
    cc_space_v_vv_chol, cholesky_mo_num, tmp_k, cholesky_mo_num, 1.d0, &
    g_vir, nV*nV)
  deallocate(tmp_k)

  allocate(tmp_vo(cholesky_mo_num,nV,nO))
  call dgemm('N','T',cholesky_mo_num*nV, nO, nV, 1.d0, &
    cc_space_v_vv_chol, cholesky_mo_num*nV, t1, nO, 0.d0, tmp_vo, cholesky_mo_num*nV)

  allocate(tmp_vo2(cholesky_mo_num,nO,nV))
  do beta=1,nV
    do i=1,nO
      do k=1,cholesky_mo_num
        tmp_vo2(k,i,beta) = -tmp_vo(k,beta,i)
      enddo
    enddo
  enddo
  deallocate(tmp_vo)

  do beta = 1, nV
    do a = 1, nV
      g_vir(a,beta) = g_vir(a,beta) + H_vv(a,beta)
    enddo
  enddo

  call dgemm('T','N', nV, nV, nO*cholesky_mo_num, 1.d0, &
     cc_space_v_ov_chol, cholesky_mo_num*nO, &
     tmp_vo2, cholesky_mo_num*nO, 1.d0, g_vir, nV)

end

! J1

subroutine compute_J1_chol(nO,nV,t1,t2,v_ovvo,v_ovoo,v_vvoo,J1)
  implicit none

  integer, intent(in)           :: nO,nV
  double precision, intent(in)  :: t1(nO, nV)
  double precision, intent(in)  :: t2(nO, nO, nV, nV)
  double precision, intent(in)  :: v_ovvo(nO,nV,nV,nO), v_ovoo(nO,nV,nO,nO)
  double precision, intent(in)  :: v_vvoo(nV,nV,nO,nO)
  double precision, intent(out) :: J1(nO, nV, nV, nO)

  integer :: a,tmp_a,b,k,l,c,d,tmp_c,tmp_d,i,j,u,v, beta, gam

  double precision, allocatable :: X_ovoo(:,:,:,:), Y_ovov(:,:,:,:)
  allocate(X_ovoo(nO,nV,nO,nO),Y_ovov(nO,nV,nO,nV))

  !$omp parallel &
  !$omp shared(nO,nV,J1,v_ovvo,v_ovoo,X_ovoo) &
  !$omp private(i,j,a,u,beta) &
  !$omp default(none)
  do i = 1, nO
    !$omp do
    do beta = 1, nV
      do a = 1, nV
        do u = 1, nO
          J1(u,a,beta,i) = v_ovvo(u,a,beta,i)
        enddo
      enddo
    enddo
    !$omp end do nowait
  enddo

  !$omp do collapse(2)
  do j = 1, nO
    do i = 1, nO
      do a = 1, nV
        do u = 1, nO
          X_ovoo(u,a,i,j) = v_ovoo(u,a,j,i)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  call dgemm('N','N',nO*nV*nO,nV,nO, &
            -1d0, X_ovoo, size(X_ovoo,1) * size(X_ovoo,2) * size(X_ovoo,3), &
                  t1    , size(t1,1), &
             0d0, Y_ovov, size(Y_ovov,1) * size(Y_ovov,2) * size(Y_ovov,3))

  !$omp parallel &
  !$omp shared(nO,nV,J1,Y_ovov) &
  !$omp private(i,beta,a,u) &
  !$omp default(none)
  do i = 1, nO
    !$omp do
    do beta = 1, nV
      do a = 1, nV
        do u = 1, nO
          J1(u,a,beta,i) = J1(u,a,beta,i) + Y_ovov(u,a,i,beta)
        enddo
      enddo
    enddo
    !$omp end do nowait
  enddo
  !$omp end parallel
  deallocate(X_ovoo)

  double precision, allocatable :: tmp_cc(:,:,:), J1_tmp(:,:,:,:)
  allocate(tmp_cc(cholesky_mo_num,nV,nO), J1_tmp(nV,nO,nV,nO))

  call dgemm('N','T', cholesky_mo_num*nV, nO, nV, 1.d0, &
      cc_space_v_vv_chol, cholesky_mo_num*nV, &
      t1, nO, &
      0.d0, tmp_cc, cholesky_mo_num*nV)

  call dgemm('T','N', nV*nO, nV*nO, cholesky_mo_num, 1.d0, &
      tmp_cc, cholesky_mo_num, cc_space_v_vo_chol, cholesky_mo_num, &
      0.d0, J1_tmp, nV*nO)

  deallocate(tmp_cc)

  do i=1,nO
    do b=1,nV
      do a=1,nV
        do u=1,nO
          J1(u,a,b,i) = J1(u,a,b,i) + J1_tmp(b,u,a,i)
        enddo
      enddo
    enddo
  enddo

  deallocate(J1_tmp)

  !- cc_space_v_vvoo(a,b,i,j) * (0.5d0 * t2(u,j,b,beta) + t1(u,b) * t1(j,beta)) &
  double precision, allocatable :: X_voov(:,:,:,:), Z_ovvo(:,:,:,:)
  allocate(X_voov(nV,nO,nO,nV), Z_ovvo(nO,nV,nV,nO))
  !$omp parallel &
  !$omp shared(nO,nV,t2,t1,Y_ovov,X_voov,v_vvoo) &
  !$omp private(i,beta,a,u,b,j) &
  !$omp default(none)
  !$omp do
  do b = 1, nV
    do j = 1, nO
      do beta = 1, nV
        do u = 1, nO
          Y_ovov(u,beta,j,b) = 0.5d0 * t2(u,j,b,beta) + t1(u,b) * t1(j,beta)
        enddo
      enddo
    enddo
  enddo
  !$omp end do nowait

  !$omp do
  do b = 1, nV
    do j = 1, nO
      do i = 1, nO
        do a = 1, nV
          X_voov(a,i,j,b) = v_vvoo(a,b,i,j)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  call dgemm('N','T',nO*nV,nV*nO,nO*nV, &
             -1d0, Y_ovov, size(Y_ovov,1) * size(Y_ovov,2), &
                   X_voov, size(X_voov,1) * size(X_voov,2), &
              0d0, Z_ovvo, size(Z_ovvo,1) * size(Z_ovvo,2))
  deallocate(X_voov)

  double precision, allocatable :: X_ovvo(:,:,:,:), Y_vovo(:,:,:,:)
  allocate(X_ovvo(nO,nV,nV,nO),Y_vovo(nV,nO,nV,nO))
  !$omp parallel &
  !$omp shared(nO,nV,J1,Z_ovvo,t2,Y_vovo,v_vvoo,X_ovvo) &
  !$omp private(i,beta,a,u,j,b) &
  !$omp default(none)
  do i = 1, nO
    !$omp do
    do beta = 1, nV
      do a = 1, nV
        do u = 1, nO
          J1(u,a,beta,i) = J1(u,a,beta,i) + Z_ovvo(u,beta,a,i)
        enddo
      enddo
    enddo
    !$omp end do nowait
  enddo

  !+ 0.5d0 * (2d0 * cc_space_v_vvoo(a,b,i,j) - cc_space_v_vvoo(b,a,i,j)) * t2(u,j,beta,b)
  do j = 1, nO
    !$omp do
    do b = 1, nV
      do i = 1, nO
        do a = 1, nV
          Y_vovo(a,i,b,j) = 0.5d0 * (2d0 * v_vvoo(a,b,i,j) - v_vvoo(b,a,i,j))
        enddo
      enddo
    enddo
    !$omp end do nowait
  enddo

  do j = 1, nO
    !$omp do
    do b = 1, nV
      do beta = 1, nV
        do u = 1, nO
          X_ovvo(u,beta,b,j) = t2(u,j,beta,b)
        enddo
      enddo
    enddo
    !$omp end do nowait
  enddo
  !$omp end parallel

  call dgemm('N','T',nO*nV,nV*nO,nV*nO, &
             1d0, X_ovvo, size(X_ovvo,1) * size(X_ovvo,2), &
                  Y_vovo, size(Y_vovo,1) * size(Y_vovo,2), &
             0d0, Z_ovvo, size(Z_ovvo,1) * size(Z_ovvo,2))

  !$omp parallel &
  !$omp shared(nO,nV,J1,Z_ovvo) &
  !$omp private(i,beta,a,u) &
  !$omp default(none)
  do i = 1, nO
    !$omp do
    do beta = 1, nV
      do a = 1, nV
        do u = 1, nO
          J1(u,a,beta,i) = J1(u,a,beta,i) + Z_ovvo(u,beta,a,i)
        enddo
      enddo
    enddo
    !$omp end do nowait
  enddo
  !$omp end parallel

  deallocate(X_ovvo,Z_ovvo,Y_ovov)

end

! K1

subroutine compute_K1_chol(nO,nV,t1,t2,v_ovoo,v_vvoo,v_ovov,K1)

  implicit none

  integer, intent(in)           :: nO,nV
  double precision, intent(in)  :: t1(nO, nV)
  double precision, intent(in)  :: t2(nO, nO, nV, nV)
  double precision, intent(in)  :: v_vvoo(nV,nV,nO,nO), v_ovov(nO,nV,nO,nV)
  double precision, intent(in)  :: v_ovoo(nO,nV,nO,nO)
  double precision, intent(out) :: K1(nO, nV, nO, nV)

  double precision, allocatable :: X(:,:,:,:), Y(:,:,:,:), Z(:,:,:,:)

  integer :: a,tmp_a,b,k,l,c,d,tmp_c,tmp_d,i,j,u,v, beta, gam

  allocate(X(nV,nO,nV,nO),Y(nO,nV,nV,nO),Z(nO,nV,nV,nO))

  !$omp parallel &
  !$omp shared(nO,nV,K1,X,Y,v_vvoo,v_ovov,t1,t2) &
  !$omp private(i,beta,a,u,j,b) &
  !$omp default(none)
  !$omp do
  do beta = 1, nV
    do i = 1, nO
      do a = 1, nV
        do u = 1, nO
          K1(u,a,i,beta) = v_ovov(u,a,i,beta)
        enddo
      enddo
    enddo
  enddo
  !$omp end do nowait

  do i = 1, nO
    !$omp do
    do a = 1, nV
      do j = 1, nO
        do b = 1, nV
          X(b,j,a,i) = - v_vvoo(b,a,i,j)
        enddo
      enddo
    enddo
    !$omp end do nowait
  enddo

  do j = 1, nO
    !$omp do
    do b = 1, nV
      do beta = 1, nV
        do u = 1, nO
          Y(u,beta,b,j) = 0.5d0 * t2(u,j,b,beta) + t1(u,b) * t1(j,beta)
        enddo
      enddo
    enddo
    !$omp end do
  enddo
  !$omp end parallel

  call dgemm('N','N',nO*nV*nO,nV,nO, &
            -1d0, v_ovoo, size(v_ovoo,1) * size(v_ovoo,2) * size(v_ovoo,3), &
                  t1    , size(t1,1), &
            1d0, K1    , size(K1,1) * size(K1,2) * size(K1,3))

  double precision, allocatable :: K1tmp(:,:,:,:), t1v(:,:,:)
  allocate(K1tmp(nO,nO,nV,nV), t1v(cholesky_mo_num,nO,nO))

  call dgemm('N','T', cholesky_mo_num*nO, nO, nV, 1.d0, &
    cc_space_v_ov_chol, cholesky_mo_num*nO, t1, nO, 0.d0, &
    t1v, cholesky_mo_num*nO)

  call dgemm('T','N', nO*nO, nV*nV, cholesky_mo_num, 1.d0, &
    t1v, cholesky_mo_num, cc_space_v_vv_chol, cholesky_mo_num, 0.d0, &
    K1tmp, nO*nO)

  deallocate(t1v)
  ! Y(u,beta,b,j) * X(b,j,a,i) = Z(u,beta,a,i)
  call dgemm('N','N',nV*nO,nO*nV,nV*nO, &
             1d0, Y, size(Y,1) * size(Y,2), &
                  X, size(X,1) * size(X,2), &
             0d0, Z, size(Z,1) * size(Z,2))

  !$omp parallel &
  !$omp shared(nO,nV,K1,Z,K1tmp) &
  !$omp private(i,beta,a,u) &
  !$omp default(none)
  !$omp do
   do beta = 1, nV
    do i = 1, nO
      do a = 1, nV
        do u = 1, nO
          K1(u,a,i,beta) = K1(u,a,i,beta) + K1tmp(u,i,a,beta) + Z(u,beta,a,i)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  deallocate(K1tmp,X,Y,Z)

end
