! H_oo

subroutine compute_H_oo_chol(nO,nV,tau_x,d_cc_space_f_oo, &
    d_cc_space_v_ov_chol,d_cc_space_v_vo_chol,H_oo)
  use gpu
  implicit none

  integer, intent(in)           :: nO,nV
  type(gpu_double2), intent(in)    :: d_cc_space_f_oo
  type(gpu_double3), intent(in)    :: d_cc_space_v_ov_chol, d_cc_space_v_vo_chol
  type(gpu_double4), intent(in)    :: tau_x
  type(gpu_double2), intent(out)   :: H_oo

  integer :: a,b,i,j,u,k

  type(gpu_double3) :: tau_kau, tmp_vov, tmp_ovv

  call gpu_allocate(tau_kau, cholesky_mo_num, nV, nO)

  type(gpu_blas) :: blas


  !$OMP PARALLEL  &
  !$OMP DEFAULT(SHARED) &
  !$OMP PRIVATE(blas,u,b,tmp_vov,tmp_ovv)

  !$OMP SINGLE
  !$OMP TASK
  call gpu_copy(d_cc_space_f_oo, H_oo)
  !$OMP END TASK
  !$OMP END SINGLE

  call gpu_allocate(tmp_ovv, nO, nV, nV)
  call gpu_allocate(tmp_vov, nV, nO, nV)

  call gpu_blas_create(blas)

  !$OMP DO
  do u=1,nO
    call gpu_dgeam(blas, 'N', 'N', 1, nO*nV*nV, 1.d0, &
           tau_x%f(u,1,1,1), nO, 0.d0, tau_x%f(1,1,1,1), nO, tmp_ovv%f(1,1,1), 1)
    do b=1,nV
      call gpu_dgeam(blas, 'T', 'T', nV, nO, 1.d0, &
           tmp_ovv%f(1,1,b), nO, 0.d0, &
           tmp_ovv%f(1,1,b), nO, tmp_vov%f(1,1,b), nV)
    enddo
    call gpu_dgemm(blas, 'N','T',cholesky_mo_num,nV,nO*nV,1.d0, &
      d_cc_space_v_ov_chol%f(1,1,1), cholesky_mo_num, tmp_vov%f(1,1,1), nV, &
      0.d0, tau_kau%f(1,1,u), cholesky_mo_num)
  enddo
  !$OMP END DO

  call gpu_blas_destroy(blas)

  call gpu_deallocate(tmp_vov)
  call gpu_deallocate(tmp_ovv)

  !$OMP TASKWAIT
  !$OMP END PARALLEL

  call gpu_dgemm(blas_handle, 'T', 'N', nO, nO, cholesky_mo_num*nV, 1.d0, &
    tau_kau%f(1,1,1), cholesky_mo_num*nV,  d_cc_space_v_vo_chol%f(1,1,1), cholesky_mo_num*nV, &
    1.d0, H_oo%f(1,1), nO)

  call gpu_synchronize()
  call gpu_deallocate(tau_kau)
end

! H_vv

subroutine compute_H_vv_chol(nO,nV,tau_x,d_cc_space_f_vv, &
         d_cc_space_v_ov_chol,H_vv)
  use gpu
  implicit none

  integer, intent(in)              :: nO,nV
  type(gpu_double2), intent(in)    :: d_cc_space_f_vv
  type(gpu_double3), intent(in)    :: d_cc_space_v_ov_chol
  type(gpu_double4), intent(in)    :: tau_x
  type(gpu_double2), intent(out)   :: H_vv

  integer :: a,b,i,j,u,k, beta

  type(gpu_double3) :: tau_kia, tmp_oov

  call gpu_allocate(tau_kia, cholesky_mo_num, nO, nV)

  type(gpu_blas) :: blas

  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED) &
  !$OMP PRIVATE(a,b,tmp_oov,blas)

  !$OMP SINGLE
  !$OMP TASK
  call gpu_copy(d_cc_space_f_vv, H_vv)
  !$OMP END TASK
  !$OMP END SINGLE

  call gpu_blas_create(blas)
  call gpu_allocate(tmp_oov, nO, nO, nV)

  !$OMP DO
  do a = 1, nV
    do b=1,nV
      call gpu_dgeam(blas, 'N', 'N', nO, nO, 1.d0, &
        tau_x%f(1,1,a,b), nO, 0.d0, &
        tau_x%f(1,1,a,b), nO, tmp_oov%f(1,1,b), nO)
    enddo
    call gpu_dgemm(blas, 'N', 'T', cholesky_mo_num, nO, nO*nV, 1.d0, &
      d_cc_space_v_ov_chol%f(1,1,1), cholesky_mo_num, tmp_oov%f(1,1,1), nO, &
      0.d0, tau_kia%f(1,1,a), cholesky_mo_num)
  enddo
  !$OMP END DO

  call gpu_blas_destroy(blas)

  call gpu_deallocate(tmp_oov)
  !$OMP TASKWAIT
  !$OMP END PARALLEL

  call gpu_dgemm(blas_handle, 'T', 'N', nV, nV, cholesky_mo_num*nO, -1.d0, &
    tau_kia%f(1,1,1), cholesky_mo_num*nO,  d_cc_space_v_ov_chol%f(1,1,1), cholesky_mo_num*nO, &
    1.d0, H_vv%f(1,1), nV)

  call gpu_synchronize()
  call gpu_deallocate(tau_kia)
end

! H_vo
subroutine compute_H_vo_chol(nO,nV,t1,d_cc_space_f_vo, &
         d_cc_space_v_ov_chol,d_cc_space_v_vo_chol, H_vo)
  use gpu
  implicit none

  integer, intent(in)            :: nO,nV
  type(gpu_double2), intent(in)  :: t1, d_cc_space_f_vo
  type(gpu_double3), intent(in)  :: d_cc_space_v_ov_chol, d_cc_space_v_vo_chol
  type(gpu_double2), intent(out) :: H_vo

  integer :: a,b,i,j,u,k

  type(gpu_double1) :: tmp_k
  type(gpu_double3) :: tmp, tmp2

  call gpu_copy(d_cc_space_f_vo, H_vo)

  call gpu_allocate(tmp_k, cholesky_mo_num)

  call gpu_dgemm(blas_handle, 'N', 'N', cholesky_mo_num, 1, nO*nV, 2.d0, &
     d_cc_space_v_ov_chol%f(1,1,1), cholesky_mo_num, &
     t1%f(1,1), nO*nV, 0.d0, tmp_k%f(1), cholesky_mo_num)

  call gpu_dgemm(blas_handle, 'T', 'N', nV*nO, 1, cholesky_mo_num, 1.d0, &
      d_cc_space_v_vo_chol%f(1,1,1), cholesky_mo_num, tmp_k%f(1), cholesky_mo_num, 1.d0, &
      H_vo%f(1,1), nV*nO)

  call gpu_deallocate(tmp_k)


  call gpu_allocate(tmp,  cholesky_mo_num, nO, nO)

  call gpu_dgemm(blas_handle, 'N', 'T', cholesky_mo_num*nO, nO, nV, 1.d0, &
    d_cc_space_v_ov_chol%f(1,1,1), cholesky_mo_num*nO, t1%f(1,1), nO, 0.d0, tmp%f(1,1,1), cholesky_mo_num*nO)

  call gpu_allocate(tmp2, cholesky_mo_num, nO, nO)

  type(gpu_stream) :: stream(nO)
  do i=1,nO
    call gpu_stream_create(stream(i))
  enddo

  !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i,j)
  do i=1,nO
    do j=1,nO
      call gpu_set_stream(blas_handle,stream(j))
      call gpu_dgeam(blas_handle, 'N', 'N', cholesky_mo_num, 1, 1.d0, &
        tmp%f(1,i,j), cholesky_mo_num, 0.d0, &
        tmp%f(1,i,j), cholesky_mo_num, tmp2%f(1,j,i), cholesky_mo_num)
    enddo
  enddo
  !$OMP END PARALLEL DO

  call gpu_set_stream(blas_handle,gpu_default_stream)
  call gpu_synchronize()

  do i=1,nO
    call gpu_stream_destroy(stream(i))
  enddo
  call gpu_deallocate(tmp)

  call gpu_dgemm(blas_handle, 'T','N', nV, nO, cholesky_mo_num*nO, -1.d0, &
    d_cc_space_v_ov_chol%f(1,1,1), cholesky_mo_num*nO, tmp2%f(1,1,1), cholesky_mo_num*nO, &
    1.d0, H_vo%f(1,1), nV)

  call gpu_synchronize()
  call gpu_deallocate(tmp2)
end

! R1

subroutine compute_r1_space_chol(nO,nV,t1,t2,tau,H_oo,H_vv,H_vo,r1,max_r1,d_cc_space_f_ov,d_cc_space_f_vo, &
    d_cc_space_v_voov, d_cc_space_v_ovov, d_cc_space_v_oovo, d_cc_space_v_vo_chol, d_cc_space_v_vv_chol)
  use gpu
  implicit none

  ! in
  integer, intent(in)           :: nO, nV
  type(gpu_double2), intent(in) :: t1, H_oo, H_vo, H_vv, d_cc_space_f_ov,d_cc_space_f_vo
  type(gpu_double3), intent(in) :: d_cc_space_v_vo_chol, d_cc_space_v_vv_chol
  type(gpu_double4), intent(in) :: t2, tau, d_cc_space_v_voov, d_cc_space_v_ovov, d_cc_space_v_oovo

  ! out
  type(gpu_double2), intent(out) :: r1
  double precision, intent(out)  :: max_r1

  ! internal
  integer                       :: u,i,j,beta,a,b

  type(gpu_stream) :: stream(nV)

  do a=1,nV
    call gpu_stream_create(stream(a))
  enddo

  type(gpu_double2) :: X_oo
  call gpu_allocate(X_oo,nO,nO)

  call gpu_copy(d_cc_space_f_ov, r1)

  call gpu_set_stream(blas_handle, stream(1))
  call gpu_dgemm(blas_handle, 'N','N', nO, nV, nV, &
             1d0, t1%f(1,1)  , size(t1%f,1), &
                  H_vv%f(1,1), size(H_vv%f,1), &
             1d0, r1%f(1,1)  , size(r1%f,1))

  call gpu_dgemm(blas_handle, 'N','N', nO, nV, nO, &
             -1d0, H_oo%f(1,1), size(H_oo%f,1), &
                   t1%f(1,1)  , size(t1%f,1), &
              1d0, r1%f(1,1), size(r1%f,1))

  call gpu_set_stream(blas_handle, stream(nV))
  call gpu_dgemm(blas_handle, 'N','N', nO, nO, nV, &
             -2d0, t1%f(1,1), size(t1%f,1), &
                   d_cc_space_f_vo%f(1,1), size(d_cc_space_f_vo%f,1), &
              0d0, X_oo%f(1,1), size(X_oo%f,1))

  call gpu_synchronize()
  call gpu_set_stream(blas_handle, gpu_default_stream)

  call gpu_dgemm(blas_handle, 'T','N', nO, nV, nO, &
             1d0, X_oo%f(1,1), size(X_oo%f,2), &
                  t1%f(1,1)  , size(t1%f,1), &
             1d0, r1%f(1,1)  , size(r1%f,1))



  type(gpu_double4) :: X_voov
  call gpu_allocate(X_voov, nV, nO, nO, nV)

  do i=1,nO
    do beta=1,nV
      call gpu_set_stream(blas_handle, stream(beta))
      call gpu_dgeam(blas_handle, 'T', 'T', nV, nO,  -1.d0, t2%f(1,i,1,beta), &
         nO*nO, t1%f(i,beta), t1%f(1,1), nO, X_voov%f(1,i,1,beta), nV*nO)
    enddo
  enddo

  do beta=1,nV
    call gpu_set_stream(blas_handle, stream(beta))
    call gpu_dgeam(blas_handle, 'N', 'T', nV, nO*nO,  1.d0, X_voov%f(1,1,1,beta), &
         nV, 2.d0, t2%f(1,1,1,beta), nO*nO, X_voov%f(1,1,1,beta), nV)
  enddo

  call gpu_synchronize()
  call gpu_deallocate(X_oo)

  call gpu_set_stream(blas_handle, gpu_default_stream)

  call gpu_dgemv(blas_handle, 'T', nV*nO, nO*nV, &
             1d0, X_voov%f(1,1,1,1), size(X_voov%f,1) * size(X_voov%f,2), &
                  H_vo%f(1,1)  , 1, &
             1d0, r1%f(1,1)    , 1)

  type(gpu_double4) :: X_ovov
  call gpu_allocate(X_ovov, nO, nV, nO, nV)

  do beta = 1, nV
    call gpu_set_stream(blas_handle, stream(beta))
    do u=1,nO
      call gpu_dgeam(blas_handle, 'N', 'T', nO, nV, -1.d0, d_cc_space_v_ovov%f(1,1,u,beta), &
      nO, 2.d0, d_cc_space_v_voov%f(1,u,1,beta), nV*nO, X_ovov%f(1,1,u,beta), nO)
    enddo
  enddo

  call gpu_set_stream(blas_handle, gpu_default_stream)
  call gpu_synchronize()
  call gpu_deallocate(X_voov)

  call gpu_dgemv(blas_handle, 'T', nO*nV, nO*nV, &
             1d0, X_ovov%f(1,1,1,1), size(X_ovov%f,1) * size(X_ovov%f,2), &
                  t1%f(1,1), 1, &
             1d0, r1%f(1,1), 1)


  integer :: iblock, block_size, nVmax
  type(gpu_double4) :: W_vvov, W_vvov_tmp, T_vvoo

  block_size = 16
  call gpu_allocate(T_vvoo, nV,nV,nO,nO)

  call gpu_dgeam(blas_handle, 'T', 'N', nV*nV, nO*nO, 1.d0, tau%f(1,1,1,1), &
    nO*nO, 0.d0, T_vvoo%f(1,1,1,1), nV*nV, T_vvoo%f(1,1,1,1), nV*nV)

  call gpu_allocate(W_vvov,nV, nV,nO,block_size)
  call gpu_allocate(W_vvov_tmp, nV,nO,nV,block_size)

  do iblock = 1, nV, block_size
    nVmax = min(block_size,nV-iblock+1)

    call gpu_dgemm(blas_handle, 'T','N', nV*nO, nV*nVmax, cholesky_mo_num, 1.d0, &
      d_cc_space_v_vo_chol%f(1,1,1) , cholesky_mo_num, &
      d_cc_space_v_vv_chol%f(1,1,iblock), cholesky_mo_num, &
      0.d0, W_vvov_tmp%f(1,1,1,1), nV*nO)

    call gpu_synchronize()
    do b=1,nV
      call gpu_set_stream(blas_handle, stream(b))
      do i=1,nO
        call gpu_dgeam(blas_handle, 'N', 'N', nV, nVmax,  2.d0, W_vvov_tmp%f(1,i,b,1), &
         nV*nO*nV, 0.d0, W_vvov_tmp%f(1,i,b,1), nV*nO*nV, W_vvov%f(1,b,i,1), nV*nV*nO)
      enddo
    enddo

    call gpu_synchronize()

    do beta = 1,  nVmax
      call gpu_set_stream(blas_handle, stream(beta))
      call gpu_dgeam(blas_handle, 'N', 'T', nV, nV*nO,  1.d0, W_vvov%f(1,1,1,beta), &
         nV, -1.d0, W_vvov_tmp%f(1,1,1,beta), nV*nO, W_vvov%f(1,1,1,beta), nV)
    enddo
    call gpu_synchronize()

    call gpu_dgemm(blas_handle, 'T','N',nO,nVmax,nO*nV*nV, &
             1d0, T_vvoo%f(1,1,1,1), nV*nV*nO, &
                  W_vvov%f(1,1,1,1), nO*nV*nV, &
             1d0, r1%f(1,iblock), nO)
  enddo

  call gpu_deallocate(X_ovov)

  type(gpu_double4) :: W_oovo
  call gpu_allocate(W_oovo, nO,nO,nV,nO)

  do u = 1, nO
    do a = 1, nV
      call gpu_set_stream(blas_handle, stream(a))
      call gpu_dgeam(blas_handle, 'N', 'T', nO, nO,  2.d0, d_cc_space_v_oovo%f(1,1,a,u), &
         nO, -1.d0, d_cc_space_v_oovo%f(1,1,a,u), nO, W_oovo%f(1,1,a,u), nO)
    enddo
  enddo

  call gpu_set_stream(blas_handle, gpu_default_stream)
  call gpu_synchronize()

  call gpu_deallocate(W_vvov)
  call gpu_deallocate(T_vvoo)

  ! Change the sign for consistency with the code in spin orbitals
  call gpu_dgemm(blas_handle, 'T','N', nO, nV, nO*nO*nV, &
              1d0, W_oovo%f(1,1,1,1), size(W_oovo%f,1) * size(W_oovo%f,2) * size(W_oovo%f,3), &
                   tau%f(1,1,1,1), size(tau%f,1) * size(tau%f,2) * size(tau%f,3), &
             -1d0, r1%f(1,1), size(r1%f,1))

  call gpu_synchronize()
  call gpu_deallocate(W_oovo)

  max_r1 = 0d0
  do a = 1, nV
    do i = 1, nO
      max_r1 = max(dabs(r1%f(i,a)), max_r1)
    enddo
  enddo

  do a=1,nV
    call gpu_stream_destroy(stream(a))
  enddo

end


! R2

subroutine compute_r2_space_chol(nO,nV,t1,t2,tau,H_oo,H_vv, &
    d_cc_space_v_oovv, d_cc_space_v_vooo, d_cc_space_v_oooo, d_cc_space_v_oovo, d_cc_space_v_ovvo, d_cc_space_v_ovoo, &
    d_cc_space_v_ovov, d_cc_space_v_vvoo, d_cc_space_v_oo_chol, d_cc_space_v_ov_chol, d_cc_space_v_vo_chol, d_cc_space_v_vv_chol, &
    d_cc_space_f_vo, &
    r2,max_r2)
  use gpu
  implicit none

  ! in
  integer, intent(in)            :: nO, nV
  type(gpu_double2), intent(in)  :: t1, H_oo, H_vv, d_cc_space_f_vo
  type(gpu_double4), intent(in)  :: t2, tau, d_cc_space_v_oovv
  type(gpu_double4), intent(in)  :: d_cc_space_v_vooo, d_cc_space_v_oooo
  type(gpu_double4), intent(in)  :: d_cc_space_v_vvoo, d_cc_space_v_oovo
  type(gpu_double4), intent(in)  :: d_cc_space_v_ovvo, d_cc_space_v_ovoo
  type(gpu_double4), intent(in)  :: d_cc_space_v_ovov
  type(gpu_double3), intent(in)  :: d_cc_space_v_oo_chol, d_cc_space_v_ov_chol
  type(gpu_double3), intent(in)  :: d_cc_space_v_vo_chol, d_cc_space_v_vv_chol

  ! out
  double precision, intent(out)  :: max_r2
  type(gpu_double4), intent(out) :: r2

  ! internal
  integer                       :: u,v,i,j,beta,gam,a,b
  double precision              :: max_r2_local

  type(gpu_stream) :: stream(nV)

  call set_multiple_levels_omp(.False.)

  call gpu_copy(d_cc_space_v_oovv, r2)

  type(gpu_double4) :: A1
  call gpu_allocate(A1,nO,nO,nO,nO)
  call compute_A1_chol(nO,nV,t1,t2,tau,d_cc_space_v_vooo, &
      d_cc_space_v_oooo, d_cc_space_v_vvoo, A1)

  call gpu_dgemm(blas_handle, 'N','N',nO*nO,nV*nV,nO*nO, &
             1d0, A1%f(1,1,1,1), size(A1%f,1) * size(A1%f,2), &
                  tau%f(1,1,1,1), size(tau%f,1) * size(tau%f,2), &
             1d0, r2%f(1,1,1,1), size(r2%f,1) * size(r2%f,2))

  call gpu_deallocate(A1)

  integer :: block_size, iblock, k
  block_size = 16
  type(gpu_double3) :: tmp_cc, B1, tmpB1
  type(gpu_double2) :: tmp_cc2

  call gpu_allocate(tmp_cc,cholesky_mo_num,nV,nV)
  call gpu_dgemm(blas_handle, 'N','N', cholesky_mo_num*nV, nV, nO, 1.d0, &
      d_cc_space_v_vo_chol%f(1,1,1), cholesky_mo_num*nV, t1%f(1,1), nO, 0.d0, tmp_cc%f(1,1,1), cholesky_mo_num*nV)

  call set_multiple_levels_omp(.False.)
  call gpu_synchronize()

  type(gpu_blas) :: blas

  !$OMP PARALLEL PRIVATE(gam, iblock, B1, tmpB1, tmp_cc2, beta, b, a, blas)
  call gpu_allocate(B1,nV,nV,block_size)
  call gpu_allocate(tmpB1,nV,block_size,nV)
  call gpu_allocate(tmp_cc2,cholesky_mo_num,nV)

  call gpu_blas_create(blas)

  !$OMP DO
  do gam = 1, nV

    call gpu_dgeam(blas, 'N', 'N', cholesky_mo_num, nV, 1.d0, d_cc_space_v_vv_chol%f(1,1,gam), &
         cholesky_mo_num, -1.d0, tmp_cc%f(1,1,gam), cholesky_mo_num, tmp_cc2%f(1,1), cholesky_mo_num)

    do iblock = 1, nV, block_size

        call gpu_dgemm(blas, 'T', 'N', nV*min(block_size, nV-iblock+1), nV, cholesky_mo_num, &
                -1.d0, tmp_cc%f(1,1,iblock), cholesky_mo_num, &
                d_cc_space_v_vv_chol%f(1,1,gam), cholesky_mo_num, &
                0.d0, tmpB1%f(1,1,1), nV*block_size)

        call gpu_dgemm(blas, 'T','N', nV*min(block_size, nV-iblock+1), nV, cholesky_mo_num, &
                1.d0, d_cc_space_v_vv_chol%f(1,1,iblock), cholesky_mo_num, &
                tmp_cc2%f(1,1), cholesky_mo_num, &
                1.d0, tmpB1%f(1,1,1), nV*block_size)

        do beta = iblock, min(nV, iblock+block_size-1)
          call gpu_dgeam(blas, 'N', 'N', nV, nV, 1.d0, tmpB1%f(1,beta-iblock+1,1), &
             nV*block_size, 0.d0, B1%f(1,1,beta-iblock+1), nV, B1%f(1,1,beta-iblock+1), nV)
        enddo

        call gpu_dgemm(blas, 'N','N',nO*nO,min(block_size, nV-iblock+1),nV*nV, &
              1d0, tau%f(1,1,1,1), size(tau%f,1) * size(tau%f,2), &
                   B1%f(1,1,1) , size(B1%f ,1) * size(B1%f ,2), &
              1d0, r2%f(1,1,iblock,gam),  size(r2%f ,1) * size(r2%f ,2))
      enddo

  enddo
  !$OMP ENDDO

  call gpu_blas_destroy(blas)

  call gpu_deallocate(B1)
  call gpu_deallocate(tmpB1)
  call gpu_deallocate(tmp_cc2)
  !$OMP END PARALLEL

  call gpu_deallocate(tmp_cc)

  type(gpu_double4) :: X_oovv
  call gpu_allocate(X_oovv,nO,nO,nV,nV)
  call gpu_copy(t2,X_oovv)

  type(gpu_double2) :: g_occ, g_vir
  call gpu_allocate(g_vir,nV,nV)
  call gpu_allocate(g_occ,nO,nO)
  call compute_g_vir_chol(nO,nV,t1,t2,H_vv,d_cc_space_f_vo, &
    d_cc_space_v_ov_chol, d_cc_space_v_vv_chol, g_vir)
  call compute_g_occ_chol(nO,nV,t1%f,t2%f,H_oo%f,g_occ%f)

  type(gpu_double4) :: Y_oovv
  call gpu_allocate(Y_oovv,nO,nO,nV,nV)

  call gpu_dgemm(blas_handle, 'N','N',nO*nO*nV,nV,nV, &
             1d0, X_oovv%f(1,1,1,1), size(X_oovv%f,1) * size(X_oovv%f,2) * size(X_oovv%f,3), &
                  g_vir%f(1,1), size(g_vir%f,1), &
             0d0, Y_oovv%f(1,1,1,1), size(Y_oovv%f,1) * size(Y_oovv%f,2) * size(Y_oovv%f,3))

  call gpu_dgemm(blas_handle, 'N','N',nO,nO*nV*nV,nO, &
             -1d0, g_occ%f(1,1), size(g_occ%f,1), &
                  t2%f(1,1,1,1)    , size(t2%f,1), &
             1d0, Y_oovv%f(1,1,1,1), size(Y_oovv%f,1))

  call gpu_dgemm(blas_handle, 'N','N',nO*nO*nV,nV,nO, &
            -1d0, d_cc_space_v_oovo%f(1,1,1,1), size(cc_space_v_oovo,1) * size(cc_space_v_oovo,2) * size(cc_space_v_oovo,3), &
                  t1%f(1,1) , size(t1%f,1), &
             1d0, Y_oovv%f(1,1,1,1), size(Y_oovv%f,1) * size(Y_oovv%f,2) * size(Y_oovv%f,3))


  call gpu_dgeam(blas_handle, 'N', 'N', nO*nO, nV*nV, 1.d0, Y_oovv%f(1,1,1,1), &
         nO*nO, 1.d0, r2%f(1,1,1,1), nO*nO, r2%f(1,1,1,1), nO*nO)

  call gpu_synchronize()
  call gpu_deallocate(X_oovv)

  call gpu_deallocate(g_vir)
  call gpu_deallocate(g_occ)

  type(gpu_double4) :: X_vovo, Y_oovo
  call gpu_allocate(X_vovo,nV,nO,nV,nO)

  do a=1,nV
    call gpu_stream_create(stream(a))
  enddo

  do gam = 1, nV
    call gpu_set_stream(blas_handle, stream(gam))
    do beta = 1, nV
      call gpu_dgeam(blas_handle, 'N', 'T', nO, nO, 1.d0, r2%f(1,1,beta,gam), &
           nO, 1.d0, Y_oovv%f(1,1,gam,beta), nO, r2%f(1,1,beta,gam), nO)
    enddo
  enddo

  do i = 1, nO
    do gam = 1, nV
      call gpu_set_stream(blas_handle, stream(gam))
      call gpu_dgeam(blas_handle, 'T', 'N', nV, nO, 1.d0, d_cc_space_v_ovvo%f(1,1,gam,i), &
           nO, 0.d0, X_vovo%f(1,1,gam,i), nV, X_vovo%f(1,1,gam,i), nV)
    enddo
  enddo

  do a=1,nV
    call gpu_stream_destroy(stream(a))
  enddo
  call gpu_set_stream(blas_handle, gpu_default_stream)



  call gpu_allocate(Y_oovo,nO,nO,nV,nO)

  !$OMP PARALLEL PRIVATE(blas, iblock, gam, X_vovv)
  call gpu_blas_create(blas)
  type(gpu_double4) :: X_vovv
  call gpu_allocate(X_vovv,nV,nO,nV,block_size)
  !$OMP DO
  do iblock = 1, nV, block_size
    do gam = iblock, min(nV, iblock+block_size-1)
      call gpu_dgemm(blas, 'T','N',nV, nO*nV, cholesky_mo_num, 1.d0, &
        d_cc_space_v_vv_chol%f(1,1,gam), cholesky_mo_num, d_cc_space_v_ov_chol%f(1,1,1), &
        cholesky_mo_num, 0.d0, X_vovv%f(1,1,1,gam-iblock+1), nV)

    enddo

    call gpu_dgemm(blas, 'N','N', nO, &
             nO*nV*min(block_size, nV-iblock+1),nV, &
             1.d0, t1%f(1,1)    , size(t1%f,1), &
             X_vovv%f(1,1,1,1), size(X_vovv%f,1), &
             0d0, Y_oovv%f(1,1,1,iblock), size(Y_oovv%f,1))
  enddo
  !$OMP END DO

  call gpu_blas_destroy(blas)
  call gpu_deallocate(X_vovv)
  !$OMP END PARALLEL

  call gpu_dgemm(blas_handle, 'N','N',nO,nO*nV*nO,nV, &
             1d0, t1%f(1,1), size(t1%f,1), &
                  X_vovo%f(1,1,1,1), size(X_vovo%f,1), &
             0d0, Y_oovo%f(1,1,1,1), size(Y_oovo%f,1))

  call gpu_dgemm(blas_handle, 'N','N',nO*nO*nV, nV, nO, &
            -1d0, Y_oovo%f(1,1,1,1), size(Y_oovo%f,1) * size(Y_oovo%f,2) * size(Y_oovo%f,3), &
                  t1%f(1,1)    , size(t1%f,1), &
             1d0, Y_oovv%f(1,1,1,1), size(Y_oovv%f,1) * size(Y_oovv%f,2) * size(Y_oovv%f,3))

  call gpu_synchronize()
  call gpu_deallocate(X_vovo)
  call gpu_deallocate(Y_oovo)

  do a=1,nV
    call gpu_stream_create(stream(a))
  enddo

  do gam = 1, nV
    call gpu_set_stream(blas_handle, stream(gam))
    do beta = 1, nV
      call gpu_dgeam(blas_handle, 'T', 'N', nO, nO, 1.d0, Y_oovv%f(1,1,beta,gam), &
           nO, 1.d0, r2%f(1,1,beta,gam), nO, r2%f(1,1,beta,gam), nO)
    enddo
    do j=1,nO
      call gpu_dgeam(blas_handle, 'N', 'N', nO, nV, 1.d0, r2%f(1,j,1,gam), &
           nO*nO, 1.d0, Y_oovv%f(1,j,gam,1), nO*nO*nV, r2%f(1,j,1,gam), nO*nO)
    enddo
  enddo

  call gpu_set_stream(blas_handle, gpu_default_stream)


  call gpu_synchronize()
  call gpu_deallocate(Y_oovv)

  type(gpu_double4) :: X_ovvo
  type(gpu_double3) :: tcc, tcc2
  call gpu_allocate(tcc2,cholesky_mo_num,nV,nO)
  call gpu_allocate(X_ovvo,nO,nV,nV,nO)
  call gpu_allocate(tcc,cholesky_mo_num,nO,nV)

  call gpu_dgemm(blas_handle, 'N','T', cholesky_mo_num*nV, nO, nV, 1.d0, &
     d_cc_space_v_vv_chol%f(1,1,1), cholesky_mo_num*nV, t1%f(1,1), nO, &
     0.d0, tcc2%f(1,1,1), cholesky_mo_num*nV)

  call gpu_dgemm(blas_handle, 'N','N', cholesky_mo_num*nO, nV, nO, 1.d0, &
     d_cc_space_v_oo_chol%f(1,1,1), cholesky_mo_num*nO, t1%f(1,1), nO, &
     0.d0, tcc%f(1,1,1), cholesky_mo_num*nO)

  call gpu_dgemm(blas_handle, 'T','N', nO*nV, nV*nO, cholesky_mo_num, 1.d0, &
              tcc%f(1,1,1), cholesky_mo_num, tcc2%f(1,1,1), cholesky_mo_num, 0.d0, &
              X_ovvo%f(1,1,1,1), nO*nV)

  call gpu_synchronize()


  do gam = 1, nV
    call gpu_set_stream(blas_handle, stream(gam))
    do j=1,nO
      call gpu_dgeam(blas_handle, 'N', 'N', nO, nV, -1.d0, X_ovvo%f(1,1,gam,j), &
           nO, 1.d0, r2%f(1,j,1,gam), nO*nO, r2%f(1,j,1,gam), nO*nO)
    enddo
    do beta = 1, nV
      call gpu_dgeam(blas_handle, 'T', 'N', nO, nO, -1.d0, X_ovvo%f(1,gam,beta,1), &
           nO*nV*nV, 1.d0, r2%f(1,1,beta,gam), nO, r2%f(1,1,beta,gam), nO)
    enddo
  enddo

  call gpu_set_stream(blas_handle, gpu_default_stream)

  call gpu_synchronize
  call gpu_deallocate(tcc)
  call gpu_deallocate(tcc2)
  call gpu_deallocate(X_ovvo)


  type(gpu_double4) :: J1, K1
  type(gpu_double4) :: Y_voov, Z_ovov


  call gpu_allocate(J1,nO,nV,nV,nO)
  call compute_J1_chol(nO,nV,t1,t2,d_cc_space_v_ovvo,d_cc_space_v_ovoo, &
       d_cc_space_v_vvoo,d_cc_space_v_vo_chol,d_cc_space_v_vv_chol,J1)

  call gpu_allocate(K1,nO,nV,nO,nV)
  call compute_K1_chol(nO,nV,t1,t2,d_cc_space_v_ovoo,d_cc_space_v_vvoo, &
       d_cc_space_v_ovov,d_cc_space_v_ov_chol,d_cc_space_v_vv_chol,K1)


  call gpu_allocate(X_ovvo,nO,nV,nV,nO)
  call gpu_allocate(Y_voov,nV,nO,nO,nV)

  do a=1, nV
    call gpu_set_stream(blas_handle, stream(a))
    do i=1, nO
      call gpu_dgeam(blas_handle, 'N', 'N', nO, nV, 1.d0, J1%f(1,a,1,i), &
           nO*nV, -0.5d0, K1%f(1,a,i,1), nO*nV*nO, X_ovvo%f(1,1,a,i), nO)
      call gpu_dgeam(blas_handle, 'T', 'T', nV, nO, 2.d0, t2%f(1,i,1,a), &
           nO*nO, -1.d0, t2%f(1,i,a,1), nO*nO*nV, Y_voov%f(1,1,i,a), nV)
    enddo
  enddo

  call gpu_allocate(Z_ovov,nO,nV,nO,nV)

  call gpu_synchronize()
  call gpu_deallocate(J1)
  call gpu_set_stream(blas_handle, gpu_default_stream)


  call gpu_dgemm(blas_handle, 'N','N', nO*nV,nO*nV,nV*nO, &
             1d0, X_ovvo%f(1,1,1,1), size(X_ovvo%f,1) * size(X_ovvo%f,2), &
                  Y_voov%f(1,1,1,1), size(Y_voov%f,1) * size(Y_voov%f,2), &
             0d0, Z_ovov%f(1,1,1,1), size(Z_ovov%f,1) * size(Z_ovov%f,2))

  call gpu_synchronize()
  call gpu_deallocate(Y_voov)
  call gpu_deallocate(X_ovvo)

  type(gpu_double4) :: Y_ovov, X_ovov
  call gpu_allocate(X_ovov,nO,nV,nO,nV)
  call gpu_allocate(Y_ovov,nO,nV,nO,nV)

  do a=1, nV
    call gpu_set_stream(blas_handle, stream(a))
    do j=1,nO
      call gpu_dgeam(blas_handle, 'N', 'N', nO, nV, 1.d0, t2%f(1,j,1,a), &
           nO*nO, 0.d0, t2%f(1,j,1,a), nO*nO, Y_ovov%f(1,a,j,1), nO*nV*nO)
    enddo
    do beta=1, nV
      call gpu_dgeam(blas_handle, 'T', 'T', nO, nO, 0.5d0, K1%f(1,a,1,beta), &
           nO*nV, 0.d0, K1%f(1,a,1,beta), nO*nV, X_ovov%f(1,a,1,beta), nO*nV)
    enddo
  enddo
  call gpu_set_stream(blas_handle, gpu_default_stream)

  call gpu_synchronize()

  call gpu_dgemm(blas_handle, 'T','N',nO*nV,nO*nV,nO*nV, &
            -1d0, X_ovov%f(1,1,1,1), size(X_ovov%f,1) * size(X_ovov%f,2), &
                  Y_ovov%f(1,1,1,1), size(Y_ovov%f,1) * size(Y_ovov%f,2), &
             1d0, Z_ovov%f(1,1,1,1), size(Z_ovov%f,1) * size(Z_ovov%f,2))

  call gpu_synchronize()

  do gam=1, nV
    call gpu_set_stream(blas_handle, stream(gam))
    do j=1,nO
      call gpu_dgeam(blas_handle, 'N', 'N', nO, nV, 1.d0, r2%f(1,j,1,gam), &
           nO*nO, 1.d0, Z_ovov%f(1,1,j,gam), nO, r2%f(1,j,1,gam), nO*nO)
      call gpu_dgeam(blas_handle, 'N', 'N', nO, nV, 1.d0, K1%f(1,1,j,gam), &
           nO, 0.d0, K1%f(1,1,j,gam), nO, X_ovov%f(1,gam,j,1), nO*nV*nO)
      call gpu_dgeam(blas_handle, 'N', 'N', nO, nO, 1.d0, t2%f(1,j,1,gam), &
           nO*nO, 0.d0, t2%f(1,j,1,gam), nO*nO, Y_ovov%f(1,gam,j,1), nO*nV*nO)
    enddo
    do beta=1, nV
      call gpu_dgeam(blas_handle, 'N', 'T', nO, nO, 1.d0, r2%f(1,1,beta,gam), &
           nO, 1.d0, Z_ovov%f(1,gam,1,beta), nO*nV, r2%f(1,1,beta,gam), nO)
    enddo
  enddo

  call gpu_set_stream(blas_handle, gpu_default_stream)

  call gpu_deallocate(K1)

  call gpu_dgemm(blas_handle, 'N','N',nO*nV,nO*nV,nO*nV, &
             1d0, X_ovov%f(1,1,1,1), size(X_ovov%f,1) * size(X_ovov%f,2), &
                  Y_ovov%f(1,1,1,1), size(Y_ovov%f,1) * size(Y_ovov%f,2), &
             0d0, Z_ovov%f(1,1,1,1), size(Z_ovov%f,1) * size(Z_ovov%f,2))

  call gpu_synchronize()

  call gpu_deallocate(X_ovov)
  call gpu_deallocate(Y_ovov)

  ! Change the sign for consistency with the code in spin orbitals
  do gam = 1, nV
    call gpu_set_stream(blas_handle, stream(gam))
    do j=1,nO
      call gpu_dgeam(blas_handle, 'N', 'N', nO, nV,  1.d0, r2%f(1,j,1,gam), &
           nO*nO, -1.d0, Z_ovov%f(1,gam,j,1), nO*nV*nO, r2%f(1,j,1,gam), nO*nO)
    enddo
    do beta = 1, nV
      call gpu_dgeam(blas_handle, 'N', 'T', nO, nO, -1.d0, r2%f(1,1,beta,gam), &
           nO, 1.d0, Z_ovov%f(1,beta,1,gam), nO*nV, r2%f(1,1,beta,gam), nO)
    enddo
  enddo

  call gpu_deallocate(Z_ovov)

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
          max_r2_local = max(r2%f(i,j,a,b), max_r2_local)
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

subroutine compute_A1_chol(nO,nV,t1,t2,tau,d_cc_space_v_vooo, &
  d_cc_space_v_oooo, d_cc_space_v_vvoo, A1)
  use gpu
  implicit none

  integer, intent(in)            :: nO,nV
  type(gpu_double2), intent(in)  :: t1
  type(gpu_double4), intent(in)  :: t2, tau
  type(gpu_double4), intent(in)  :: d_cc_space_v_vooo, d_cc_space_v_oooo, d_cc_space_v_vvoo
  type(gpu_double4), intent(out) :: A1

  integer :: a,tmp_a,b,k,l,c,d,tmp_c,tmp_d,i,j,u,v, beta

  type(gpu_double4) :: Y_oooo
  call gpu_allocate(Y_oooo,nO,nO,nO,nO)

  ! A1(u,v,i,j) = cc_space_v_oooo(u,v,i,j)
  ! A1(u,v,i,j) += cc_space_v_ovoo(u,a,i,j) * t1(v,a) &

  call gpu_dgemm(blas_handle, 'N','N', nO, nO*nO*nO, nV, &
             1d0, t1%f(1,1)  , size(t1%f,1), &
                  d_cc_space_v_vooo%f(1,1,1,1), size(d_cc_space_v_vooo%f,1), &
             0d0, Y_oooo%f(1,1,1,1), size(Y_oooo%f,1))

  type(gpu_stream) :: stream(nO)

  do i=1, nO
    call gpu_stream_create(stream(i))
  enddo

  call gpu_synchronize()

  do j = 1, nO
    call gpu_set_stream(blas_handle, stream(j))
    do i = 1, nO
      call gpu_dgeam(blas_handle, 'N', 'T', nO, nO, 1.d0, d_cc_space_v_oooo%f(1,1,i,j), &
           nO, 1.d0, Y_oooo%f(1,1,j,i), nO, A1%f(1,1,i,j), nO)
    enddo
    call gpu_dgeam(blas_handle, 'N', 'N', nO, nO*nO, 1.d0, A1%f(1,1,1,j), &
         nO, 1.d0, Y_oooo%f(1,1,1,j), nO, A1%f(1,1,1,j), nO)
  enddo

  call gpu_set_stream(blas_handle, gpu_default_stream)
  do i=1, nO
    call gpu_stream_destroy(stream(i))
  enddo

  call gpu_deallocate(Y_oooo)

  ! A1(u,v,i,j) += cc_space_v_vvoo(a,b,i,j) * tau(u,v,a,b)
  call gpu_dgemm(blas_handle, 'N','N', nO*nO, nO*nO, nV*nV, &
             1d0, tau%f(1,1,1,1), size(tau%f,1) * size(tau%f,2), &
                  d_cc_space_v_vvoo%f(1,1,1,1), size(d_cc_space_v_vvoo%f,1) * size(d_cc_space_v_vvoo%f,2), &
             1d0, A1%f(1,1,1,1), size(A1%f,1) * size(A1%f,2))
  call gpu_synchronize()

end

! g_occ

subroutine compute_g_occ_chol(nO,nV,t1,t2,H_oo,g_occ)
  use gpu

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

subroutine compute_g_vir_chol(nO,nV,t1,t2,H_vv,d_cc_space_f_vo, &
  d_cc_space_v_ov_chol, d_cc_space_v_vv_chol, g_vir)
  use gpu

  implicit none

  integer, intent(in)           :: nO,nV
  type(gpu_double2), intent(in)  :: t1, H_vv, d_cc_space_f_vo
  type(gpu_double3), intent(in)  :: d_cc_space_v_ov_chol, d_cc_space_v_vv_chol
  type(gpu_double4), intent(in)  :: t2
  type(gpu_double2), intent(out) :: g_vir

  integer :: a,tmp_a,b,k,l,c,d,tmp_c,tmp_d,i,j,u,v, beta, gam

  type(gpu_stream) :: stream(max(nO,4))

  do i=1,max(nO,4)
    call gpu_stream_create(stream(i))
  enddo

  call gpu_set_stream(blas_handle, stream(1))
  call gpu_dgemm(blas_handle, 'N','N',nV,nV,nO, &
             -1d0, d_cc_space_f_vo%f(1,1) , size(d_cc_space_f_vo%f,1), &
                   t1%f(1,1)   , size(t1%f,1), &
              0d0, g_vir%f(1,1), size(g_vir%f,1))

  type(gpu_double1) :: tmp_k
  type(gpu_double3) :: tmp_vo, tmp_vo2

  call gpu_allocate(tmp_k,cholesky_mo_num)

  call gpu_set_stream(blas_handle, stream(2))
  call gpu_dgemm(blas_handle, 'N','N', cholesky_mo_num, 1, nO*nV, 1.d0, &
    d_cc_space_v_ov_chol%f(1,1,1), cholesky_mo_num, t1%f(1,1), nO*nV, 0.d0, tmp_k%f(1), cholesky_mo_num)

  call gpu_dgemm(blas_handle, 'T','N', nV*nV, 1, cholesky_mo_num, 2.d0, &
    d_cc_space_v_vv_chol%f(1,1,1), cholesky_mo_num, tmp_k%f(1), cholesky_mo_num, 1.d0, &
    g_vir%f(1,1), nV*nV)

  call gpu_set_stream(blas_handle, stream(3))
  call gpu_allocate(tmp_vo,cholesky_mo_num,nV,nO)

  call gpu_dgemm(blas_handle, 'N','T',cholesky_mo_num*nV, nO, nV, 1.d0, &
    d_cc_space_v_vv_chol%f(1,1,1), cholesky_mo_num*nV, t1%f(1,1), nO, 0.d0, tmp_vo%f(1,1,1), cholesky_mo_num*nV)

  call gpu_allocate(tmp_vo2,cholesky_mo_num,nO,nV)

  call gpu_synchronize()
  call gpu_deallocate(tmp_k)

  do i=1,nO
    call gpu_set_stream(blas_handle, stream(i))
    call gpu_dgeam(blas_handle, 'N', 'N', cholesky_mo_num, nV, -1.d0, tmp_vo%f(1,1,i), &
         cholesky_mo_num, 0.d0, tmp_vo%f(1,1,i), cholesky_mo_num, tmp_vo2%f(1,i,1), cholesky_mo_num*nO)
  enddo

  call gpu_set_stream(blas_handle, gpu_default_stream)

  do i=1,max(nO,4)
    call gpu_stream_destroy(stream(i))
  enddo
  call gpu_deallocate(tmp_vo)

  call gpu_dgeam(blas_handle, 'N', 'N', nV, nV, 1.d0, g_vir%f(1,1), &
         nV, 1.d0, H_vv%f(1,1), nV, g_vir%f(1,1), nV)

  call gpu_dgemm(blas_handle, 'T','N', nV, nV, nO*cholesky_mo_num, 1.d0, &
     d_cc_space_v_ov_chol%f(1,1,1), cholesky_mo_num*nO, &
     tmp_vo2%f(1,1,1), cholesky_mo_num*nO, 1.d0, g_vir%f(1,1), nV)

  call gpu_synchronize()
  call gpu_deallocate(tmp_vo2)

end

! J1

subroutine compute_J1_chol(nO,nV,t1,t2,v_ovvo,v_ovoo,v_vvoo,d_cc_space_v_vo_chol,d_cc_space_v_vv_chol,J1)
  use gpu
  implicit none

  integer, intent(in)            :: nO,nV
  type(gpu_double2), intent(in)  :: t1
  type(gpu_double4), intent(in)  :: t2, v_ovvo, v_ovoo, v_vvoo
  type(gpu_double4), intent(out) :: J1
  type(gpu_double3), intent(out) :: d_cc_space_v_vo_chol,d_cc_space_v_vv_chol

  integer :: a,tmp_a,b,k,l,c,d,tmp_c,tmp_d,i,j,u,v, beta, gam

  type(gpu_double4) :: X_ovoo, Y_ovov

  call gpu_allocate(X_ovoo,nO,nV,nO,nO)

  type(gpu_stream) :: stream(nV)


  do i=1,nO
    call gpu_stream_create(stream(i))
  enddo

  do j = 1, nO
    call gpu_set_stream(blas_handle, stream(j))
    do i = 1, nO
      call gpu_dgeam(blas_handle, 'N', 'N', nO, nV, 1.d0, v_ovoo%f(1,1,j,i), &
         nO, 0.d0, X_ovoo%f(1,1,i,j), nO, X_ovoo%f(1,1,i,j), nO)
    enddo
  enddo

  call gpu_set_stream(blas_handle, gpu_default_stream)

  do i=1,nO
    call gpu_stream_destroy(stream(i))
  enddo

  call gpu_allocate(Y_ovov,nO,nV,nO,nV)

  call gpu_dgemm(blas_handle, 'N','N',nO*nV*nO,nV,nO, &
            -1d0, X_ovoo%f(1,1,1,1), size(X_ovoo%f,1) * size(X_ovoo%f,2) * size(X_ovoo%f,3), &
                  t1%f(1,1)    , size(t1%f,1), &
             0d0, Y_ovov%f(1,1,1,1), size(Y_ovov%f,1) * size(Y_ovov%f,2) * size(Y_ovov%f,3))


  call gpu_copy(v_ovvo, J1)

  call gpu_synchronize()

  do a=1,nV
    call gpu_stream_create(stream(a))
  enddo

  do i = 1, nO
    do beta = 1, nV
      call gpu_set_stream(blas_handle, stream(beta))
      call gpu_dgeam(blas_handle, 'N', 'N', nO, nV, 1.d0, J1%f(1,1,beta,i), &
         nO, 1.d0, Y_ovov%f(1,1,i,beta), nO, J1%f(1,1,beta,i), nO)
    enddo
  enddo

  call gpu_allocate(tmp_cc,cholesky_mo_num,nV,nO)
  call gpu_allocate(J1_tmp,nV,nO,nV,nO)

  call gpu_set_stream(blas_handle, gpu_default_stream)

  type(gpu_double4) :: J1_tmp
  type(gpu_double3) :: tmp_cc

  call gpu_dgemm(blas_handle, 'N','T', cholesky_mo_num*nV, nO, nV, 1.d0, &
      d_cc_space_v_vv_chol%f(1,1,1), cholesky_mo_num*nV, &
      t1%f(1,1), nO, &
      0.d0, tmp_cc%f(1,1,1), cholesky_mo_num*nV)

  call gpu_dgemm(blas_handle, 'T','N', nV*nO, nV*nO, cholesky_mo_num, 1.d0, &
      tmp_cc%f(1,1,1), cholesky_mo_num, d_cc_space_v_vo_chol%f(1,1,1), cholesky_mo_num, &
      0.d0, J1_tmp%f(1,1,1,1), nV*nO)


  call gpu_deallocate(X_ovoo)

  call gpu_synchronize()
  call gpu_deallocate(tmp_cc)

  do i = 1, nO
    do a = 1, nV
      call gpu_set_stream(blas_handle, stream(a))
      call gpu_dgeam(blas_handle, 'N', 'T', nO, nV, 1.d0, J1%f(1,a,1,i), &
         nO*nV, 1.d0, J1_tmp%f(1,1,a,i), nV, J1%f(1,a,1,i), nO*nV)
    enddo
  enddo

  type(gpu_double4) :: X_voov, Z_ovvo

  call gpu_allocate(X_voov,nV,nO,nO,nV)
  call gpu_allocate(Z_ovvo,nO,nV,nV,nO)

  do j = 1, nO
    do beta = 1, nV
      call gpu_set_stream(blas_handle, stream(beta))
      call gpu_dgeam(blas_handle, 'N', 'N', nO, nV, 0.5d0, t2%f(1,j,1,beta), &
         nO*nO, t1%f(j,beta), t1%f(1,1), nO, Y_ovov%f(1,beta,j,1), nO*nV*nO)
    enddo
  enddo

  do b = 1, nV
    call gpu_set_stream(blas_handle, stream(b))
    call gpu_dgeam(blas_handle, 'N', 'N', nV, nO*nO, 1.d0, v_vvoo%f(1,b,1,1), &
         nV*nV, 0.d0, X_voov%f(1,1,1,b), nV, X_voov%f(1,1,1,b), nV)
  enddo

  call gpu_set_stream(blas_handle, gpu_default_stream)

  call gpu_synchronize()
  call gpu_deallocate(J1_tmp)

  call gpu_dgemm(blas_handle, 'N','T',nO*nV,nV*nO,nO*nV, &
             -1d0, Y_ovov%f(1,1,1,1), size(Y_ovov%f,1) * size(Y_ovov%f,2), &
                   X_voov%f(1,1,1,1), size(X_voov%f,1) * size(X_voov%f,2), &
              0d0, Z_ovvo%f(1,1,1,1), size(Z_ovvo%f,1) * size(Z_ovvo%f,2))

  call gpu_synchronize()

  do i = 1, nO
    do a = 1, nV
      call gpu_set_stream(blas_handle, stream(a))
      call gpu_dgeam(blas_handle, 'N', 'N', nO, nV, 1.d0, J1%f(1,a,1,i), &
          nO*nV, 1.d0, Z_ovvo%f(1,1,a,i), nO, J1%f(1,a,1,i), nO*nV)
    enddo
  enddo

  type(gpu_double4) :: X_ovvo, Y_vovo
  call gpu_allocate(Y_vovo,nV,nO,nV,nO)

  do j = 1, nO
    do i = 1, nO
      call gpu_set_stream(blas_handle, stream(i))
      call gpu_dgeam(blas_handle, 'N', 'T', nV, nV, 1.d0, v_vvoo%f(1,1,i,j), &
          nV, -0.5d0, v_vvoo%f(1,1,i,j), nV, Y_vovo%f(1,i,1,j), nO*nV)
    enddo
  enddo

  call gpu_allocate(X_ovvo,nO,nV,nV,nO)

  do j = 1, nO
    do b = 1, nV
      call gpu_set_stream(blas_handle, stream(b))
      call gpu_dgeam(blas_handle, 'N', 'N', nO, nV, 1.d0, t2%f(1,j,1,b), &
          nO*nO, 0.d0, t2%f(1,j,1,b), nO*nO, X_ovvo%f(1,1,b,j), nO)
    enddo
  enddo

  call gpu_set_stream(blas_handle, gpu_default_stream)
  call gpu_synchronize()
  call gpu_deallocate(X_voov)

  call gpu_dgemm(blas_handle, 'N','T',nO*nV,nV*nO,nV*nO, &
             1d0, X_ovvo%f(1,1,1,1), size(X_ovvo%f,1) * size(X_ovvo%f,2), &
                  Y_vovo%f(1,1,1,1), size(Y_vovo%f,1) * size(Y_vovo%f,2), &
             0d0, Z_ovvo%f(1,1,1,1), size(Z_ovvo%f,1) * size(Z_ovvo%f,2))

  call gpu_synchronize()

  do i = 1, nO
    do beta = 1, nV
      call gpu_set_stream(blas_handle, stream(beta))
      call gpu_dgeam(blas_handle, 'N', 'N', nO, nV, 1.d0, J1%f(1,1,beta,i), &
          nO, 1.d0, Z_ovvo%f(1,beta,1,i), nO*nV, J1%f(1,1,beta,i), nO)
    enddo
  enddo

  call gpu_set_stream(blas_handle, gpu_default_stream)
  call gpu_deallocate(Y_ovov)
  call gpu_deallocate(X_ovvo)

  do a = 1, nV
    call gpu_stream_destroy(stream(a))
  enddo

  call gpu_deallocate(Z_ovvo)

end

! K1

subroutine compute_K1_chol(nO,nV,t1,t2,v_ovoo,v_vvoo,v_ovov, &
   d_cc_space_v_ov_chol,d_cc_space_v_vv_chol,K1)
  use gpu

  implicit none

  integer, intent(in)            :: nO,nV
  type(gpu_double2), intent(in)  :: t1
  type(gpu_double4), intent(in)  :: t2, v_vvoo, v_ovov, v_ovoo
  type(gpu_double3), intent(in)  :: d_cc_space_v_ov_chol, d_cc_space_v_vv_chol
  type(gpu_double4), intent(out) :: K1

  type(gpu_double4) :: X, Y, Z

  integer :: a,tmp_a,b,k,l,c,d,tmp_c,tmp_d,i,j,u,v, beta, gam


  call gpu_copy(v_ovov, K1)

  type(gpu_stream) :: stream(nV)
  do a = 1, nV
    call gpu_stream_create(stream(a))
  enddo

  call gpu_allocate(X,nV,nO,nV,nO)

  do i = 1, nO
    do a = 1, nV
      call gpu_set_stream(blas_handle, stream(a))
      call gpu_dgeam(blas_handle, 'N', 'N', nV, nO, -1.d0, v_vvoo%f(1,a,i,1), &
          nV*nV*nO, 0.d0, v_vvoo%f(1,a,i,1), nV*nV*nO, X%f(1,1,a,i), nV)
    enddo
  enddo

  call gpu_allocate(Y,nO,nV,nV,nO)

  do j = 1, nO
    do beta = 1, nV
      call gpu_set_stream(blas_handle, stream(beta))
      call gpu_dgeam(blas_handle, 'N', 'N', nO, nV, 0.5d0, t2%f(1,j,1,beta), &
          nO*nO, t1%f(j,beta), t1%f(1,1), nO, Y%f(1,beta,1,j), nO*nV)
    enddo
  enddo

  call gpu_set_stream(blas_handle, gpu_default_stream)

  call gpu_dgemm(blas_handle, 'N','N',nO*nV*nO,nV,nO, &
            -1d0, v_ovoo%f(1,1,1,1), size(v_ovoo%f,1) * size(v_ovoo%f,2) * size(v_ovoo%f,3), &
                  t1%f(1,1)    , size(t1%f,1), &
            1d0, K1%f(1,1,1,1)    , size(K1%f,1) * size(K1%f,2) * size(K1%f,3))

  type(gpu_double4) :: K1tmp
  type(gpu_double3) :: t1v

  call gpu_allocate(t1v,cholesky_mo_num,nO,nO)

  call gpu_dgemm(blas_handle, 'N','T', cholesky_mo_num*nO, nO, nV, 1.d0, &
    d_cc_space_v_ov_chol%f(1,1,1), cholesky_mo_num*nO, t1%f(1,1), nO, 0.d0, &
    t1v%f(1,1,1), cholesky_mo_num*nO)

  call gpu_allocate(K1tmp,nO,nO,nV,nV)

  call gpu_dgemm(blas_handle, 'T','N', nO*nO, nV*nV, cholesky_mo_num, 1.d0, &
    t1v%f(1,1,1), cholesky_mo_num, d_cc_space_v_vv_chol%f(1,1,1), cholesky_mo_num, 0.d0, &
    K1tmp%f(1,1,1,1), nO*nO)

  call gpu_allocate(Z,nO,nV,nV,nO)
  call gpu_synchronize()

  ! Y(u,beta,b,j) * X(b,j,a,i) = Z(u,beta,a,i)
  call gpu_dgemm(blas_handle, 'N','N',nV*nO,nO*nV,nV*nO, &
             1d0, Y%f(1,1,1,1), size(Y%f,1) * size(Y%f,2), &
                  X%f(1,1,1,1), size(X%f,1) * size(X%f,2), &
             0d0, Z%f(1,1,1,1), size(Z%f,1) * size(Z%f,2))

  call gpu_synchronize()
  call gpu_deallocate(t1v)

   do i = 1, nO
    do beta = 1, nV
      call gpu_set_stream(blas_handle, stream(beta))
      call gpu_dgeam(blas_handle, 'N', 'N', nO, nV, 1.d0, K1%f(1,1,i,beta), &
          nO, 1.d0, K1tmp%f(1,i,1,beta), nO*nO, K1%f(1,1,i,beta), nO)
      call gpu_dgeam(blas_handle, 'N', 'N', nO, nV, 1.d0, K1%f(1,1,i,beta), &
          nO, 1.d0, Z%f(1,beta,1,i), nO*nV, K1%f(1,1,i,beta), nO)
    enddo
  enddo

  call gpu_deallocate(X)
  call gpu_deallocate(Y)

  do a = 1, nV
    call gpu_stream_destroy(stream(a))
  enddo

  call gpu_deallocate(K1tmp)
  call gpu_deallocate(Z)

end
