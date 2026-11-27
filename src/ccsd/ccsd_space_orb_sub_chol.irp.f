use sizet

! H_oo

subroutine compute_H_oo_chol(nO,nV,tau_x,d_cc_space_f_oo, &
    d_cc_space_v_ov_chol,d_cc_space_v_vo_chol,H_oo)
 use sizet
  BEGIN_DOC
  ! Computes the block H_oo used in the T1 residual, the
  ! GPU/Cholesky counterpart of compute_H_oo, called from
  ! run_ccsd_space_orb when do_mo_cholesky is set. H_oo is
  ! initialized to d_cc_space_f_oo and updated with tau_x
  ! contractions over d_cc_space_v_ov_chol and d_cc_space_v_vo_chol.
  ! The BLAS sequence reproduces the compute_H_oo formula
  ! f_oo(u,i) + sum w_oovv(i,j,a,b) * tau(u,j,a,b), the
  ! two-electron contraction formed through the Cholesky vectors
  ! (author-confirmed 2026-10-31).
  END_DOC
  use gpu
  implicit none

  integer, intent(in)           :: nO,nV
  type(gpu_double2), intent(in)    :: d_cc_space_f_oo
  type(gpu_double3), intent(in)    :: d_cc_space_v_ov_chol, d_cc_space_v_vo_chol
  type(gpu_double4), intent(in)    :: tau_x
  type(gpu_double2), intent(out)   :: H_oo

  integer(size_t) :: u, b
  integer         :: a,i,j,k

  type(gpu_double3) :: tau_kau, tmp_vov, tmp_ovv

  call gpu_allocate(tau_kau, int(cholesky_mo_num,4), int(nV,4), int(nO,4))

  type(gpu_blas) :: blas


  !$OMP PARALLEL  &
  !$OMP DEFAULT(SHARED) &
  !$OMP PRIVATE(blas,u,b,tmp_vov,tmp_ovv)

  !$OMP SINGLE
  !$OMP TASK
  call gpu_copy(d_cc_space_f_oo, H_oo)
  !$OMP END TASK
  !$OMP END SINGLE

  call gpu_allocate(tmp_ovv, int(nO,4), int(nV,4), int(nV,4))
  call gpu_allocate(tmp_vov, int(nV,4), int(nO,4), int(nV,4))

  call gpu_blas_create(blas)

  !$OMP DO
  do u=1,nO
    call gpu_dgeam(blas, 'N', 'N', 1, int(nO*nV*nV,4), 1.d0, tau_x%f(u,1,1,1) &
      , int(nO,4), 0.d0, tau_x%f(1,1,1,1), int(nO,4), tmp_ovv%f(1,1,1), 1)
    do b=1,nV
      call gpu_dgeam(blas, 'T', 'T', int(nV,4), int(nO,4), 1.d0, &
           tmp_ovv%f(1,1,b), int(nO,4), 0.d0, &
           tmp_ovv%f(1,1,b), int(nO,4), tmp_vov%f(1,1,b), int(nV,4))
    enddo
    call gpu_dgemm(blas, 'N','T',int(cholesky_mo_num,4),int(nV,4), &
      int(nO*nV,4),1.d0, &
      d_cc_space_v_ov_chol%f(1,1,1), int(cholesky_mo_num,4), tmp_vov%f(1,1,1), &
      nV, &
      0.d0, tau_kau%f(1,1,u), int(cholesky_mo_num,4))
  enddo
  !$OMP END DO

  call gpu_blas_destroy(blas)

  call gpu_deallocate(tmp_vov)
  call gpu_deallocate(tmp_ovv)

  !$OMP TASKWAIT
  !$OMP END PARALLEL

  call gpu_dgemm(blas_handle, 'T', 'N', int(nO,4), int(nO,4) &
    , int(cholesky_mo_num*nV,4), 1.d0, tau_kau%f(1,1,1) &
    , int(cholesky_mo_num*nV,4), d_cc_space_v_vo_chol%f(1,1,1) &
    , int(cholesky_mo_num*nV,4), 1.d0, H_oo%f(1,1), int(nO,4))

  call gpu_synchronize()
  call gpu_deallocate(tau_kau)
end

! H_vv

subroutine compute_H_vv_chol(nO,nV,tau_x,d_cc_space_f_vv, &
         d_cc_space_v_ov_chol,H_vv)
 use sizet
  BEGIN_DOC
  ! Computes the block H_vv used in the T1 residual, the
  ! GPU/Cholesky counterpart of compute_H_vv, called from
  ! run_ccsd_space_orb when do_mo_cholesky is set. H_vv is
  ! initialized to d_cc_space_f_vv and updated with tau_x
  ! contractions over d_cc_space_v_ov_chol.
  ! The BLAS sequence reproduces the compute_H_vv formula
  ! f_vv(a,beta) - sum w_vvoo(a,b,i,j) * tau(i,j,beta,b), the
  ! contraction formed through d_cc_space_v_ov_chol
  ! (author-confirmed 2026-10-31).
  END_DOC
  use gpu
  implicit none

  integer, intent(in)              :: nO,nV
  type(gpu_double2), intent(in)    :: d_cc_space_f_vv
  type(gpu_double3), intent(in)    :: d_cc_space_v_ov_chol
  type(gpu_double4), intent(in)    :: tau_x
  type(gpu_double2), intent(out)   :: H_vv

  integer(size_t) :: a, b
  integer         :: i,j,u,k,beta

  type(gpu_double3) :: tau_kia, tmp_oov

  call gpu_allocate(tau_kia, int(cholesky_mo_num,4), int(nO,4), int(nV,4))

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
  call gpu_allocate(tmp_oov, int(nO,4), int(nO,4), int(nV,4))

  !$OMP DO
  do a = 1, nV
    do b=1,nV
      call gpu_dgeam(blas, 'N', 'N', int(nO,4), int(nO,4), 1.d0, &
        tau_x%f(1,1,a,b), int(nO,4), 0.d0, &
        tau_x%f(1,1,a,b), int(nO,4), tmp_oov%f(1,1,b), int(nO,4))
    enddo
    call gpu_dgemm(blas, 'N', 'T', int(cholesky_mo_num,4), int(nO,4), &
      int(nO*nV,4) &
      , 1.d0, d_cc_space_v_ov_chol%f(1,1,1), int(cholesky_mo_num,4), &
      tmp_oov%f(1,1,1) &
      , nO, 0.d0, tau_kia%f(1,1,a), int(cholesky_mo_num,4))
  enddo
  !$OMP END DO

  call gpu_blas_destroy(blas)

  call gpu_deallocate(tmp_oov)
  !$OMP TASKWAIT
  !$OMP END PARALLEL

  call gpu_dgemm(blas_handle, 'T', 'N', int(nV,4), int(nV,4) &
    , int(cholesky_mo_num*nO,4), -1.d0, tau_kia%f(1,1,1) &
    , int(cholesky_mo_num*nO,4), d_cc_space_v_ov_chol%f(1,1,1) &
    , int(cholesky_mo_num*nO,4), 1.d0, H_vv%f(1,1), int(nV,4))

  call gpu_synchronize()
  call gpu_deallocate(tau_kia)
end

! H_vo
subroutine compute_H_vo_chol(nO,nV,t1,d_cc_space_f_vo, &
         d_cc_space_v_ov_chol,d_cc_space_v_vo_chol, H_vo)
 use sizet
  BEGIN_DOC
  ! Computes the block H_vo used in the T1 residual, the
  ! GPU/Cholesky counterpart of compute_H_vo, called from
  ! run_ccsd_space_orb when do_mo_cholesky is set. H_vo is
  ! initialized to d_cc_space_f_vo and updated with t1
  ! contractions over d_cc_space_v_ov_chol and d_cc_space_v_vo_chol.
  ! The BLAS sequence reproduces the compute_H_vo formula
  ! f_vo(a,i) + sum w_vvoo(a,b,i,j) * t1(j,b); the pre-GPU
  ! transpose loop (commit fa6d1419~1) is a gpu_dgeam that
  ! copies each column in transposed order (author-confirmed
  ! 2026-10-31).
  END_DOC
  use gpu
  implicit none

  integer, intent(in)            :: nO,nV
  type(gpu_double2), intent(in)  :: t1, d_cc_space_f_vo
  type(gpu_double3), intent(in)  :: d_cc_space_v_ov_chol, d_cc_space_v_vo_chol
  type(gpu_double2), intent(out) :: H_vo

  integer(size_t) :: i
  integer         :: a,b,j,u,k

  type(gpu_double1) :: tmp_k
  type(gpu_double3) :: tmp, tmp2

  call gpu_copy(d_cc_space_f_vo, H_vo)

  call gpu_allocate(tmp_k, int(cholesky_mo_num,4))

  call gpu_dgemm(blas_handle, 'N', 'N', int(cholesky_mo_num,4), 1, &
    int(nO*nV,4), 2.d0 &
    , d_cc_space_v_ov_chol%f(1,1,1), int(cholesky_mo_num,4), t1%f(1,1), nO*nV, &
    0.d0 &
    , tmp_k%f(1), int(cholesky_mo_num,4))

  call gpu_dgemm(blas_handle, 'T', 'N', int(nV*nO,4), 1, &
    int(cholesky_mo_num,4), 1.d0 &
    , d_cc_space_v_vo_chol%f(1,1,1), int(cholesky_mo_num,4), tmp_k%f(1) &
    , int(cholesky_mo_num,4), 1.d0, H_vo%f(1,1), int(nV*nO,4))

  call gpu_deallocate(tmp_k)


  call gpu_allocate(tmp,  int(cholesky_mo_num,4), int(nO,4), int(nO,4))

  call gpu_dgemm(blas_handle, 'N', 'T', int(cholesky_mo_num*nO,4), int(nO,4) &
    , int(nV,4), 1.d0, d_cc_space_v_ov_chol%f(1,1,1) &
    , int(cholesky_mo_num*nO,4), t1%f(1,1), nO, 0.d0, tmp%f(1,1,1) &
    , int(cholesky_mo_num*nO,4))

  call gpu_allocate(tmp2, int(cholesky_mo_num,4), int(nO,4), int(nO,4))

  type(gpu_stream) :: stream(nO)
  do i=1,nO
    call gpu_stream_create(stream(i))
  enddo

  !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i,j)
  do i=1,nO
    do j=1,nO
      call gpu_set_stream(blas_handle,stream(j))
      call gpu_dgeam(blas_handle, 'N', 'N', int(cholesky_mo_num,4), 1, 1.d0, &
        tmp%f(1,i,j), int(cholesky_mo_num,4), 0.d0, &
        tmp%f(1,i,j), int(cholesky_mo_num,4), tmp2%f(1,j,i), &
        int(cholesky_mo_num,4))
    enddo
  enddo
  !$OMP END PARALLEL DO

  call gpu_set_stream(blas_handle,gpu_default_stream)
  call gpu_synchronize()

  do i=1,nO
    call gpu_stream_destroy(stream(i))
  enddo
  call gpu_deallocate(tmp)

  call gpu_dgemm(blas_handle, 'T', 'N', int(nV,4), int(nO,4) &
    , int(cholesky_mo_num*nO,4), -1.d0, d_cc_space_v_ov_chol%f(1,1,1) &
    , int(cholesky_mo_num*nO,4), tmp2%f(1,1,1), int(cholesky_mo_num*nO,4), &
    1.d0 &
    , H_vo%f(1,1), int(nV,4))

  call gpu_synchronize()
  call gpu_deallocate(tmp2)
end

! R1

subroutine compute_r1_space_chol(nO,nV,t1,t2,tau,H_oo,H_vv,H_vo,r1,max_r1,d_cc_space_f_ov,d_cc_space_f_vo, &
    d_cc_space_v_voov, d_cc_space_v_ovov, d_cc_space_v_oovo, d_cc_space_v_vo_chol, d_cc_space_v_vv_chol)
 use sizet
  BEGIN_DOC
  ! Computes the residual r1 of the T1 equation, the
  ! GPU/Cholesky counterpart of compute_r1_space, called from
  ! run_ccsd_space_orb when do_mo_cholesky is set. r1 is
  ! initialized to d_cc_space_f_ov and updated from t1, t2, tau,
  ! H_oo, H_vv, H_vo, d_cc_space_f_vo, the integrals v_voov,
  ! v_ovov, v_oovo, and the Cholesky vectors v_vo_chol and
  ! v_vv_chol. r1 is sign-flipped at the end for consistency
  ! with the spin-orbital code; max_r1 is the largest |r1|.
  ! The BLAS sequence reproduces the non-Cholesky r1 formula of
  ! compute_r1_space; the pre-GPU X_voov expression (commit
  ! fa6d1419~1) is two gpu_dgeam calls on t2 and t1 slices
  ! (author-confirmed 2026-10-31).
  END_DOC
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
  integer(size_t)               :: u,i,beta,a,b
  integer                       :: j

  type(gpu_stream) :: stream(nV)

  do a=1,nV
    call gpu_stream_create(stream(a))
  enddo

  type(gpu_double2) :: X_oo
  call gpu_allocate(X_oo,int(nO,4),int(nO,4))

  call gpu_copy(d_cc_space_f_ov, r1)

  call gpu_set_stream(blas_handle, stream(1))
  call gpu_dgemm(blas_handle, 'N','N', int(nO,4), int(nV,4), int(nV,4), &
             1d0, t1%f(1,1)  , size(t1%f,1), &
                  H_vv%f(1,1), size(H_vv%f,1), &
             1d0, r1%f(1,1)  , size(r1%f,1))

  call gpu_dgemm(blas_handle, 'N','N', int(nO,4), int(nV,4), int(nO,4), &
             -1d0, H_oo%f(1,1), size(H_oo%f,1), &
                   t1%f(1,1)  , size(t1%f,1), &
              1d0, r1%f(1,1), size(r1%f,1))

  call gpu_set_stream(blas_handle, stream(nV))
  call gpu_dgemm(blas_handle, 'N','N', int(nO,4), int(nO,4), int(nV,4), &
             -2d0, t1%f(1,1), size(t1%f,1), &
                   d_cc_space_f_vo%f(1,1), size(d_cc_space_f_vo%f,1), &
              0d0, X_oo%f(1,1), size(X_oo%f,1))

  call gpu_synchronize()
  call gpu_set_stream(blas_handle, gpu_default_stream)

  call gpu_dgemm(blas_handle, 'T','N', int(nO,4), int(nV,4), int(nO,4), &
             1d0, X_oo%f(1,1), size(X_oo%f,2), &
                  t1%f(1,1)  , size(t1%f,1), &
             1d0, r1%f(1,1)  , size(r1%f,1))



  type(gpu_double4) :: X_voov
  call gpu_allocate(X_voov, int(nV,4), int(nO,4), int(nO,4), int(nV,4))

  do i=1,nO
    do beta=1,nV
      call gpu_set_stream(blas_handle, stream(beta))
      call gpu_dgeam(blas_handle, 'T', 'T', int(nV,4), int(nO,4), -1.d0 &
        , t2%f(1,i,1,beta), int(nO*nO,4), t1%f(i,beta), t1%f(1,1), int(nO,4) &
        , X_voov%f(1,i,1,beta), int(nV*nO,4))
    enddo
  enddo

  do beta=1,nV
    call gpu_set_stream(blas_handle, stream(beta))
    call gpu_dgeam(blas_handle, 'N', 'T', int(nV,4), int(nO*nO,4), 1.d0 &
      , X_voov%f(1,1,1,beta), int(nV,4), 2.d0, t2%f(1,1,1,beta), int(nO*nO,4) &
      , X_voov%f(1,1,1,beta), int(nV,4))
  enddo

  call gpu_synchronize()
  call gpu_deallocate(X_oo)

  call gpu_set_stream(blas_handle, gpu_default_stream)

  call gpu_dgemv(blas_handle, 'T', int(nV*nO,4), int(nO*nV,4), &
             1d0, X_voov%f(1,1,1,1), size(X_voov%f,1) * size(X_voov%f,2), &
                  H_vo%f(1,1)  , 1, &
             1d0, r1%f(1,1)    , 1)

  type(gpu_double4) :: X_ovov
  call gpu_allocate(X_ovov, int(nO,4), int(nV,4), int(nO,4), int(nV,4))

  do beta = 1, nV
    call gpu_set_stream(blas_handle, stream(beta))
    do u=1,nO
      call gpu_dgeam(blas_handle, 'N', 'T', int(nO,4), int(nV,4), -1.d0 &
        , d_cc_space_v_ovov%f(1,1,u,beta), int(nO,4), 2.d0 &
        , d_cc_space_v_voov%f(1,u,1,beta), int(nV*nO,4), X_ovov%f(1,1,u,beta) &
        , int(nO,4))
    enddo
  enddo

  call gpu_set_stream(blas_handle, gpu_default_stream)
  call gpu_synchronize()
  call gpu_deallocate(X_voov)

  call gpu_dgemv(blas_handle, 'T', int(nO*nV,4), int(nO*nV,4), &
             1d0, X_ovov%f(1,1,1,1), size(X_ovov%f,1) * size(X_ovov%f,2), &
                  t1%f(1,1), 1, &
             1d0, r1%f(1,1), 1)


  integer :: iblock, block_size, nVmax
  type(gpu_double4) :: W_vvov, W_vvov_tmp, T_vvoo

  block_size = 16
  call gpu_allocate(T_vvoo, int(nV,4),int(nV,4),int(nO,4),int(nO,4))

  call gpu_dgeam(blas_handle, 'T', 'N', int(nV*nV,4), int(nO*nO,4), 1.d0 &
    , tau%f(1,1,1,1), int(nO*nO,4), 0.d0, T_vvoo%f(1,1,1,1), int(nV*nV,4) &
    , T_vvoo%f(1,1,1,1), int(nV*nV,4))

  call gpu_allocate(W_vvov,int(nV,4), int(nV,4),int(nO,4),block_size)
  call gpu_allocate(W_vvov_tmp, int(nV,4),int(nO,4),int(nV,4),block_size)

  do iblock = 1, nV, block_size
    nVmax = min(block_size,int(nV-iblock+1,4))

    call gpu_dgemm(blas_handle, 'T', 'N', int(nV*nO,4), int(nV*nVmax,4) &
      , int(cholesky_mo_num,4), 1.d0, d_cc_space_v_vo_chol%f(1,1,1), &
      int(cholesky_mo_num,4) &
      , d_cc_space_v_vv_chol%f(1,1,iblock), int(cholesky_mo_num,4), 0.d0 &
      , W_vvov_tmp%f(1,1,1,1), int(nV*nO,4))

    call gpu_synchronize()
    do b=1,nV
      call gpu_set_stream(blas_handle, stream(b))
      do i=1,nO
        call gpu_dgeam(blas_handle, 'N', 'N', int(nV,4), nVmax, 2.d0 &
          , W_vvov_tmp%f(1,i,b,1), int(nV*nO*nV,4), 0.d0 &
          , W_vvov_tmp%f(1,i,b,1), int(nV*nO*nV,4), W_vvov%f(1,b,i,1) &
          , int(nV*nV*nO,4))
      enddo
    enddo

    call gpu_synchronize()

    do beta = 1,  nVmax
      call gpu_set_stream(blas_handle, stream(beta))
      call gpu_dgeam(blas_handle, 'N', 'T', int(nV,4), int(nV*nO,4), 1.d0 &
        , W_vvov%f(1,1,1,beta), int(nV,4), -1.d0, W_vvov_tmp%f(1,1,1,beta) &
        , int(nV*nO,4), W_vvov%f(1,1,1,beta), int(nV,4))
    enddo
    call gpu_synchronize()

    call gpu_dgemm(blas_handle, 'T','N',int(nO,4),nVmax,int(nO*nV*nV,4), &
             1d0, T_vvoo%f(1,1,1,1), int(nV*nV*nO,4), &
                  W_vvov%f(1,1,1,1), nO*nV*nV, &
             1d0, r1%f(1,iblock), int(nO,4))
  enddo

  call gpu_deallocate(X_ovov)

  type(gpu_double4) :: W_oovo
  call gpu_allocate(W_oovo, int(nO,4),int(nO,4),int(nV,4),int(nO,4))

  do u = 1, nO
    do a = 1, nV
      call gpu_set_stream(blas_handle, stream(a))
      call gpu_dgeam(blas_handle, 'N', 'T', int(nO,4), int(nO,4), 2.d0 &
        , d_cc_space_v_oovo%f(1,1,a,u), int(nO,4), -1.d0 &
        , d_cc_space_v_oovo%f(1,1,a,u), int(nO,4), W_oovo%f(1,1,a,u), int(nO,4))
    enddo
  enddo

  call gpu_set_stream(blas_handle, gpu_default_stream)
  call gpu_synchronize()

  call gpu_deallocate(W_vvov)
  call gpu_deallocate(T_vvoo)

  ! Change the sign for consistency with the code in spin orbitals
  call gpu_dgemm(blas_handle, 'T', 'N', int(nO,4), int(nV,4), int(nO*nO*nV,4) &
    , 1d0, W_oovo%f(1,1,1,1) &
    , size(W_oovo%f,1) * size(W_oovo%f,2) * size(W_oovo%f,3), tau%f(1,1,1,1) &
    , size(tau%f,1) * size(tau%f,2) * size(tau%f,3), -1d0, r1%f(1,1) &
    , size(r1%f,1))

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
 use sizet
  BEGIN_DOC
  ! Computes the residual r2 of the T2 equation, the
  ! GPU/Cholesky counterpart of compute_r2_space, called from
  ! run_ccsd_space_orb when do_mo_cholesky is set. r2 is
  ! initialized to the integral d_cc_space_v_oovv and updated
  ! from t1, t2, tau, the global integrals, and the Cholesky
  ! vectors, using compute_A1_chol, compute_g_occ_chol,
  ! compute_g_vir_chol, compute_J1_chol and compute_K1_chol.
  ! r2 is sign-flipped at the end for consistency with the
  ! spin-orbital code; max_r2 is the largest absolute value of r2
  ! over its indices, i.e. max_r2 = max |r2|, consistent with
  ! max_r1 (dabs) and with the spin-orbital reference. The
  ! dabs is taken in the per-thread leaf max_r2_local (see
  ! compute_r2_space), so it cannot be moved to the reduction:
  ! max(|a|,|b|) != |max(a,b)|. The missing-dabs defect was fixed
  ! in issue #396 (2026-10-31): the dabs was added to the leaf
  ! accumulation.
  END_DOC
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
  integer(size_t)               :: i,j,beta,gam,a,b
  integer                       :: u,v
  double precision              :: max_r2_local

  type(gpu_stream) :: stream(nV)


  call gpu_copy(d_cc_space_v_oovv, r2)

  type(gpu_double4) :: A1
  call gpu_allocate(A1,int(nO,4),int(nO,4),int(nO,4),int(nO,4))
  call compute_A1_chol(nO,nV,t1,t2,tau,d_cc_space_v_vooo, &
      d_cc_space_v_oooo, d_cc_space_v_vvoo, A1)

  call gpu_dgemm(blas_handle, 'N','N',int(nO*nO,4),int(nV*nV,4),int(nO*nO,4), &
             1d0, A1%f(1,1,1,1), size(A1%f,1) * size(A1%f,2), &
                  tau%f(1,1,1,1), size(tau%f,1) * size(tau%f,2), &
             1d0, r2%f(1,1,1,1), size(r2%f,1) * size(r2%f,2))

  call gpu_deallocate(A1)

  integer(size_t) :: iblock
  integer         :: block_size, k
  block_size = 16
  type(gpu_double3) :: tmp_cc, B1, tmpB1
  type(gpu_double2) :: tmp_cc2

  call gpu_allocate(tmp_cc,int(cholesky_mo_num,4),int(nV,4),int(nV,4))
  call gpu_dgemm(blas_handle, 'N', 'N', int(cholesky_mo_num*nV,4), int(nV,4) &
    , int(nO,4), 1.d0, d_cc_space_v_vo_chol%f(1,1,1) &
    , int(cholesky_mo_num*nV,4), t1%f(1,1), nO, 0.d0, tmp_cc%f(1,1,1) &
    , int(cholesky_mo_num*nV,4))

  call gpu_synchronize()

  type(gpu_blas) :: blas

  !$OMP PARALLEL PRIVATE(gam, iblock, B1, tmpB1, tmp_cc2, beta, b, a, blas)
  call gpu_allocate(B1,int(nV,4),int(nV,4),block_size)
  call gpu_allocate(tmpB1,int(nV,4),block_size,int(nV,4))
  call gpu_allocate(tmp_cc2,int(cholesky_mo_num,4),int(nV,4))

  call gpu_blas_create(blas)

  !$OMP DO
  do gam = 1, nV

    call gpu_dgeam(blas, 'N', 'N', int(cholesky_mo_num,4), int(nV,4), 1.d0 &
      , d_cc_space_v_vv_chol%f(1,1,gam), int(cholesky_mo_num,4), -1.d0 &
      , tmp_cc%f(1,1,gam), int(cholesky_mo_num,4), tmp_cc2%f(1,1), &
      int(cholesky_mo_num,4))

    do iblock = 1, nV, block_size

        call gpu_dgemm(blas, 'T', 'N' &
          , int(nV*min(block_size,int(nV-iblock+1,4)),4), int(nV,4) &
          , int(cholesky_mo_num,4), -1.d0, tmp_cc%f(1,1,iblock), &
          int(cholesky_mo_num,4) &
          , d_cc_space_v_vv_chol%f(1,1,gam), int(cholesky_mo_num,4), 0.d0 &
          , tmpB1%f(1,1,1), int(nV*block_size,4))

        call gpu_dgemm(blas, 'T', 'N' &
          , int(nV*min(block_size,int(nV-iblock+1,4)),4), int(nV,4) &
          , int(cholesky_mo_num,4), 1.d0, d_cc_space_v_vv_chol%f(1,1,iblock) &
          , int(cholesky_mo_num,4), tmp_cc2%f(1,1), int(cholesky_mo_num,4), &
          1.d0 &
          , tmpB1%f(1,1,1), int(nV*block_size,4))

        do beta = iblock, min(nV, iblock+block_size-1)
          call gpu_dgeam(blas, 'N', 'N', int(nV,4), int(nV,4), 1.d0 &
            , tmpB1%f(1,beta-iblock+1,1), int(nV*block_size,4), 0.d0 &
            , B1%f(1,1,beta-iblock+1), int(nV,4), B1%f(1,1,beta-iblock+1) &
            , int(nV,4))
        enddo

        call gpu_dgemm(blas, 'N', 'N', int(nO*nO,4) &
          , int(min(block_size,int(nV-iblock+1,4)),4), int(nV*nV,4), 1d0 &
          , tau%f(1,1,1,1), size(tau%f,1) * size(tau%f,2), B1%f(1,1,1) &
          , size(B1%f ,1) * size(B1%f ,2), 1d0, r2%f(1,1,iblock,gam) &
          , size(r2%f ,1) * size(r2%f ,2))
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
  call gpu_allocate(X_oovv,int(nO,4),int(nO,4),int(nV,4),int(nV,4))
  call gpu_copy(t2,X_oovv)

  type(gpu_double2) :: g_occ, g_vir
  call gpu_allocate(g_vir,int(nV,4),int(nV,4))
  call gpu_allocate(g_occ,int(nO,4),int(nO,4))
  call compute_g_vir_chol(nO,nV,t1,t2,H_vv,d_cc_space_f_vo, &
    d_cc_space_v_ov_chol, d_cc_space_v_vv_chol, g_vir)
  call compute_g_occ_chol(nO,nV,t1,t2,H_oo, &
    d_cc_space_f_vo, d_cc_space_v_ov_chol, d_cc_space_v_oo_chol, d_cc_space_v_ovoo, g_occ)

  type(gpu_double4) :: Y_oovv
  call gpu_allocate(Y_oovv,int(nO,4),int(nO,4),int(nV,4),int(nV,4))

  call gpu_dgemm(blas_handle, 'N', 'N', int(nO*nO*nV,4), int(nV,4), int(nV,4) &
    , 1d0, X_oovv%f(1,1,1,1) &
    , size(X_oovv%f,1) * size(X_oovv%f,2) * size(X_oovv%f,3), g_vir%f(1,1) &
    , size(g_vir%f,1), 0d0, Y_oovv%f(1,1,1,1) &
    , size(Y_oovv%f,1) * size(Y_oovv%f,2) * size(Y_oovv%f,3))

  call gpu_dgemm(blas_handle, 'N','N',int(nO,4),int(nO*nV*nV,4),int(nO,4), &
             -1d0, g_occ%f(1,1), size(g_occ%f,1), &
                  t2%f(1,1,1,1)    , size(t2%f,1), &
             1d0, Y_oovv%f(1,1,1,1), size(Y_oovv%f,1))

  call gpu_dgemm(blas_handle, 'N', 'N', int(nO*nO*nV,4), int(nV,4), int(nO,4) &
    , -1d0, d_cc_space_v_oovo%f(1,1,1,1) &
    , size(cc_space_v_oovo,1) * size(cc_space_v_oovo,2) * &
    size(cc_space_v_oovo,3), t1%f(1,1), size(t1%f,1), 1d0, Y_oovv%f(1,1,1,1) &
    , size(Y_oovv%f,1) * size(Y_oovv%f,2) * size(Y_oovv%f,3))


  call gpu_dgeam(blas_handle, 'N', 'N', int(nO*nO,4), int(nV*nV,4), 1.d0 &
    , Y_oovv%f(1,1,1,1), int(nO*nO,4), 1.d0, r2%f(1,1,1,1), int(nO*nO,4) &
    , r2%f(1,1,1,1), int(nO*nO,4))

  call gpu_synchronize()
  call gpu_deallocate(X_oovv)

  call gpu_deallocate(g_vir)
  call gpu_deallocate(g_occ)

  type(gpu_double4) :: X_vovo, Y_oovo
  call gpu_allocate(X_vovo,int(nV,4),int(nO,4),int(nV,4),int(nO,4))

  do a=1,nV
    call gpu_stream_create(stream(a))
  enddo

  do gam = 1, nV
    call gpu_set_stream(blas_handle, stream(gam))
    do beta = 1, nV
      call gpu_dgeam(blas_handle, 'N', 'T', int(nO,4), int(nO,4), 1.d0 &
        , r2%f(1,1,beta,gam), int(nO,4), 1.d0, Y_oovv%f(1,1,gam,beta) &
        , int(nO,4), r2%f(1,1,beta,gam), int(nO,4))
    enddo
  enddo

  do i = 1, nO
    do gam = 1, nV
      call gpu_set_stream(blas_handle, stream(gam))
      call gpu_dgeam(blas_handle, 'T', 'N', int(nV,4), int(nO,4), 1.d0 &
        , d_cc_space_v_ovvo%f(1,1,gam,i), int(nO,4), 0.d0, X_vovo%f(1,1,gam,i) &
        , int(nV,4), X_vovo%f(1,1,gam,i), int(nV,4))
    enddo
  enddo

  do a=1,nV
    call gpu_stream_destroy(stream(a))
  enddo
  call gpu_set_stream(blas_handle, gpu_default_stream)



  call gpu_allocate(Y_oovo,int(nO,4),int(nO,4),int(nV,4),int(nO,4))

  !$OMP PARALLEL PRIVATE(blas, iblock, gam, X_vovv)
  call gpu_blas_create(blas)
  type(gpu_double4) :: X_vovv
  call gpu_allocate(X_vovv,int(nV,4),int(nO,4),int(nV,4),block_size)
  !$OMP DO
  do iblock = 1, nV, block_size
    do gam = iblock, min(nV, iblock+block_size-1)
      call gpu_dgemm(blas, 'T', 'N', int(nV,4), int(nO*nV,4), &
        int(cholesky_mo_num,4) &
        , 1.d0, d_cc_space_v_vv_chol%f(1,1,gam), int(cholesky_mo_num,4) &
        , d_cc_space_v_ov_chol%f(1,1,1), int(cholesky_mo_num,4), 0.d0 &
        , X_vovv%f(1,1,1,gam-iblock+1), int(nV,4))

    enddo

    call gpu_dgemm(blas, 'N','N', int(nO,4), &
             int(nO*nV*min(block_size,int(nV-iblock+1,4)),4),int(nV,4), &
             1.d0, t1%f(1,1)    , size(t1%f,1), &
             X_vovv%f(1,1,1,1), size(X_vovv%f,1), &
             0d0, Y_oovv%f(1,1,1,iblock), size(Y_oovv%f,1))
  enddo
  !$OMP END DO

  call gpu_blas_destroy(blas)
  call gpu_deallocate(X_vovv)
  !$OMP END PARALLEL

  call gpu_dgemm(blas_handle, 'N','N',int(nO,4),int(nO*nV*nO,4),int(nV,4), &
             1d0, t1%f(1,1), size(t1%f,1), &
                  X_vovo%f(1,1,1,1), size(X_vovo%f,1), &
             0d0, Y_oovo%f(1,1,1,1), size(Y_oovo%f,1))

  call gpu_dgemm(blas_handle, 'N', 'N', int(nO*nO*nV,4), int(nV,4), int(nO,4) &
    , -1d0, Y_oovo%f(1,1,1,1) &
    , size(Y_oovo%f,1) * size(Y_oovo%f,2) * size(Y_oovo%f,3), t1%f(1,1) &
    , size(t1%f,1), 1d0, Y_oovv%f(1,1,1,1) &
    , size(Y_oovv%f,1) * size(Y_oovv%f,2) * size(Y_oovv%f,3))

  call gpu_synchronize()
  call gpu_deallocate(X_vovo)
  call gpu_deallocate(Y_oovo)

  do a=1,nV
    call gpu_stream_create(stream(a))
  enddo

  do gam = 1, nV
    call gpu_set_stream(blas_handle, stream(gam))
    do beta = 1, nV
      call gpu_dgeam(blas_handle, 'T', 'N', int(nO,4), int(nO,4), 1.d0 &
        , Y_oovv%f(1,1,beta,gam), int(nO,4), 1.d0, r2%f(1,1,beta,gam) &
        , int(nO,4), r2%f(1,1,beta,gam), int(nO,4))
    enddo
    do j=1,nO
      call gpu_dgeam(blas_handle, 'N', 'N', int(nO,4), int(nV,4), 1.d0 &
        , r2%f(1,j,1,gam), int(nO*nO,4), 1.d0, Y_oovv%f(1,j,gam,1) &
        , int(nO*nO*nV,4), r2%f(1,j,1,gam), int(nO*nO,4))
    enddo
  enddo

  call gpu_set_stream(blas_handle, gpu_default_stream)


  call gpu_synchronize()
  call gpu_deallocate(Y_oovv)

  type(gpu_double4) :: X_ovvo
  type(gpu_double3) :: tcc, tcc2
  call gpu_allocate(tcc2,int(cholesky_mo_num,4),int(nV,4),int(nO,4))
  call gpu_allocate(X_ovvo,int(nO,4),int(nV,4),int(nV,4),int(nO,4))
  call gpu_allocate(tcc,int(cholesky_mo_num,4),int(nO,4),int(nV,4))

  call gpu_dgemm(blas_handle, 'N', 'T', int(cholesky_mo_num*nV,4), int(nO,4) &
    , int(nV,4), 1.d0, d_cc_space_v_vv_chol%f(1,1,1) &
    , int(cholesky_mo_num*nV,4), t1%f(1,1), nO, 0.d0, tcc2%f(1,1,1) &
    , int(cholesky_mo_num*nV,4))

  call gpu_dgemm(blas_handle, 'N', 'N', int(cholesky_mo_num*nO,4), int(nV,4) &
    , int(nO,4), 1.d0, d_cc_space_v_oo_chol%f(1,1,1) &
    , int(cholesky_mo_num*nO,4), t1%f(1,1), nO, 0.d0, tcc%f(1,1,1) &
    , int(cholesky_mo_num*nO,4))

  call gpu_dgemm(blas_handle, 'T', 'N', int(nO*nV,4), int(nV*nO,4) &
    , int(cholesky_mo_num,4), 1.d0, tcc%f(1,1,1), int(cholesky_mo_num,4), &
    tcc2%f(1,1,1) &
    , int(cholesky_mo_num,4), 0.d0, X_ovvo%f(1,1,1,1), int(nO*nV,4))

  call gpu_synchronize()


  do gam = 1, nV
    call gpu_set_stream(blas_handle, stream(gam))
    do j=1,nO
      call gpu_dgeam(blas_handle, 'N', 'N', int(nO,4), int(nV,4), -1.d0 &
        , X_ovvo%f(1,1,gam,j), int(nO,4), 1.d0, r2%f(1,j,1,gam), int(nO*nO,4) &
        , r2%f(1,j,1,gam), int(nO*nO,4))
    enddo
    do beta = 1, nV
      call gpu_dgeam(blas_handle, 'T', 'N', int(nO,4), int(nO,4), -1.d0 &
        , X_ovvo%f(1,gam,beta,1), int(nO*nV*nV,4), 1.d0, r2%f(1,1,beta,gam) &
        , int(nO,4), r2%f(1,1,beta,gam), int(nO,4))
    enddo
  enddo

  call gpu_set_stream(blas_handle, gpu_default_stream)

  call gpu_synchronize
  call gpu_deallocate(tcc)
  call gpu_deallocate(tcc2)
  call gpu_deallocate(X_ovvo)


  type(gpu_double4) :: J1, K1
  type(gpu_double4) :: Y_voov, Z_ovov


  call gpu_allocate(J1,int(nO,4),int(nV,4),int(nV,4),int(nO,4))
  call compute_J1_chol(nO,nV,t1,t2,d_cc_space_v_ovvo,d_cc_space_v_ovoo, &
       d_cc_space_v_vvoo,d_cc_space_v_vo_chol,d_cc_space_v_vv_chol,J1)

  call gpu_allocate(K1,int(nO,4),int(nV,4),int(nO,4),int(nV,4))
  call compute_K1_chol(nO,nV,t1,t2,d_cc_space_v_ovoo,d_cc_space_v_vvoo, &
       d_cc_space_v_ovov,d_cc_space_v_ov_chol,d_cc_space_v_vv_chol,K1)


  call gpu_allocate(X_ovvo,int(nO,4),int(nV,4),int(nV,4),int(nO,4))
  call gpu_allocate(Y_voov,int(nV,4),int(nO,4),int(nO,4),int(nV,4))

  do a=1, nV
    call gpu_set_stream(blas_handle, stream(a))
    do i=1, nO
      call gpu_dgeam(blas_handle, 'N', 'N', int(nO,4), int(nV,4), 1.d0 &
        , J1%f(1,a,1,i), int(nO*nV,4), -0.5d0, K1%f(1,a,i,1), int(nO*nV*nO,4) &
        , X_ovvo%f(1,1,a,i), int(nO,4))
      call gpu_dgeam(blas_handle, 'T', 'T', int(nV,4), int(nO,4), 2.d0 &
        , t2%f(1,i,1,a), int(nO*nO,4), -1.d0, t2%f(1,i,a,1), int(nO*nO*nV,4) &
        , Y_voov%f(1,1,i,a), int(nV,4))
    enddo
  enddo

  call gpu_allocate(Z_ovov,int(nO,4),int(nV,4),int(nO,4),int(nV,4))

  call gpu_synchronize()
  call gpu_deallocate(J1)
  call gpu_set_stream(blas_handle, gpu_default_stream)


  call gpu_dgemm(blas_handle, 'N','N', int(nO*nV,4),int(nO*nV,4),int(nV*nO,4), &
             1d0, X_ovvo%f(1,1,1,1), size(X_ovvo%f,1) * size(X_ovvo%f,2), &
                  Y_voov%f(1,1,1,1), size(Y_voov%f,1) * size(Y_voov%f,2), &
             0d0, Z_ovov%f(1,1,1,1), size(Z_ovov%f,1) * size(Z_ovov%f,2))

  call gpu_synchronize()
  call gpu_deallocate(Y_voov)
  call gpu_deallocate(X_ovvo)

  type(gpu_double4) :: Y_ovov, X_ovov
  call gpu_allocate(X_ovov,int(nO,4),int(nV,4),int(nO,4),int(nV,4))
  call gpu_allocate(Y_ovov,int(nO,4),int(nV,4),int(nO,4),int(nV,4))

  do a=1, nV
    call gpu_set_stream(blas_handle, stream(a))
    do j=1,nO
      call gpu_dgeam(blas_handle, 'N', 'N', int(nO,4), int(nV,4), 1.d0 &
        , t2%f(1,j,1,a), int(nO*nO,4), 0.d0, t2%f(1,j,1,a), int(nO*nO,4) &
        , Y_ovov%f(1,a,j,1), int(nO*nV*nO,4))
    enddo
    do beta=1, nV
      call gpu_dgeam(blas_handle, 'T', 'T', int(nO,4), int(nO,4), 0.5d0 &
        , K1%f(1,a,1,beta), int(nO*nV,4), 0.d0, K1%f(1,a,1,beta), int(nO*nV,4) &
        , X_ovov%f(1,a,1,beta), int(nO*nV,4))
    enddo
  enddo
  call gpu_set_stream(blas_handle, gpu_default_stream)

  call gpu_synchronize()

  call gpu_dgemm(blas_handle, 'T','N',int(nO*nV,4),int(nO*nV,4),int(nO*nV,4), &
            -1d0, X_ovov%f(1,1,1,1), size(X_ovov%f,1) * size(X_ovov%f,2), &
                  Y_ovov%f(1,1,1,1), size(Y_ovov%f,1) * size(Y_ovov%f,2), &
             1d0, Z_ovov%f(1,1,1,1), size(Z_ovov%f,1) * size(Z_ovov%f,2))

  call gpu_synchronize()

  do gam=1, nV
    call gpu_set_stream(blas_handle, stream(gam))
    do j=1,nO
      call gpu_dgeam(blas_handle, 'N', 'N', int(nO,4), int(nV,4), 1.d0 &
        , r2%f(1,j,1,gam), int(nO*nO,4), 1.d0, Z_ovov%f(1,1,j,gam), int(nO,4) &
        , r2%f(1,j,1,gam), int(nO*nO,4))
      call gpu_dgeam(blas_handle, 'N', 'N', int(nO,4), int(nV,4), 1.d0 &
        , K1%f(1,1,j,gam), int(nO,4), 0.d0, K1%f(1,1,j,gam), int(nO,4) &
        , X_ovov%f(1,gam,j,1), int(nO*nV*nO,4))
      call gpu_dgeam(blas_handle, 'N', 'N', int(nO,4), int(nO,4), 1.d0 &
        , t2%f(1,j,1,gam), int(nO*nO,4), 0.d0, t2%f(1,j,1,gam), int(nO*nO,4) &
        , Y_ovov%f(1,gam,j,1), int(nO*nV*nO,4))
    enddo
    do beta=1, nV
      call gpu_dgeam(blas_handle, 'N', 'T', int(nO,4), int(nO,4), 1.d0 &
        , r2%f(1,1,beta,gam), int(nO,4), 1.d0, Z_ovov%f(1,gam,1,beta) &
        , int(nO*nV,4), r2%f(1,1,beta,gam), int(nO,4))
    enddo
  enddo

  call gpu_set_stream(blas_handle, gpu_default_stream)

  call gpu_deallocate(K1)

  call gpu_dgemm(blas_handle, 'N','N',int(nO*nV,4),int(nO*nV,4),int(nO*nV,4), &
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
      call gpu_dgeam(blas_handle, 'N', 'N', int(nO,4), int(nV,4), 1.d0 &
        , r2%f(1,j,1,gam), int(nO*nO,4), -1.d0, Z_ovov%f(1,gam,j,1) &
        , int(nO*nV*nO,4), r2%f(1,j,1,gam), int(nO*nO,4))
    enddo
    do beta = 1, nV
      call gpu_dgeam(blas_handle, 'N', 'T', int(nO,4), int(nO,4), -1.d0 &
        , r2%f(1,1,beta,gam), int(nO,4), 1.d0, Z_ovov%f(1,beta,1,gam) &
        , int(nO*nV,4), r2%f(1,1,beta,gam), int(nO,4))
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
          max_r2_local = max(dabs(r2%f(i,j,a,b)), max_r2_local)
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
 use sizet
  BEGIN_DOC
  ! Computes A1(u,v,i,j), an intermediate of the T2 residual,
  ! the GPU/Cholesky counterpart of compute_A1, called from
  ! compute_r2_space_chol. A1 is built from d_cc_space_v_oooo,
  ! d_cc_space_v_vooo and d_cc_space_v_vvoo contracted with t1
  ! and tau, as sketched in the in-code comments; the v_ovoo
  ! term of compute_A1 is formed from d_cc_space_v_vooo with a
  ! transposed index.
  ! The BLAS formulation reproduces the four-term formula of
  ! compute_A1 (author-confirmed 2026-10-31).
  END_DOC
  use gpu
  implicit none

  integer, intent(in)            :: nO,nV
  type(gpu_double2), intent(in)  :: t1
  type(gpu_double4), intent(in)  :: t2, tau
  type(gpu_double4), intent(in)  :: d_cc_space_v_vooo, d_cc_space_v_oooo, d_cc_space_v_vvoo
  type(gpu_double4), intent(out) :: A1

  integer(size_t) :: i,j
  integer         :: a,tmp_a,b,k,l,c,d,tmp_c,tmp_d,u,v,beta

  type(gpu_double4) :: Y_oooo
  call gpu_allocate(Y_oooo,int(nO,4),int(nO,4),int(nO,4),int(nO,4))

  ! A1(u,v,i,j) = cc_space_v_oooo(u,v,i,j)
  ! A1(u,v,i,j) += cc_space_v_ovoo(u,a,i,j) * t1(v,a) &

  call gpu_dgemm(blas_handle, 'N','N', int(nO,4), int(nO*nO*nO,4), int(nV,4), &
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
      call gpu_dgeam(blas_handle, 'N', 'T', int(nO,4), int(nO,4), 1.d0 &
        , d_cc_space_v_oooo%f(1,1,i,j), int(nO,4), 1.d0, Y_oooo%f(1,1,j,i) &
        , int(nO,4), A1%f(1,1,i,j), int(nO,4))
    enddo
    call gpu_dgeam(blas_handle, 'N', 'N', int(nO,4), int(nO*nO,4), 1.d0 &
      , A1%f(1,1,1,j), int(nO,4), 1.d0, Y_oooo%f(1,1,1,j), int(nO,4) &
      , A1%f(1,1,1,j), int(nO,4))
  enddo

  call gpu_set_stream(blas_handle, gpu_default_stream)
  do i=1, nO
    call gpu_stream_destroy(stream(i))
  enddo

  call gpu_deallocate(Y_oooo)

  ! A1(u,v,i,j) += cc_space_v_vvoo(a,b,i,j) * tau(u,v,a,b)
  call gpu_dgemm(blas_handle, 'N', 'N', int(nO*nO,4), int(nO*nO,4) &
    , int(nV*nV,4), 1d0, tau%f(1,1,1,1), size(tau%f,1) * size(tau%f,2) &
    , d_cc_space_v_vvoo%f(1,1,1,1) &
    , size(d_cc_space_v_vvoo%f,1) * size(d_cc_space_v_vvoo%f,2), 1d0 &
    , A1%f(1,1,1,1), size(A1%f,1) * size(A1%f,2))
  call gpu_synchronize()

end

! g_occ

subroutine compute_g_occ_chol(nO,nV,t1,t2,H_oo, &
  d_cc_space_f_vo, d_cc_space_v_ov_chol, d_cc_space_v_oo_chol, d_cc_space_v_ovoo, g_occ)
  BEGIN_DOC
  ! Computes g_occ(u,i), an occupied-space intermediate of the
  ! T2 residual, the GPU/Cholesky counterpart of compute_g_occ,
  ! called from compute_r2_space_chol. g_occ is initialized to
  ! H_oo and updated with the t1 and d_cc_space_f_vo product, a
  ! t1 contraction over the Cholesky vectors d_cc_space_v_ov_chol
  ! and d_cc_space_v_oo_chol, and a d_cc_space_v_ovoo contraction
  ! with t1.
  ! The BLAS sequence reproduces the compute_g_occ formula
  ! H_oo(u,i) + sum f_vo(a,i) * t1(u,a) + sum (2 *
  ! v_ovoo(u,a,i,j) - v_ovoo(u,a,j,i)) * t1(j,a)
  ! (author-confirmed 2026-10-31).
  END_DOC
  use gpu

  implicit none

  integer, intent(in)           :: nO,nV
  type(gpu_double2), intent(in)  :: t1, H_oo, d_cc_space_f_vo
  type(gpu_double3), intent(in)  :: d_cc_space_v_ov_chol, d_cc_space_v_oo_chol
  type(gpu_double4), intent(in)  :: t2, d_cc_space_v_ovoo
  type(gpu_double2), intent(out) :: g_occ

  call gpu_copy(H_oo, g_occ)

  call gpu_dgemm(blas_handle, 'N','N',int(nO,4),int(nO,4),int(nV,4), &
             1d0, t1%f(1,1), size(t1%f,1), &
                  d_cc_space_f_vo%f(1,1), size(d_cc_space_f_vo%f,1), &
             1d0, g_occ%f(1,1), size(g_occ%f,1))

  type(gpu_double1) :: X
  call gpu_allocate(X,int(cholesky_mo_num,4))

  call gpu_dgemv(blas_handle, 'N',int(cholesky_mo_num,4),int(nO*nV,4),2.d0, &
    d_cc_space_v_ov_chol%f(1,1,1), int(cholesky_mo_num,4), &
    t1%f(1,1), 1, 0.d0, X%f(1), 1)

  call gpu_dgemv(blas_handle, 'T',int(cholesky_mo_num,4),int(nO*nO,4),1.d0, &
    d_cc_space_v_oo_chol%f(1,1,1), int(cholesky_mo_num,4), &
    X%f(1), 1, 1.d0, g_occ%f(1,1), 1)

  call gpu_dgemv(blas_handle, 'T',int(nO*nV,4),int(nO*nO,4),-1.d0, &
    d_cc_space_v_ovoo%f(1,1,1,1), int(nO*nV,4), &
    t1%f(1,1), 1, 1.d0, g_occ%f(1,1), 1)

  call gpu_synchronize()
  call gpu_deallocate(X)

end

! g_vir

subroutine compute_g_vir_chol(nO,nV,t1,t2,H_vv,d_cc_space_f_vo, &
  d_cc_space_v_ov_chol, d_cc_space_v_vv_chol, g_vir)
 use sizet
  BEGIN_DOC
  ! Computes g_vir(a,beta), a virtual-space intermediate of the
  ! T2 residual, the GPU/Cholesky counterpart of compute_g_vir,
  ! called from compute_r2_space_chol. g_vir is built from
  ! -d_cc_space_f_vo : t1, t1 contractions over the Cholesky
  ! vectors d_cc_space_v_ov_chol and d_cc_space_v_vv_chol, and
  ! H_vv.
  ! The BLAS sequence reproduces the compute_g_vir formula
  ! (H_vv - f_vo : t1 + the v_vvvo terms); the pre-GPU signed
  ! transpose of tmp_vo (commit fa6d1419~1) is one gpu_dgeam
  ! (author-confirmed 2026-10-31).
  END_DOC
  use gpu

  implicit none

  integer, intent(in)           :: nO,nV
  type(gpu_double2), intent(in)  :: t1, H_vv, d_cc_space_f_vo
  type(gpu_double3), intent(in)  :: d_cc_space_v_ov_chol, d_cc_space_v_vv_chol
  type(gpu_double4), intent(in)  :: t2
  type(gpu_double2), intent(out) :: g_vir

  integer(size_t) :: i
  integer         :: a,tmp_a,b,k,l,c,d,tmp_c,tmp_d,j,u,v,beta,gam

  type(gpu_stream) :: stream(max(nO,4))

  do i=1,max(nO,4)
    call gpu_stream_create(stream(i))
  enddo

  call gpu_set_stream(blas_handle, stream(1))
  call gpu_dgemm(blas_handle, 'N','N',int(nV,4),int(nV,4),int(nO,4), &
             -1d0, d_cc_space_f_vo%f(1,1) , size(d_cc_space_f_vo%f,1), &
                   t1%f(1,1)   , size(t1%f,1), &
              0d0, g_vir%f(1,1), size(g_vir%f,1))

  type(gpu_double1) :: tmp_k
  type(gpu_double3) :: tmp_vo, tmp_vo2

  call gpu_allocate(tmp_k,int(cholesky_mo_num,4))

  call gpu_set_stream(blas_handle, stream(2))
  call gpu_dgemm(blas_handle, 'N', 'N', int(cholesky_mo_num,4), 1, &
    int(nO*nV,4), 1.d0 &
    , d_cc_space_v_ov_chol%f(1,1,1), int(cholesky_mo_num,4), t1%f(1,1), nO*nV, &
    0.d0 &
    , tmp_k%f(1), int(cholesky_mo_num,4))

  call gpu_dgemm(blas_handle, 'T', 'N', int(nV*nV,4), 1, &
    int(cholesky_mo_num,4), 2.d0 &
    , d_cc_space_v_vv_chol%f(1,1,1), int(cholesky_mo_num,4), tmp_k%f(1) &
    , int(cholesky_mo_num,4), 1.d0, g_vir%f(1,1), int(nV*nV,4))

  call gpu_set_stream(blas_handle, stream(3))
  call gpu_allocate(tmp_vo,int(cholesky_mo_num,4),int(nV,4),int(nO,4))

  call gpu_dgemm(blas_handle, 'N', 'T', int(cholesky_mo_num*nV,4), int(nO,4) &
    , int(nV,4), 1.d0, d_cc_space_v_vv_chol%f(1,1,1) &
    , int(cholesky_mo_num*nV,4), t1%f(1,1), nO, 0.d0, tmp_vo%f(1,1,1) &
    , int(cholesky_mo_num*nV,4))

  call gpu_allocate(tmp_vo2,int(cholesky_mo_num,4),int(nO,4),int(nV,4))

  call gpu_synchronize()
  call gpu_deallocate(tmp_k)

  do i=1,nO
    call gpu_set_stream(blas_handle, stream(i))
    call gpu_dgeam(blas_handle, 'N', 'N', int(cholesky_mo_num,4), int(nV,4), &
      -1.d0 &
      , tmp_vo%f(1,1,i), int(cholesky_mo_num,4), 0.d0, tmp_vo%f(1,1,i) &
      , int(cholesky_mo_num,4), tmp_vo2%f(1,i,1), int(cholesky_mo_num*nO,4))
  enddo

  call gpu_set_stream(blas_handle, gpu_default_stream)

  do i=1,max(nO,4)
    call gpu_stream_destroy(stream(i))
  enddo
  call gpu_deallocate(tmp_vo)

  call gpu_dgeam(blas_handle, 'N', 'N', int(nV,4), int(nV,4), 1.d0 &
    , g_vir%f(1,1), int(nV,4), 1.d0, H_vv%f(1,1), int(nV,4), g_vir%f(1,1) &
    , int(nV,4))

  call gpu_dgemm(blas_handle, 'T', 'N', int(nV,4), int(nV,4) &
    , int(nO*cholesky_mo_num,4), 1.d0, d_cc_space_v_ov_chol%f(1,1,1) &
    , int(cholesky_mo_num*nO,4), tmp_vo2%f(1,1,1), int(cholesky_mo_num*nO,4), &
    1.d0 &
    , g_vir%f(1,1), int(nV,4))

  call gpu_synchronize()
  call gpu_deallocate(tmp_vo2)

end

! J1

subroutine compute_J1_chol(nO,nV,t1,t2,v_ovvo,v_ovoo,v_vvoo,d_cc_space_v_vo_chol,d_cc_space_v_vv_chol,J1)
 use sizet
  BEGIN_DOC
  ! Computes J1(u,a,beta,i), an intermediate of the T2 residual,
  ! the GPU/Cholesky counterpart of compute_J1, called from
  ! compute_r2_space_chol. J1 is initialized to d_cc_space_v_ovvo
  ! and updated from d_cc_space_v_ovoo, d_cc_space_v_vvoo, t1,
  ! t2 and the Cholesky vectors d_cc_space_v_vo_chol and
  ! d_cc_space_v_vv_chol.
  ! The BLAS sequence reproduces the compute_J1 formula; the
  ! pre-GPU explicit expressions for Y_ovov and X_voov (commit
  ! fa6d1419~1) are gpu_dgeam calls on t2 and v_vvoo slices
  ! (author-confirmed 2026-10-31).
  END_DOC
  use gpu
  implicit none

  integer, intent(in)            :: nO,nV
  type(gpu_double2), intent(in)  :: t1
  type(gpu_double4), intent(in)  :: t2, v_ovvo, v_ovoo, v_vvoo
  type(gpu_double3), intent(in) :: d_cc_space_v_vo_chol,d_cc_space_v_vv_chol
  type(gpu_double4), intent(out) :: J1

  integer(size_t) :: a,b,i,j,beta
  integer         :: tmp_a,k,l,c,d,tmp_c,tmp_d,u,v,gam

  type(gpu_double4) :: X_ovoo, Y_ovov

  call gpu_allocate(X_ovoo,int(nO,4),int(nV,4),int(nO,4),int(nO,4))

  type(gpu_stream) :: stream(nV)


  do i=1,nO
    call gpu_stream_create(stream(i))
  enddo

  do j = 1, nO
    call gpu_set_stream(blas_handle, stream(j))
    do i = 1, nO
      call gpu_dgeam(blas_handle, 'N', 'N', int(nO,4), int(nV,4), 1.d0 &
        , v_ovoo%f(1,1,j,i), int(nO,4), 0.d0, X_ovoo%f(1,1,i,j), int(nO,4) &
        , X_ovoo%f(1,1,i,j), int(nO,4))
    enddo
  enddo

  call gpu_set_stream(blas_handle, gpu_default_stream)

  do i=1,nO
    call gpu_stream_destroy(stream(i))
  enddo

  call gpu_allocate(Y_ovov,int(nO,4),int(nV,4),int(nO,4),int(nV,4))

  call gpu_dgemm(blas_handle, 'N', 'N', int(nO*nV*nO,4), int(nV,4), int(nO,4) &
    , -1d0, X_ovoo%f(1,1,1,1) &
    , size(X_ovoo%f,1) * size(X_ovoo%f,2) * size(X_ovoo%f,3), t1%f(1,1) &
    , size(t1%f,1), 0d0, Y_ovov%f(1,1,1,1) &
    , size(Y_ovov%f,1) * size(Y_ovov%f,2) * size(Y_ovov%f,3))


  call gpu_copy(v_ovvo, J1)

  call gpu_synchronize()

  do a=1,nV
    call gpu_stream_create(stream(a))
  enddo

  do i = 1, nO
    do beta = 1, nV
      call gpu_set_stream(blas_handle, stream(beta))
      call gpu_dgeam(blas_handle, 'N', 'N', int(nO,4), int(nV,4), 1.d0 &
        , J1%f(1,1,beta,i), int(nO,4), 1.d0, Y_ovov%f(1,1,i,beta), int(nO,4) &
        , J1%f(1,1,beta,i), int(nO,4))
    enddo
  enddo

  call gpu_allocate(tmp_cc,int(cholesky_mo_num,4),int(nV,4),int(nO,4))
  call gpu_allocate(J1_tmp,int(nV,4),int(nO,4),int(nV,4),int(nO,4))

  call gpu_set_stream(blas_handle, gpu_default_stream)

  type(gpu_double4) :: J1_tmp
  type(gpu_double3) :: tmp_cc

  call gpu_dgemm(blas_handle, 'N', 'T', int(cholesky_mo_num*nV,4), int(nO,4) &
    , int(nV,4), 1.d0, d_cc_space_v_vv_chol%f(1,1,1) &
    , int(cholesky_mo_num*nV,4), t1%f(1,1), nO, 0.d0, tmp_cc%f(1,1,1) &
    , int(cholesky_mo_num*nV,4))

  call gpu_dgemm(blas_handle, 'T', 'N', int(nV*nO,4), int(nV*nO,4) &
    , int(cholesky_mo_num,4), 1.d0, tmp_cc%f(1,1,1), int(cholesky_mo_num,4) &
    , d_cc_space_v_vo_chol%f(1,1,1), int(cholesky_mo_num,4), 0.d0, &
    J1_tmp%f(1,1,1,1) &
    , int(nV*nO,4))


  call gpu_deallocate(X_ovoo)

  call gpu_synchronize()
  call gpu_deallocate(tmp_cc)

  do i = 1, nO
    do a = 1, nV
      call gpu_set_stream(blas_handle, stream(a))
      call gpu_dgeam(blas_handle, 'N', 'T', int(nO,4), int(nV,4), 1.d0 &
        , J1%f(1,a,1,i), int(nO*nV,4), 1.d0, J1_tmp%f(1,1,a,i), int(nV,4) &
        , J1%f(1,a,1,i), int(nO*nV,4))
    enddo
  enddo

  type(gpu_double4) :: X_voov, Z_ovvo

  call gpu_allocate(X_voov,int(nV,4),int(nO,4),int(nO,4),int(nV,4))
  call gpu_allocate(Z_ovvo,int(nO,4),int(nV,4),int(nV,4),int(nO,4))

  do j = 1, nO
    do beta = 1, nV
      call gpu_set_stream(blas_handle, stream(beta))
      call gpu_dgeam(blas_handle, 'N', 'N', int(nO,4), int(nV,4), 0.5d0 &
        , t2%f(1,j,1,beta), int(nO*nO,4), t1%f(j,beta), t1%f(1,1), int(nO,4) &
        , Y_ovov%f(1,beta,j,1), int(nO*nV*nO,4))
    enddo
  enddo

  do b = 1, nV
    call gpu_set_stream(blas_handle, stream(b))
    call gpu_dgeam(blas_handle, 'N', 'N', int(nV,4), int(nO*nO,4), 1.d0 &
      , v_vvoo%f(1,b,1,1), int(nV*nV,4), 0.d0, X_voov%f(1,1,1,b), int(nV,4) &
      , X_voov%f(1,1,1,b), int(nV,4))
  enddo

  call gpu_set_stream(blas_handle, gpu_default_stream)

  call gpu_synchronize()
  call gpu_deallocate(J1_tmp)

  call gpu_dgemm(blas_handle, 'N','T',int(nO*nV,4),int(nV*nO,4),int(nO*nV,4), &
             -1d0, Y_ovov%f(1,1,1,1), size(Y_ovov%f,1) * size(Y_ovov%f,2), &
                   X_voov%f(1,1,1,1), size(X_voov%f,1) * size(X_voov%f,2), &
              0d0, Z_ovvo%f(1,1,1,1), size(Z_ovvo%f,1) * size(Z_ovvo%f,2))

  call gpu_synchronize()

  do i = 1, nO
    do a = 1, nV
      call gpu_set_stream(blas_handle, stream(a))
      call gpu_dgeam(blas_handle, 'N', 'N', int(nO,4), int(nV,4), 1.d0 &
        , J1%f(1,a,1,i), int(nO*nV,4), 1.d0, Z_ovvo%f(1,1,a,i), int(nO,4) &
        , J1%f(1,a,1,i), int(nO*nV,4))
    enddo
  enddo

  type(gpu_double4) :: X_ovvo, Y_vovo
  call gpu_allocate(Y_vovo,int(nV,4),int(nO,4),int(nV,4),int(nO,4))

  do j = 1, nO
    do i = 1, nO
      call gpu_set_stream(blas_handle, stream(i))
      call gpu_dgeam(blas_handle, 'N', 'T', int(nV,4), int(nV,4), 1.d0 &
        , v_vvoo%f(1,1,i,j), int(nV,4), -0.5d0, v_vvoo%f(1,1,i,j), int(nV,4) &
        , Y_vovo%f(1,i,1,j), int(nO*nV,4))
    enddo
  enddo

  call gpu_allocate(X_ovvo,int(nO,4),int(nV,4),int(nV,4),int(nO,4))

  do j = 1, nO
    do b = 1, nV
      call gpu_set_stream(blas_handle, stream(b))
      call gpu_dgeam(blas_handle, 'N', 'N', int(nO,4), int(nV,4), 1.d0 &
        , t2%f(1,j,1,b), int(nO*nO,4), 0.d0, t2%f(1,j,1,b), int(nO*nO,4) &
        , X_ovvo%f(1,1,b,j), int(nO,4))
    enddo
  enddo

  call gpu_set_stream(blas_handle, gpu_default_stream)
  call gpu_synchronize()
  call gpu_deallocate(X_voov)

  call gpu_dgemm(blas_handle, 'N','T',int(nO*nV,4),int(nV*nO,4),int(nV*nO,4), &
             1d0, X_ovvo%f(1,1,1,1), size(X_ovvo%f,1) * size(X_ovvo%f,2), &
                  Y_vovo%f(1,1,1,1), size(Y_vovo%f,1) * size(Y_vovo%f,2), &
             0d0, Z_ovvo%f(1,1,1,1), size(Z_ovvo%f,1) * size(Z_ovvo%f,2))

  call gpu_synchronize()

  do i = 1, nO
    do beta = 1, nV
      call gpu_set_stream(blas_handle, stream(beta))
      call gpu_dgeam(blas_handle, 'N', 'N', int(nO,4), int(nV,4), 1.d0 &
        , J1%f(1,1,beta,i), int(nO,4), 1.d0, Z_ovvo%f(1,beta,1,i) &
        , int(nO*nV,4), J1%f(1,1,beta,i), int(nO,4))
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
 use sizet
  BEGIN_DOC
  ! Computes K1(u,a,i,beta), an intermediate of the T2 residual,
  ! the GPU/Cholesky counterpart of compute_K1, called from
  ! compute_r2_space_chol. K1 is initialized to d_cc_space_v_ovov
  ! and updated from d_cc_space_v_ovoo, d_cc_space_v_vvoo, t1,
  ! t2 and the Cholesky vectors d_cc_space_v_ov_chol and
  ! d_cc_space_v_vv_chol.
  ! The BLAS sequence reproduces the compute_K1 formula; the
  ! pre-GPU explicit loops for Y and for the K1tmp accumulation
  ! (commit fa6d1419~1) are gpu_dgeam calls that a standard
  ! BLAS reading of the arguments does not reproduce
  ! (author-confirmed 2026-10-31).
  END_DOC
  use gpu

  implicit none

  integer, intent(in)            :: nO,nV
  type(gpu_double2), intent(in)  :: t1
  type(gpu_double4), intent(in)  :: t2, v_vvoo, v_ovov, v_ovoo
  type(gpu_double3), intent(in)  :: d_cc_space_v_ov_chol, d_cc_space_v_vv_chol
  type(gpu_double4), intent(out) :: K1

  type(gpu_double4) :: X, Y, Z

  integer(size_t) :: a,i,beta
  integer         :: tmp_a,b,k,l,c,d,tmp_c,tmp_d,j,u,v,gam


  call gpu_copy(v_ovov, K1)

  type(gpu_stream) :: stream(nV)
  do a = 1, nV
    call gpu_stream_create(stream(a))
  enddo

  call gpu_allocate(X,int(nV,4),int(nO,4),int(nV,4),int(nO,4))
  call gpu_allocate(Y,int(nO,4),int(nV,4),int(nV,4),int(nO,4))

  do a = 1, nV
    call gpu_set_stream(blas_handle, stream(a))
    do i = 1, nO
      call gpu_dgeam(blas_handle, 'N', 'N', int(nV,4), int(nO,4), -1.d0 &
        , v_vvoo%f(1,a,i,1), int(nV*nV*nO,4), 0.d0, v_vvoo%f(1,a,i,1) &
        , int(nV*nV*nO,4), X%f(1,1,a,i), int(nV,4))
      call gpu_dgeam(blas_handle, 'N', 'N', int(nO,4), int(nV,4), 0.5d0 &
        , t2%f(1,i,1,a), int(nO*nO,4), t1%f(i,a), t1%f(1,1), int(nO,4) &
        , Y%f(1,a,1,i), int(nO*nV,4))
    enddo
  enddo

  call gpu_set_stream(blas_handle, gpu_default_stream)

  call gpu_dgemm(blas_handle, 'N', 'N', int(nO*nV*nO,4), int(nV,4), int(nO,4) &
    , -1d0, v_ovoo%f(1,1,1,1) &
    , size(v_ovoo%f,1) * size(v_ovoo%f,2) * size(v_ovoo%f,3), t1%f(1,1) &
    , size(t1%f,1), 1d0, K1%f(1,1,1,1) &
    , size(K1%f,1) * size(K1%f,2) * size(K1%f,3))

  type(gpu_double4) :: K1tmp
  type(gpu_double3) :: t1v

  call gpu_allocate(t1v,int(cholesky_mo_num,4),int(nO,4),int(nO,4))

  call gpu_dgemm(blas_handle, 'N', 'T', int(cholesky_mo_num*nO,4), int(nO,4) &
    , int(nV,4), 1.d0, d_cc_space_v_ov_chol%f(1,1,1) &
    , int(cholesky_mo_num*nO,4), t1%f(1,1), nO, 0.d0, t1v%f(1,1,1) &
    , int(cholesky_mo_num*nO,4))

  call gpu_allocate(K1tmp,int(nO,4),int(nO,4),int(nV,4),int(nV,4))

  call gpu_dgemm(blas_handle, 'T', 'N', int(nO*nO,4), int(nV*nV,4) &
    , int(cholesky_mo_num,4), 1.d0, t1v%f(1,1,1), int(cholesky_mo_num,4) &
    , d_cc_space_v_vv_chol%f(1,1,1), int(cholesky_mo_num,4), 0.d0, &
    K1tmp%f(1,1,1,1) &
    , int(nO*nO,4))

  call gpu_allocate(Z,int(nO,4),int(nV,4),int(nV,4),int(nO,4))
  call gpu_synchronize()

  ! Y(u,beta,b,j) * X(b,j,a,i) = Z(u,beta,a,i)
  call gpu_dgemm(blas_handle, 'N','N',int(nV*nO,4),int(nO*nV,4),int(nV*nO,4), &
             1d0, Y%f(1,1,1,1), size(Y%f,1) * size(Y%f,2), &
                  X%f(1,1,1,1), size(X%f,1) * size(X%f,2), &
             0d0, Z%f(1,1,1,1), size(Z%f,1) * size(Z%f,2))

  call gpu_synchronize()
  call gpu_deallocate(t1v)

  do beta = 1, nV
   call gpu_set_stream(blas_handle, stream(beta))
   do i = 1, nO
      call gpu_dgeam(blas_handle, 'N', 'N', int(nO,4), int(nV,4), 1.d0 &
        , K1%f(1,1,i,beta), int(nO,4), 1.d0, K1tmp%f(1,i,1,beta), int(nO*nO,4) &
        , K1%f(1,1,i,beta), int(nO,4))
      call gpu_dgeam(blas_handle, 'N', 'N', int(nO,4), int(nV,4), 1.d0 &
        , K1%f(1,1,i,beta), int(nO,4), 1.d0, Z%f(1,beta,1,i), int(nO*nV,4) &
        , K1%f(1,1,i,beta), int(nO,4))
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
