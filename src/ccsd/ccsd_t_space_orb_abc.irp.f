! Main

subroutine ccsd_par_t_space_v3(nO,nV,t1,t2,f_o,f_v,v_vvvo,v_vvoo,v_vooo,energy)

  implicit none

  integer, intent(in)           :: nO,nV
  double precision, intent(in)  :: t1(nO,nV), f_o(nO), f_v(nV)
  double precision, intent(in)  :: t2(nO,nO,nV,nV)
  double precision, intent(in)  :: v_vvvo(nV,nV,nV,nO), v_vvoo(nV,nV,nO,nO), v_vooo(nV,nO,nO,nO)
  double precision, intent(out) :: energy

  double precision, allocatable :: W(:,:,:,:,:,:)
  double precision, allocatable :: V(:,:,:,:,:,:)
  double precision, allocatable :: W_abc(:,:,:), V_abc(:,:,:)
  double precision, allocatable :: W_cab(:,:,:), W_cba(:,:,:)
  double precision, allocatable :: W_bca(:,:,:)
  double precision, allocatable :: X_vovv(:,:,:,:), X_ooov(:,:,:,:), X_oovv(:,:,:,:)
  double precision, allocatable :: T_voov(:,:,:,:), T_oovv(:,:,:,:)
  integer                       :: i,j,k,l,a,b,c,d
  double precision              :: e,ta,tb, delta, delta_abc

  allocate(X_vovv(nV,nO,nV,nV), X_ooov(nO,nO,nO,nV), X_oovv(nO,nO,nV,nV))
  allocate(T_voov(nV,nO,nO,nV),T_oovv(nO,nO,nV,nV))

  call set_multiple_levels_omp(.False.)

  ! Temporary arrays
  !$OMP PARALLEL &
  !$OMP SHARED(nO,nV,T_voov,T_oovv,X_vovv,X_ooov,X_oovv, &
  !$OMP t1,t2,v_vvvo,v_vooo,v_vvoo) &
  !$OMP PRIVATE(a,b,c,d,i,j,k,l) &
  !$OMP DEFAULT(NONE)

  !v_vvvo(b,a,d,i) * t2(k,j,c,d) &
  !X_vovv(d,i,b,a,i) * T_voov(d,j,c,k)

  !$OMP DO collapse(3)
  do i = 1, nO
    do a = 1, nV
      do b = 1, nV
        do d = 1, nV
          X_vovv(d,i,b,a) = v_vvvo(b,a,d,i)
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO nowait

  !$OMP DO collapse(3)
  do j = 1, nO
    do k = 1, nO
      do c = 1, nV
        do d = 1, nV
          T_voov(d,k,j,c) = t2(k,j,c,d)
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO nowait

  !v_vooo(c,j,k,l) * t2(i,l,a,b) &
  !X_ooov(l,j,k,c) * T_oovv(l,i,a,b) &

  !$OMP DO collapse(3)
  do c = 1, nV
    do k = 1, nO
      do j = 1, nO
        do l = 1, nO
           X_ooov(l,j,k,c) = v_vooo(c,j,k,l)
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO nowait

  !$OMP DO collapse(3)
  do i = 1, nO
    do b = 1, nV
      do a = 1, nV
        do l = 1, nO
          T_oovv(l,i,a,b) = t2(i,l,a,b)
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO nowait

  !X_oovv(j,k,b,c) * T1_vo(a,i) &

  !$OMP DO collapse(3)
  do c = 1, nV
    do b = 1, nV
      do j = 1, nO
        do k = 1, nO
          X_oovv(j,k,b,c) = v_vvoo(b,c,j,k)
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO nowait

  !$OMP END PARALLEL

  energy = 0d0
  !$OMP PARALLEL                                               &
  !$OMP PRIVATE(a,b,c,W_abc,W_cab,W_bca,W_cba,V_abc)           &
  !$OMP PRIVATE(i,j,k,e,delta,delta_abc)                       &
  !$OMP DEFAULT(SHARED)
  allocate(W_abc(nO,nO,nO), W_cab(nO,nO,nO), V_abc(nO,nO,nO), &
           W_bca(nO,nO,nO), W_cba(nO,nO,nO) )
  !$OMP DO
  do c = 1, nV
    do b = 1, nV
      do a = 1, nV
        e = 0d0
        delta_abc = f_v(a) + f_v(b) + f_v(c)
        call form_w_abc(nO,nV,a,b,c,T_voov,T_oovv,X_vovv,X_ooov,W_abc,W_cba,W_bca,W_cab)
        call form_v_abc(nO,nV,a,b,c,t1,X_oovv,W_abc,V_abc,W_cba)
        do i = 1, nO
          do j = 1, nO
            do k = 1, nO
              delta = 1d0 / (f_o(i) + f_o(j) + f_o(k) - delta_abc)
              e = e + (4d0 * W_abc(i,j,k) + W_bca(i,j,k) + W_cab(i,j,k))&
                  * V_abc(i,j,k) * delta
            enddo
          enddo
        enddo
        !$OMP CRITICAL
        energy = energy + e
        !$OMP END CRITICAL
      enddo
    enddo
  enddo
  !$OMP END DO

  deallocate(W_abc,V_abc,W_cab,W_bca,W_cba)
  !$OMP END PARALLEL

  energy = energy / 3.d0

  deallocate(X_vovv,X_ooov,T_voov,T_oovv)
end


subroutine form_w_abc(nO,nV,a,b,c,T_voov,T_oovv,X_vovv,X_ooov,W_abc,W_cba,W_bca,W_cab)

  implicit none

  integer, intent(in)           :: nO,nV,a,b,c
  double precision, intent(in)  :: T_voov(nV,nO,nO,nV), T_oovv(nO,nO,nV,nV)
  double precision, intent(in)  :: X_vovv(nV,nO,nV,nV), X_ooov(nO,nO,nO,nV)
  double precision, intent(out) :: W_abc(nO,nO,nO)
  double precision, intent(out) :: W_cba(nO,nO,nO)
  double precision, intent(out) :: W_bca(nO,nO,nO)
  double precision, intent(out) :: W_cab(nO,nO,nO)

  integer :: l,i,j,k,d
  double precision, allocatable, dimension(:,:,:,:) :: W_ikj
  double precision, allocatable :: X(:,:,:,:)

  allocate(W_ikj(nO,nO,nO,4))
  allocate(X(nV,nO,nO,3))

  do k=1,nO
    do i=1,nO
      do d=1,nV
        X(d,i,k,1) = T_voov(d,k,i,c)
        X(d,i,k,2) = T_voov(d,k,i,b)
        X(d,i,k,3) = T_voov(d,k,i,a)
      enddo
    enddo
  enddo

!   X_vovv(d,i,c,a) * T_voov(d,j,k,b) : i jk

  call dgemm('T','N', nO, nO*nO, nV, 1.d0, X_vovv(1,1,c,a), nV, T_voov(1,1,1,b), nV, 0.d0, W_abc, nO)
  call dgemm('T','N', nO, nO*nO, nV, 1.d0, X_vovv(1,1,a,c), nV, T_voov(1,1,1,b), nV, 0.d0, W_cba, nO)
  call dgemm('T','N', nO, nO*nO, nV, 1.d0, X_vovv(1,1,a,b), nV, T_voov(1,1,1,c), nV, 0.d0, W_bca, nO)
  call dgemm('T','N', nO, nO*nO, nV, 1.d0, X_vovv(1,1,b,c), nV, T_voov(1,1,1,a), nV, 0.d0, W_cab, nO)

!   T_voov(d,i,j,a) * X_vovv(d,k,b,c) : ij k

  call dgemm('T','N', nO*nO, nO, nV, 1.d0, T_voov(1,1,1,a), nV, X_vovv(1,1,b,c), nV, 1.d0, W_abc, nO*nO)
  call dgemm('T','N', nO*nO, nO, nV, 1.d0, T_voov(1,1,1,c), nV, X_vovv(1,1,b,a), nV, 1.d0, W_cba, nO*nO)
  call dgemm('T','N', nO*nO, nO, nV, 1.d0, T_voov(1,1,1,b), nV, X_vovv(1,1,c,a), nV, 1.d0, W_bca, nO*nO)
  call dgemm('T','N', nO*nO, nO, nV, 1.d0, T_voov(1,1,1,c), nV, X_vovv(1,1,a,b), nV, 1.d0, W_cab, nO*nO)


!   X_vovv(d,k,a,c) * T_voov(d,j,i,b) : k ji

  call dgemm('T','N', nO*nO, nO, nV, 1.d0, X(1,1,1,2), nV, X_vovv(1,1,a,c), nV, 1.d0, W_abc, nO*nO)
  call dgemm('T','N', nO*nO, nO, nV, 1.d0, X(1,1,1,2), nV, X_vovv(1,1,c,a), nV, 1.d0, W_cba, nO*nO)
  call dgemm('T','N', nO*nO, nO, nV, 1.d0, X(1,1,1,1), nV, X_vovv(1,1,b,a), nV, 1.d0, W_bca, nO*nO)
  call dgemm('T','N', nO*nO, nO, nV, 1.d0, X(1,1,1,3), nV, X_vovv(1,1,c,b), nV, 1.d0, W_cab, nO*nO)

!   X_vovv(d,i,b,a) * T_voov(d,k,j,c) : i kj

  call dgemm('T','N', nO, nO*nO, nV, 1.d0, X_vovv(1,1,b,a), nV, X(1,1,1,1), nV, 1.d0, W_abc, nO)
  call dgemm('T','N', nO, nO*nO, nV, 1.d0, X_vovv(1,1,b,c), nV, X(1,1,1,3), nV, 1.d0, W_cba, nO)
  call dgemm('T','N', nO, nO*nO, nV, 1.d0, X_vovv(1,1,c,b), nV, X(1,1,1,3), nV, 1.d0, W_bca, nO)
  call dgemm('T','N', nO, nO*nO, nV, 1.d0, X_vovv(1,1,a,c), nV, X(1,1,1,2), nV, 1.d0, W_cab, nO)

!  T_voov(d,k,i,c) * X_vovv(d,j,a,b) : ki j

  call dgemm('T','N', nO*nO, nO, nV, 1.d0, X(1,1,1,1), nV, X_vovv(1,1,a,b), nV, 0.d0, W_ikj(1,1,1,1), nO*nO)
  call dgemm('T','N', nO*nO, nO, nV, 1.d0, X(1,1,1,3), nV, X_vovv(1,1,c,b), nV, 0.d0, W_ikj(1,1,1,2), nO*nO)
  call dgemm('T','N', nO*nO, nO, nV, 1.d0, X(1,1,1,3), nV, X_vovv(1,1,b,c), nV, 0.d0, W_ikj(1,1,1,3), nO*nO)
  call dgemm('T','N', nO*nO, nO, nV, 1.d0, X(1,1,1,2), nV, X_vovv(1,1,c,a), nV, 0.d0, W_ikj(1,1,1,4), nO*nO)

!   T_voov(d,i,k,a) * X_vovv(d,j,c,b) : ik j
  call dgemm('T','N', nO*nO, nO, nV, 1.d0, T_voov(1,1,1,a), nV, X_vovv(1,1,c,b), nV, 1.d0, W_ikj(1,1,1,1), nO*nO)
  call dgemm('T','N', nO*nO, nO, nV, 1.d0, T_voov(1,1,1,c), nV, X_vovv(1,1,a,b), nV, 1.d0, W_ikj(1,1,1,2), nO*nO)
  call dgemm('T','N', nO*nO, nO, nV, 1.d0, T_voov(1,1,1,b), nV, X_vovv(1,1,a,c), nV, 1.d0, W_ikj(1,1,1,3), nO*nO)
  call dgemm('T','N', nO*nO, nO, nV, 1.d0, T_voov(1,1,1,c), nV, X_vovv(1,1,b,a), nV, 1.d0, W_ikj(1,1,1,4), nO*nO)

  deallocate(X)

  allocate(X(nO,nO,nO,3))

  do k=1,nO
    do j=1,nO
      do l=1,nO
        X(l,j,k,1) = X_ooov(l,k,j,b)
        X(l,j,k,2) = X_ooov(l,k,j,a)
        X(l,j,k,3) = X_ooov(l,k,j,c)
      enddo
    enddo
  enddo


!   - T_oovv(l,i,a,b) * X_ooov(l,j,k,c) : i jk
  call dgemm('T','N', nO, nO*nO, nO, -1.d0, T_oovv(1,1,a,b), nO, X_ooov(1,1,1,c), nO, 1.d0, W_abc, nO)
  call dgemm('T','N', nO, nO*nO, nO, -1.d0, T_oovv(1,1,c,b), nO, X_ooov(1,1,1,a), nO, 1.d0, W_cba, nO)
  call dgemm('T','N', nO, nO*nO, nO, -1.d0, T_oovv(1,1,b,c), nO, X_ooov(1,1,1,a), nO, 1.d0, W_bca, nO)
  call dgemm('T','N', nO, nO*nO, nO, -1.d0, T_oovv(1,1,c,a), nO, X_ooov(1,1,1,b), nO, 1.d0, W_cab, nO)

!   - T_oovv(l,i,a,c) * X_ooov(l,k,j,b) : i kj
  call dgemm('T','N', nO, nO*nO, nO, -1.d0, T_oovv(1,1,a,c), nO, X(1,1,1,1), nO, 1.d0, W_abc, nO)
  call dgemm('T','N', nO, nO*nO, nO, -1.d0, T_oovv(1,1,c,a), nO, X(1,1,1,1), nO, 1.d0, W_cba, nO)
  call dgemm('T','N', nO, nO*nO, nO, -1.d0, T_oovv(1,1,b,a), nO, X(1,1,1,3), nO, 1.d0, W_bca, nO)
  call dgemm('T','N', nO, nO*nO, nO, -1.d0, T_oovv(1,1,c,b), nO, X(1,1,1,2), nO, 1.d0, W_cab, nO)

!   - X_ooov(l,i,j,b) * T_oovv(l,k,c,a) : ij k
  call dgemm('T','N', nO*nO, nO, nO, -1.d0, X_ooov(1,1,1,b), nO, T_oovv(1,1,c,a), nO, 1.d0, W_abc, nO*nO)
  call dgemm('T','N', nO*nO, nO, nO, -1.d0, X_ooov(1,1,1,b), nO, T_oovv(1,1,a,c), nO, 1.d0, W_cba, nO*nO)
  call dgemm('T','N', nO*nO, nO, nO, -1.d0, X_ooov(1,1,1,c), nO, T_oovv(1,1,a,b), nO, 1.d0, W_bca, nO*nO)
  call dgemm('T','N', nO*nO, nO, nO, -1.d0, X_ooov(1,1,1,a), nO, T_oovv(1,1,b,c), nO, 1.d0, W_cab, nO*nO)

!   - X_ooov(l,j,i,a) * T_oovv(l,k,c,b) : ji k
  call dgemm('T','N', nO*nO, nO, nO, -1.d0, X(1,1,1,2), nO, T_oovv(1,1,c,b), nO, 1.d0, W_abc, nO*nO)
  call dgemm('T','N', nO*nO, nO, nO, -1.d0, X(1,1,1,3), nO, T_oovv(1,1,a,b), nO, 1.d0, W_cba, nO*nO)
  call dgemm('T','N', nO*nO, nO, nO, -1.d0, X(1,1,1,1), nO, T_oovv(1,1,a,c), nO, 1.d0, W_bca, nO*nO)
  call dgemm('T','N', nO*nO, nO, nO, -1.d0, X(1,1,1,3), nO, T_oovv(1,1,b,a), nO, 1.d0, W_cab, nO*nO)

!   - X_ooov(l,k,i,a) * T_oovv(l,j,b,c) : ki j
  call dgemm('T','N', nO*nO, nO, nO, -1.d0, X(1,1,1,2), nO, T_oovv(1,1,b,c), nO, 1.d0, W_ikj(1,1,1,1), nO*nO)
  call dgemm('T','N', nO*nO, nO, nO, -1.d0, X(1,1,1,3), nO, T_oovv(1,1,b,a), nO, 1.d0, W_ikj(1,1,1,2), nO*nO)
  call dgemm('T','N', nO*nO, nO, nO, -1.d0, X(1,1,1,1), nO, T_oovv(1,1,c,a), nO, 1.d0, W_ikj(1,1,1,3), nO*nO)
  call dgemm('T','N', nO*nO, nO, nO, -1.d0, X(1,1,1,3), nO, T_oovv(1,1,a,b), nO, 1.d0, W_ikj(1,1,1,4), nO*nO)

!   - X_ooov(l,i,k,c) * T_oovv(l,j,b,a) : ik j
  call dgemm('T','N', nO*nO, nO, nO, -1.d0, X_ooov(1,1,1,c), nO, T_oovv(1,1,b,a), nO, 1.d0, W_ikj(1,1,1,1), nO*nO)
  call dgemm('T','N', nO*nO, nO, nO, -1.d0, X_ooov(1,1,1,a), nO, T_oovv(1,1,b,c), nO, 1.d0, W_ikj(1,1,1,2), nO*nO)
  call dgemm('T','N', nO*nO, nO, nO, -1.d0, X_ooov(1,1,1,a), nO, T_oovv(1,1,c,b), nO, 1.d0, W_ikj(1,1,1,3), nO*nO)
  call dgemm('T','N', nO*nO, nO, nO, -1.d0, X_ooov(1,1,1,b), nO, T_oovv(1,1,a,c), nO, 1.d0, W_ikj(1,1,1,4), nO*nO)

  do k=1,nO
    do j=1,nO
      do i=1,nO
        W_abc(i,j,k) = W_abc(i,j,k) + W_ikj(i,k,j,1)
        W_cba(i,j,k) = W_cba(i,j,k) + W_ikj(i,k,j,2)
        W_bca(i,j,k) = W_bca(i,j,k) + W_ikj(i,k,j,3)
        W_cab(i,j,k) = W_cab(i,j,k) + W_ikj(i,k,j,4)
      enddo
    enddo
  enddo

  deallocate(X,W_ikj)
end


! V_abc

subroutine form_v_abc(nO,nV,a,b,c,T_ov,X_oovv,W_abc,V_abc,W_cba)

implicit none

  integer, intent(in)           :: nO,nV,a,b,c
  double precision, intent(in)  :: T_ov(nO,nV)
  double precision, intent(in)  :: X_oovv(nO,nO,nV,nV)
  double precision, intent(in)  :: W_abc(nO,nO,nO), W_cba(nO,nO,nO)
  double precision, intent(out) :: V_abc(nO,nO,nO)

  integer :: i,j,k

  do k = 1, nO
    do j = 1, nO
      do i = 1, nO
        V_abc(i,j,k) = W_abc(i,j,k) - W_cba(i,j,k) &
           + X_oovv(j,k,b,c) * T_ov(i,a) &
           + X_oovv(i,k,a,c) * T_ov(j,b) &
           + X_oovv(i,j,a,b) * T_ov(k,c) &
           - X_oovv(j,k,b,a) * T_ov(i,c) &
           - X_oovv(i,k,c,a) * T_ov(j,b) &
           - X_oovv(i,j,c,b) * T_ov(k,a)
      enddo
    enddo
  enddo

end

