! Dumb way

subroutine ccsd_par_t_space(nO,nV,t1,t2,energy)

  implicit none

  integer, intent(in)           :: nO,nV
  double precision, intent(in)  :: t1(nO, nV)
  double precision, intent(in)  :: t2(nO, nO, nV, nV)
  double precision, intent(out) :: energy

  double precision, allocatable :: W(:,:,:,:,:,:)
  double precision, allocatable :: V(:,:,:,:,:,:)
  integer :: i,j,k,a,b,c

  allocate(W(nO,nO,nO,nV,nV,nV))
  allocate(V(nO,nO,nO,nV,nV,nV))

  call form_w(nO,nV,t2,W)
  call form_v(nO,nV,t1,W,V)

  energy = 0d0
  do c = 1, nV
    do b = 1, nV
      do a = 1, nV
        do k = 1, nO
          do j = 1, nO
            do i = 1, nO
              energy = energy + (4d0 * W(i,j,k,a,b,c) + W(i,j,k,b,c,a) + W(i,j,k,c,a,b)) * (V(i,j,k,a,b,c) - V(i,j,k,c,b,a)) / (cc_space_f_o(i) + cc_space_f_o(j) + cc_space_f_o(k) - cc_space_f_v(a) - cc_space_f_v(b) - cc_space_f_v(c))  !delta_ooovvv(i,j,k,a,b,c)
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  energy = energy / 3d0

  deallocate(V,W)
end

subroutine form_w(nO,nV,t2,W)

  implicit none

  integer, intent(in)           :: nO,nV
  double precision, intent(in)  :: t2(nO, nO, nV, nV)
  double precision, intent(out) :: W(nO, nO, nO, nV, nV, nV)

  integer :: i,j,k,l,a,b,c,d

  W = 0d0
  do c = 1, nV
    print*,'W:',c,'/',nV
    do b = 1, nV
      do a = 1, nV
        do k = 1, nO
          do j = 1, nO
            do i = 1, nO

              do d = 1, nV
                W(i,j,k,a,b,c) = W(i,j,k,a,b,c) &
                ! chem (bd|ai)
                ! phys <ba|di>
                + cc_space_v_vvvo(b,a,d,i) * t2(k,j,c,d) &
                + cc_space_v_vvvo(c,a,d,i) * t2(j,k,b,d) & ! bc kj
                + cc_space_v_vvvo(a,c,d,k) * t2(j,i,b,d) & ! prev ac ik
                + cc_space_v_vvvo(b,c,d,k) * t2(i,j,a,d) & ! prev ab ij
                + cc_space_v_vvvo(c,b,d,j) * t2(i,k,a,d) & ! prev bc kj
                + cc_space_v_vvvo(a,b,d,j) * t2(k,i,c,d) ! prev ac ik
              enddo

              do l = 1, nO
                W(i,j,k,a,b,c) = W(i,j,k,a,b,c) &
                ! chem (ck|jl)
                ! phys <cj|kl>
                - cc_space_v_vooo(c,j,k,l) * t2(i,l,a,b) &
                - cc_space_v_vooo(b,k,j,l) * t2(i,l,a,c) & ! bc kj
                - cc_space_v_vooo(b,i,j,l) * t2(k,l,c,a) & ! prev ac ik
                - cc_space_v_vooo(a,j,i,l) * t2(k,l,c,b) & ! prev ab ij
                - cc_space_v_vooo(a,k,i,l) * t2(j,l,b,c) & ! prev bc kj
                - cc_space_v_vooo(c,i,k,l) * t2(j,l,b,a) ! prev ac ik
              enddo

            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

end

subroutine form_v(nO,nV,t1,w,v)

implicit none

  integer, intent(in)           :: nO,nV
  double precision, intent(in)  :: t1(nO, nV)
  double precision, intent(in)  :: W(nO, nO, nO, nV, nV, nV)
  double precision, intent(out) :: V(nO, nO, nO, nV, nV, nV)

  integer :: i,j,k,a,b,c

  V = 0d0
  do c = 1, nV
    do b = 1, nV
      do a = 1, nV
        do k = 1, nO
          do j = 1, nO
            do i = 1, nO
              V(i,j,k,a,b,c) = V(i,j,k,a,b,c) + W(i,j,k,a,b,c) &
              + cc_space_v_vvoo(b,c,j,k) * t1(i,a) &
              + cc_space_v_vvoo(a,c,i,k) * t1(j,b) &
              + cc_space_v_vvoo(a,b,i,j) * t1(k,c)
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

end

! Main

subroutine ccsd_par_t_space_v2(nO,nV,t1,t2,f_o,f_v,v_vvvo,v_vvoo,v_vooo,energy)

  implicit none

  integer, intent(in)           :: nO,nV
  double precision, intent(in)  :: t1(nO,nV), f_o(nO), f_v(nV)
  double precision, intent(in)  :: t2(nO,nO,nV,nV)
  double precision, intent(in)  :: v_vvvo(nV,nV,nV,nO), v_vvoo(nV,nV,nO,nO), v_vooo(nV,nO,nO,nO)
  double precision, intent(out) :: energy

  double precision, allocatable :: W(:,:,:,:,:,:)
  double precision, allocatable :: V(:,:,:,:,:,:)
  double precision, allocatable :: W_ijk(:,:,:), V_ijk(:,:,:)
  double precision, allocatable :: X_vvvo(:,:,:,:), X_ovoo(:,:,:,:), X_vvoo(:,:,:,:)
  double precision, allocatable :: T_vvoo(:,:,:,:), T_ovvo(:,:,:,:), T_vo(:,:)
  integer                       :: i,j,k,l,a,b,c,d
  double precision              :: e,ta,tb, delta, delta_ijk

  !allocate(W(nV,nV,nV,nO,nO,nO))
  !allocate(V(nV,nV,nV,nO,nO,nO))
  allocate(W_ijk(nV,nV,nV), V_ijk(nV,nV,nV))
  allocate(X_vvvo(nV,nV,nV,nO), X_ovoo(nO,nV,nO,nO), X_vvoo(nV,nV,nO,nO))
  allocate(T_vvoo(nV,nV,nO,nO), T_ovvo(nO,nV,nV,nO), T_vo(nV,nO))

  ! Temporary arrays
  !$OMP PARALLEL &
  !$OMP SHARED(nO,nV,T_vvoo,T_ovvo,T_vo,X_vvvo,X_ovoo,X_vvoo, &
  !$OMP t1,t2,v_vvvo,v_vooo,v_vvoo) &
  !$OMP PRIVATE(a,b,c,d,i,j,k,l) &
  !$OMP DEFAULT(NONE)

  !v_vvvo(b,a,d,i) * t2(k,j,c,d) &
  !X_vvvo(d,b,a,i) * T_vvoo(d,c,k,j)

  !$OMP DO collapse(3)
  do i = 1, nO
    do a = 1, nV
      do b = 1, nV
        do d = 1, nV
          X_vvvo(d,b,a,i) = v_vvvo(b,a,d,i)
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
          T_vvoo(d,c,k,j) = t2(k,j,c,d)
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO nowait

  !v_vooo(c,j,k,l) * t2(i,l,a,b) &
  !X_ovoo(l,c,j,k) * T_ovvo(l,a,b,i) &

  !$OMP DO collapse(3)
  do k = 1, nO
    do j = 1, nO
      do c = 1, nV
        do l = 1, nO
           X_ovoo(l,c,j,k) = v_vooo(c,j,k,l)
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
          T_ovvo(l,a,b,i) = t2(i,l,a,b)
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO nowait

  !v_vvoo(b,c,j,k) * t1(i,a) &
  !X_vvoo(b,c,k,j) * T1_vo(a,i) &

  !$OMP DO collapse(3)
  do j = 1, nO
    do k = 1, nO
      do c = 1, nV
        do b = 1, nV
          X_vvoo(b,c,k,j) = v_vvoo(b,c,j,k)
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO nowait

  !$OMP DO collapse(1)
  do i = 1, nO
    do a = 1, nV
      T_vo(a,i) = t1(i,a)
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  call wall_time(ta)
  energy = 0d0
  do i = 1, nO
    do j = 1, nO
      do k = 1, nO
        delta_ijk = f_o(i) + f_o(j) + f_o(k)
        call form_w_ijk(nO,nV,i,j,k,T_vvoo,T_ovvo,X_vvvo,X_ovoo,W_ijk)
        call form_v_ijk(nO,nV,i,j,k,T_vo,X_vvoo,W_ijk,V_ijk)
        !$OMP PARALLEL &
        !$OMP SHARED(energy,nV,i,j,k,W_ijk,V_ijk,f_o,f_v,delta_ijk) &
        !$OMP PRIVATE(a,b,c,e,delta) &
        !$OMP DEFAULT(NONE)
        e = 0d0
        !$OMP DO
        do c = 1, nV
          do b = 1, nV
            do a = 1, nV
              delta = 1d0 / (delta_ijk - f_v(a) - f_v(b) - f_v(c))
              !energy = energy + (4d0 * W(i,j,k,a,b,c) + W(i,j,k,b,c,a) + W(i,j,k,c,a,b)) * (V(i,j,k,a,b,c) - V(i,j,k,c,b,a)) / (cc_space_f_o(i) + cc_space_f_o(j) + cc_space_f_o(k) - cc_space_f_v(a) - cc_space_f_v(b) - cc_space_f_v(c))  !delta_ooovvv(i,j,k,a,b,c)
              e = e + (4d0 * W_ijk(a,b,c) + W_ijk(b,c,a) + W_ijk(c,a,b)) &
                       * (V_ijk(a,b,c) - V_ijk(c,b,a)) * delta
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP CRITICAL
        energy = energy + e
        !$OMP END CRITICAL
        !$OMP END PARALLEL
      enddo
    enddo
    call wall_time(tb)
    write(*,'(F12.2,A5,F12.2,A2)') dble(i)/dble(nO)*100d0, '% in ', tb - ta, ' s'
  enddo

  energy = energy / 3d0

  deallocate(W_ijk,V_ijk,X_vvvo,X_ovoo,T_vvoo,T_ovvo,T_vo)
  !deallocate(V,W)
end

! W_ijk

subroutine form_w_ijk(nO,nV,i,j,k,T_vvoo,T_ovvo,X_vvvo,X_ovoo,W)

  implicit none

  integer, intent(in)           :: nO,nV,i,j,k
  !double precision, intent(in) :: t2(nO,nO,nV,nV)
  double precision, intent(in)  :: T_vvoo(nV,nV,nO,nO), T_ovvo(nO,nV,nV,nO)
  double precision, intent(in)  :: X_vvvo(nV,nV,nV,nO), X_ovoo(nO,nV,nO,nO)
  double precision, intent(out) :: W(nV,nV,nV)!,nO,nO,nO)

  integer :: l,a,b,c,d
  double precision, allocatable, dimension(:,:,:) :: X, Y, Z

  !W = 0d0
  !do i = 1, nO
  !  do j = 1, nO
  !    do k = 1, nO

  allocate(X(nV,nV,nV))
  allocate(Y(nV,nV,nV))
  allocate(Z(nV,nV,nV))

  !$OMP PARALLEL DO
  do b = 1, nV
    do a = 1, nV
      do d = 1, nV
        Z(d,a,b) = X_vvvo(d,b,a,i)
      enddo
    enddo
  enddo
  !$OMP END PARALLEL DO

  call dgemm('T','N',nV*nV,nV,nV, 1.d0, &
     Z, nV, T_vvoo(1,1,k,j), nV, 0.d0, W, nV*nV)

  !$OMP PARALLEL DO
  do c = 1, nV
    do a = 1, nV
      do d = 1, nV
        Z(d,a,c) = X_vvvo(d,c,a,i)
      enddo
    enddo
  enddo
  !$OMP END PARALLEL DO

  call dgemm('T','N',nV*nV,nV,nV, 1.d0, &
     Z, nV, T_vvoo(1,1,j,k), nV, 0.d0, Y, nV*nV)

  call dgemm('T','N',nV*nV,nV,nV, 1.d0, &
     X_vvvo(1,1,1,k), nV, T_vvoo(1,1,j,i), nV, 1.d0, Y, nV*nV)

  call dgemm('T','N',nV,nV*nV,nV, 1.d0, &
     T_vvoo(1,1,i,j), nV, X_vvvo(1,1,1,k), nV, 1.d0, W, nV)

  call dgemm('T','N',nV,nV*nV,nV, 1.d0, &
     T_vvoo(1,1,i,k), nV, X_vvvo(1,1,1,j), nV, 1.d0, Y, nV)

  call dgemm('T','N',nV*nV,nV,nV, 1.d0, &
     X_vvvo(1,1,1,j), nV, T_vvoo(1,1,k,i), nV, 1.d0, W, nV*nV)

  deallocate(Z)


  allocate(Z(nO,nV,nV))

  call dgemm('T','N',nV*nV,nV,nO, -1.d0, &
     T_ovvo(1,1,1,i), nO, X_ovoo(1,1,j,k), nO, 1.d0, W, nV*nV)

  call dgemm('T','N',nV*nV,nV,nO, -1.d0, &
     T_ovvo(1,1,1,i), nO, X_ovoo(1,1,k,j), nO, 1.d0, Y, nV*nV)

  !$OMP PARALLEL DO
  do c = 1, nV
    do a = 1, nV
      do l = 1, nO
        Z(l,a,c) = T_ovvo(l,c,a,k)
      enddo
    enddo
  enddo
  !$OMP END PARALLEL DO

  call dgemm('T','N',nV*nV,nV,nO, -1.d0, &
     Z, nO, X_ovoo(1,1,i,j), nO, 1.d0, Y, nV*nV)

  call dgemm('T','N',nV,nV*nV,nO, -1.d0, &
     X_ovoo(1,1,j,i), nO, T_ovvo(1,1,1,k), nO, 1.d0, Y, nV)

  call dgemm('T','N',nV,nV*nV,nO, -1.d0, &
     X_ovoo(1,1,k,i), nO, T_ovvo(1,1,1,j), nO, 1.d0, W, nV)

  !$OMP PARALLEL DO
  do b = 1, nV
    do a = 1, nV
      do l = 1, nO
        Z(l,a,b) = T_ovvo(l,b,a,j)
      enddo
    enddo
  enddo
  !$OMP END PARALLEL DO

  call dgemm('T','N',nV*nV,nV,nO, -1.d0, &
     Z, nO, X_ovoo(1,1,i,k), nO, 1.d0, W, nV*nV)

  !$OMP PARALLEL DO
  do c = 1, nV
    do b = 1, nV
      do a = 1, nV
        W(a,b,c) = W(a,b,c) + Y(a,c,b)
      enddo
    enddo
  enddo
  !$OMP END PARALLEL DO

  deallocate(X,Y,Z)


!  !$OMP PARALLEL &
!  !$OMP SHARED(nO,nV,i,j,k,T_vvoo,T_ovvo,X_vvvo,X_ovoo,W) &
!  !$OMP PRIVATE(a,b,c,d,l) &
!  !$OMP DEFAULT(NONE)
!
!  !$OMP DO collapse(2)
!  do c = 1, nV
!    do b = 1, nV
!      do a = 1, nV
!        W(a,b,c) = 0.d0
!
!        do d = 1, nV
!          !W(i,j,k,a,b,c) = W(i,j,k,a,b,c) &
!          W(a,b,c) = W(a,b,c) &
!          ! chem (bd|ai)
!          ! phys <ba|di>
!          !+ cc_space_v_vvvo(b,a,d,i) * t2(k,j,c,d) &
!          !+ cc_space_v_vvvo(c,a,d,i) * t2(j,k,b,d) & ! bc kj
!          !+ cc_space_v_vvvo(a,c,d,k) * t2(j,i,b,d) & ! prev ac ik
!          !+ cc_space_v_vvvo(b,c,d,k) * t2(i,j,a,d) & ! prev ab ij
!          !+ cc_space_v_vvvo(c,b,d,j) * t2(i,k,a,d) & ! prev bc kj
!          !+ cc_space_v_vvvo(a,b,d,j) * t2(k,i,c,d) ! prev ac ik
!         + X_vvvo(d,b,a,i) * T_vvoo(d,c,k,j) &
!         + X_vvvo(d,c,a,i) * T_vvoo(d,b,j,k) & ! bc kj
!         + X_vvvo(d,a,c,k) * T_vvoo(d,b,j,i) & ! prev ac ik
!          + X_vvvo(d,b,c,k) * T_vvoo(d,a,i,j) & ! prev ab ij
!          + X_vvvo(d,c,b,j) * T_vvoo(d,a,i,k) & ! prev bc kj
!          + X_vvvo(d,a,b,j) * T_vvoo(d,c,k,i) ! prev ac ik
!        enddo
!
!      enddo
!    enddo
!  enddo
!  !$OMP END DO nowait
!
!  !$OMP DO collapse(2)
!  do c = 1, nV
!    do b = 1, nV
!      do a = 1, nV
!
!        do l = 1, nO
!          !W(i,j,k,a,b,c) = W(i,j,k,a,b,c) &
!          W(a,b,c) = W(a,b,c) &
!          ! chem (ck|jl)
!          ! phys <cj|kl>
!          !- cc_space_v_vooo(c,j,k,l) * t2(i,l,a,b) &
!          !- cc_space_v_vooo(b,k,j,l) * t2(i,l,a,c) & ! bc kj
!          !- cc_space_v_vooo(b,i,j,l) * t2(k,l,c,a) & ! prev ac ik
!          !- cc_space_v_vooo(a,j,i,l) * t2(k,l,c,b) & ! prev ab ij
!          !- cc_space_v_vooo(a,k,i,l) * t2(j,l,b,c) & ! prev bc kj
!          !- cc_space_v_vooo(c,i,k,l) * t2(j,l,b,a) ! prev ac ik
!          - T_ovvo(l,a,b,i) * X_ovoo(l,c,j,k) &
!          - T_ovvo(l,a,c,i) * X_ovoo(l,b,k,j) & ! bc kj
!          - T_ovvo(l,c,a,k) * X_ovoo(l,b,i,j) & ! prev ac ik
!          - T_ovvo(l,c,b,k) * X_ovoo(l,a,j,i) & ! prev ab ij
!          - T_ovvo(l,b,c,j) * X_ovoo(l,a,k,i) & ! prev bc kj
!          - T_ovvo(l,b,a,j) * X_ovoo(l,c,i,k) ! prev ac ik
!        enddo
!
!      enddo
!    enddo
!  enddo
!  !$OMP END DO
!  !$OMP END PARALLEL

  !    enddo
  !  enddo
  !enddo

end

! V_ijk

subroutine form_v_ijk(nO,nV,i,j,k,T_vo,X_vvoo,w,v)

implicit none

  integer, intent(in)           :: nO,nV,i,j,k
  !double precision, intent(in)  :: t1(nO,nV)
  double precision, intent(in)  :: T_vo(nV,nO)
  double precision, intent(in)  :: X_vvoo(nV,nV,nO,nO)
  double precision, intent(in)  :: W(nV,nV,nV)!,nO,nO,nO)
  double precision, intent(out) :: V(nV,nV,nV)!,nO,nO,nO)

  integer :: a,b,c

  !V = 0d0
  !do i = 1, nO
  !  do j = 1, nO
  !    do k = 1, nO

  !$OMP PARALLEL &
  !$OMP SHARED(nO,nV,i,j,k,T_vo,X_vvoo,W,V) &
  !$OMP PRIVATE(a,b,c) &
  !$OMP DEFAULT(NONE)
  !$OMP DO collapse(2)
  do c = 1, nV
    do b = 1, nV
      do a = 1, nV
        !V(i,j,k,a,b,c) = V(i,j,k,a,b,c) + W(i,j,k,a,b,c) &
        V(a,b,c) = W(a,b,c) &
        !+ cc_space_v_vvoo(b,c,j,k) * t1(i,a) &
        !+ cc_space_v_vvoo(a,c,i,k) * t1(j,b) &
        !+ cc_space_v_vvoo(a,b,i,j) * t1(k,c)
        + X_vvoo(b,c,k,j) * T_vo(a,i) &
        + X_vvoo(a,c,k,i) * T_vo(b,j) &
        + X_vvoo(a,b,j,i) * T_vo(c,k)
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  !    enddo
  !  enddo
  !enddo

end

