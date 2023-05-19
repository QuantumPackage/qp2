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
  double precision, allocatable :: W_bca(:,:,:), V_cba(:,:,:)
  double precision, allocatable :: X_vvvo(:,:,:,:), X_ovoo(:,:,:,:), X_vvoo(:,:,:,:)
  double precision, allocatable :: T_vvoo(:,:,:,:), T_ovvo(:,:,:,:), T_vo(:,:)
  integer                       :: i,j,k,l,a,b,c,d
  double precision              :: e,ta,tb, delta, delta_abc

  !allocate(W(nV,nV,nV,nO,nO,nO))
  !allocate(V(nV,nV,nV,nO,nO,nO))
  allocate(W_abc(nO,nO,nO), V_abc(nO,nO,nO), W_cab(nO,nO,nO))
  allocate(W_bca(nO,nO,nO), V_cba(nO,nO,nO), W_cba(nO,nO,nO))
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
  do c = 1, nV
    do b = 1, nV
      do a = 1, nV
        delta_abc = f_v(a) + f_v(b) + f_v(c)
        call form_w_abc(nO,nV,a,b,c,T_vvoo,T_ovvo,X_vvvo,X_ovoo,W_abc)
        call form_w_abc(nO,nV,b,c,a,T_vvoo,T_ovvo,X_vvvo,X_ovoo,W_bca)
        call form_w_abc(nO,nV,c,a,b,T_vvoo,T_ovvo,X_vvvo,X_ovoo,W_cab)
        call form_w_abc(nO,nV,c,b,a,T_vvoo,T_ovvo,X_vvvo,X_ovoo,W_cba)

        call form_v_abc(nO,nV,a,b,c,T_vo,X_vvoo,W_abc,V_abc)
        call form_v_abc(nO,nV,c,b,a,T_vo,X_vvoo,W_cba,V_cba)
        !$OMP PARALLEL                                               &
            !$OMP SHARED(energy,nO,a,b,c,W_abc,W_cab,W_bca,V_abc,V_cba,f_o,f_v,delta_abc)&
            !$OMP PRIVATE(i,j,k,e,delta)                             &
            !$OMP DEFAULT(NONE)
        e = 0d0
        !$OMP DO
        do i = 1, nO
          do j = 1, nO
            do k = 1, nO
              delta = 1d0 / (f_o(i) + f_o(j) + f_o(k) - delta_abc)
              !energy = energy + (4d0 * W(i,j,k,a,b,c) + W(i,j,k,b,c,a) + W(i,j,k,c,a,b)) * (V(i,j,k,a,b,c) - V(i,j,k,c,b,a)) / (cc_space_f_o(i) + cc_space_f_o(j) + cc_space_f_o(k) - cc_space_f_v(a) - cc_space_f_v(b) - cc_space_f_v(c))  !delta_ooovvv(i,j,k,a,b,c)
              e = e + (4d0 * W_abc(i,j,k) + W_bca(i,j,k) + W_cab(i,j,k))&
                  * (V_abc(i,j,k) - V_cba(i,j,k)) * delta
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

  deallocate(W_abc,V_abc,W_cab,V_cba,W_bca,X_vvvo,X_ovoo,T_vvoo,T_ovvo,T_vo)
  !deallocate(V,W)
end


subroutine form_w_abc(nO,nV,a,b,c,T_vvoo,T_ovvo,X_vvvo,X_ovoo,W_abc)

  implicit none

  integer, intent(in)           :: nO,nV,a,b,c
  !double precision, intent(in) :: t2(nO,nO,nV,nV)
  double precision, intent(in)  :: T_vvoo(nV,nV,nO,nO), T_ovvo(nO,nV,nV,nO)
  double precision, intent(in)  :: X_vvvo(nV,nV,nV,nO), X_ovoo(nO,nV,nO,nO)
  double precision, intent(out) :: W_abc(nO,nO,nO)

  integer :: l,i,j,k,d


  !$OMP PARALLEL &
  !$OMP SHARED(nO,nV,a,b,c,T_vvoo,T_ovvo,X_vvvo,X_ovoo,W_abc) &
  !$OMP PRIVATE(i,j,k,d,l) &
  !$OMP DEFAULT(NONE)

  !$OMP DO collapse(3)
  do k = 1, nO
    do j = 1, nO
      do i = 1, nO
        W_abc(i,j,k) = 0.d0

        do d = 1, nV
          W_abc(i,j,k) = W_abc(i,j,k) &
                 + X_vvvo(d,b,a,i) * T_vvoo(d,c,k,j) &
                 + X_vvvo(d,c,a,i) * T_vvoo(d,b,j,k) &
                 + X_vvvo(d,a,c,k) * T_vvoo(d,b,j,i) &
                 + X_vvvo(d,b,c,k) * T_vvoo(d,a,i,j) &
                 + X_vvvo(d,c,b,j) * T_vvoo(d,a,i,k) &
                 + X_vvvo(d,a,b,j) * T_vvoo(d,c,k,i)

        enddo

        do l = 1, nO
          W_abc(i,j,k) = W_abc(i,j,k) &
              - T_ovvo(l,a,b,i) * X_ovoo(l,c,j,k) &
              - T_ovvo(l,a,c,i) * X_ovoo(l,b,k,j) & ! bc kj
              - T_ovvo(l,c,a,k) * X_ovoo(l,b,i,j) & ! prev ac ik
              - T_ovvo(l,c,b,k) * X_ovoo(l,a,j,i) & ! prev ab ij
              - T_ovvo(l,b,c,j) * X_ovoo(l,a,k,i) & ! prev bc kj
              - T_ovvo(l,b,a,j) * X_ovoo(l,c,i,k) ! prev ac ik
        enddo

      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL


end


! V_abc

subroutine form_v_abc(nO,nV,a,b,c,T_vo,X_vvoo,W,V)

implicit none

  integer, intent(in)           :: nO,nV,a,b,c
  !double precision, intent(in)  :: t1(nO,nV)
  double precision, intent(in)  :: T_vo(nV,nO)
  double precision, intent(in)  :: X_vvoo(nV,nV,nO,nO)
  double precision, intent(in)  :: W(nO,nO,nO)
  double precision, intent(out) :: V(nO,nO,nO)

  integer :: i,j,k

  !$OMP PARALLEL &
  !$OMP SHARED(nO,nV,a,b,c,T_vo,X_vvoo,W,V) &
  !$OMP PRIVATE(i,j,k) &
  !$OMP DEFAULT(NONE)
  !$OMP DO collapse(2)
  do k = 1, nO
    do j = 1, nO
      do i = 1, nO
        !V(i,j,k,a,b,c) = V(i,j,k,a,b,c) + W(i,j,k,a,b,c) &
        V(i,j,k) = W(i,j,k) &
           + X_vvoo(b,c,k,j) * T_vo(a,i) &
           + X_vvoo(a,c,k,i) * T_vo(b,j) &
           + X_vvoo(a,b,j,i) * T_vo(c,k)
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

end

