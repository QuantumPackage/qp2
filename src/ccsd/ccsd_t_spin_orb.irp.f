! v1

subroutine ccsd_par_t_spin(nO,nV,t1,t2,f_o,f_v,f_ov,v_ooov,v_vvoo,v_vvvo,energy)

  implicit none

  integer, intent(in)           :: nO, nV
  double precision, intent(in)  :: t1(nO,nV), t2(nO,nO,nV,nV)
  double precision, intent(in)  :: f_o(nO), f_v(nV), f_ov(nO,nV)
  double precision, intent(in)  :: v_ooov(nO,nO,nO,nV)
  double precision, intent(in)  :: v_vvoo(nV,nV,nO,nO), v_vvvo(nV,nV,nV,nO)
  double precision, intent(out) :: energy

  double precision, allocatable :: t3(:,:,:,:,:,:), s(:,:)
  double precision              :: e_t, e_st, e_dt, delta_abc, delta
  integer                       :: i,j,k,l,m,a,b,c,d,e

  allocate(t3(nO,nO,nO,nV,nV,nV), s(nO,nV))

  t3 = 0d0
  
  ! T3
  do c = 1, nV
    do b = 1, nV
      do a = 1, nV
        delta_abc = f_v(a) + f_v(b) + f_v(c)
        do k = 1, nO
          do j = 1, nO
            do i = 1, nO
              delta = f_o(i) + f_o(j) + f_o(k) - delta_abc
              do e = 1, nV
                t3(i,j,k,a,b,c) = t3(i,j,k,a,b,c) &
                  + t2(j,k,a,e) * v_vvvo(b,c,e,i) &
                  - t2(i,k,a,e) * v_vvvo(b,c,e,j) & ! - P(ij)
                  - t2(j,i,a,e) * v_vvvo(b,c,e,k) & ! - P(ik)
                  - t2(j,k,b,e) * v_vvvo(a,c,e,i) & ! - P(ab)
                  - t2(j,k,c,e) * v_vvvo(b,a,e,i) & ! - P(ac)
                  + t2(i,k,b,e) * v_vvvo(a,c,e,j) & ! + P(ij) P(ab)
                  + t2(i,k,c,e) * v_vvvo(b,a,e,j) & ! + P(ij) P(ac)
                  + t2(j,i,b,e) * v_vvvo(a,c,e,k) & ! + P(ik) P(ab)
                  + t2(j,i,c,e) * v_vvvo(b,a,e,k)   ! + P(ik) P(ac)
              enddo
              do m = 1, nO
                t3(i,j,k,a,b,c) = t3(i,j,k,a,b,c) &
                  + t2(m,i,b,c) * v_ooov(j,k,m,a) &
                  - t2(m,j,b,c) * v_ooov(i,k,m,a) & ! - P(ij)
                  - t2(m,k,b,c) * v_ooov(j,i,m,a) & ! - P(ik)
                  - t2(m,i,a,c) * v_ooov(j,k,m,b) & ! - P(ab)
                  - t2(m,i,b,a) * v_ooov(j,k,m,c) & ! - P(ac)
                  + t2(m,j,a,c) * v_ooov(i,k,m,b) & ! + P(ij) P(ab)
                  + t2(m,j,b,a) * v_ooov(i,k,m,c) & ! + P(ij) P(ac)
                  + t2(m,k,a,c) * v_ooov(j,i,m,b) & ! + P(ik) P(ab)
                  + t2(m,k,b,a) * v_ooov(j,i,m,c)   ! + P(ik) P(ac)
              enddo
                t3(i,j,k,a,b,c) = t3(i,j,k,a,b,c) * (1d0 / delta)
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
  

  ! E_T
  e_t = 0d0
  do c = 1, nV
    do b = 1, nV
      do a = 1, nV
        delta_abc = f_v(a) + f_v(b) + f_v(c)
        do k = 1, nO
          do j = 1, nO
            do i = 1, nO
              delta = f_o(i) + f_o(j) + f_o(k) - delta_abc
              e_t = e_t + t3(i,j,k,a,b,c) * delta * t3(i,j,k,a,b,c)
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
  e_t = e_t / 36d0

  ! E_ST
  s = 0d0
  do c = 1, nV
    do b = 1, nV
      do a = 1, nV
        do k = 1, nO
          do j = 1, nO
            do i = 1, nO
              s(i,a) = s(i,a) + v_vvoo(b,c,j,k) * t3(i,j,k,a,b,c)
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  e_st = 0d0
  do a = 1, nV
    do i = 1, nO
      e_st = e_st + s(i,a) * t1(i,a)
    enddo
  enddo
  e_st = e_st * 0.25d0

  ! E_DT
  e_dt = 0d0
  do c = 1, nV
    do b = 1, nV
      do a = 1, nV
        do k = 1, nO
          do j = 1, nO
            do i = 1, nO
              e_dt = e_dt + t2(i,j,a,b) * f_ov(k,c) * t3(i,j,k,a,b,c)
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
  e_dt = e_dt * 0.25d0

  ! (T)
  !print*,e_t,e_st,e_dt
  energy = e_t + e_st + e_dt

  deallocate(t3,s)
  
end

! v2

subroutine ccsd_par_t_spin_v2(nO,nV,t1,t2,f_o,f_v,f_ov,v_ooov,v_vvoo,energy)

  implicit none

  integer, intent(in)           :: nO, nV
  double precision, intent(in)  :: t1(nO,nV), t2(nO,nO,nV,nV)
  double precision, intent(in)  :: f_o(nO), f_v(nV), f_ov(nO,nV)
  double precision, intent(in)  :: v_ooov(nO,nO,nO,nV)
  double precision, intent(in)  :: v_vvoo(nV,nV,nO,nO)
  double precision, intent(out) :: energy

  double precision, allocatable :: t3_bc(:,:,:,:), s(:,:), e_t(:), e_dt(:)
  double precision, allocatable :: A_vovv(:,:,:,:), v_vvvo(:,:,:,:)
  double precision, allocatable :: T_voov(:,:,:,:), B_ooov(:,:,:,:)
  double precision              :: e_st, delta_abc, delta, ta, tb
  integer                       :: i,j,k,l,m,a,b,c,d,e

  allocate(t3_bc(nO,nO,nO,nV), s(nO,nV), e_t(nV), e_dt(nV))
  allocate(A_vovv(nV,nO,nV,nV),v_vvvo(nV,nV,nV,nO),T_voov(nV,nO,nO,nV),B_ooov(nO,nO,nO,nV))

  call gen_v_spin(cc_nV_m,cc_nV_m,cc_nV_m,cc_nO_m, &
       cc_nV_S,cc_nV_S,cc_nV_S,cc_nO_S, &
       cc_list_vir_spin,cc_list_vir_spin,cc_list_vir_spin,cc_list_occ_spin, &
       nV,nV,nV,nO, v_vvvo)
  
  ! Init 
  s = 0d0
  e_t = 0d0
  e_st = 0d0
  e_dt = 0d0

  call wall_time(ta)
  !$OMP PARALLEL &
  !$OMP PRIVATE(i,j,k,m,a,b,c,e) &
  !$OMP SHARED(A_vovv,ta,tb,t3_bc,s,e_t,e_st,e_dt,t2,v_vvvo,v_ooov, &
  !$OMP v_vvoo,f_o,f_v,f_ov,delta,delta_abc,nO,nV,T_voov,B_ooov) &
  !$OMP DEFAULT(NONE)

  !$OMP DO collapse(3)
  do c = 1, nV
    do b = 1, nV
      do i = 1, nO
        do e = 1, nV
          A_vovv(e,i,b,c) = v_vvvo(b,c,e,i)
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO nowait

  !$omp do collapse(3)
  do a = 1, nV
    do k = 1, nO
      do j = 1, nO
        do e = 1, nV
          T_voov(e,j,k,a) = t2(j,k,a,e)
        enddo
      enddo
    enddo
  enddo
  !$omp end do nowait

  !$omp do collapse(3)
  do a = 1, nV
    do k = 1, nO
      do j = 1, nO
        do m = 1, nO
          B_ooov(m,j,k,a) = v_ooov(j,k,m,a)
        enddo
      enddo
    enddo
  enddo
  !$omp end do

  do c = 1, nV
    do b = 1, nV

      ! T3(:,:,:,:,b,c)
      ! Init
      !$OMP DO collapse(3)
      do a = 1, nV
        do k = 1, nO
          do j = 1, nO
            do i = 1, nO
              t3_bc(i,j,k,a) = 0d0
            enddo
          enddo
        enddo
      enddo
      !$OMP END DO
      
      !$OMP DO collapse(3)
      do a = 1, nV
        do k = 1, nO
          do j = 1, nO
            do i = 1, nO
              do e = 1, nV
                t3_bc(i,j,k,a) = t3_bc(i,j,k,a) &
                   !+ t2(j,k,a,e) * v_vvvo(b,c,e,i) &
                   !- t2(i,k,a,e) * v_vvvo(b,c,e,j) & ! - P(ij)
                   !- t2(j,i,a,e) * v_vvvo(b,c,e,k) & ! - P(ik)
                   !- t2(j,k,b,e) * v_vvvo(a,c,e,i) & ! - P(ab)
                   !- t2(j,k,c,e) * v_vvvo(b,a,e,i) & ! - P(ac)
                   !+ t2(i,k,b,e) * v_vvvo(a,c,e,j) & ! + P(ij) P(ab)
                   !+ t2(i,k,c,e) * v_vvvo(b,a,e,j) & ! + P(ij) P(ac)
                   !+ t2(j,i,b,e) * v_vvvo(a,c,e,k) & ! + P(ik) P(ab)
                   !+ t2(j,i,c,e) * v_vvvo(b,a,e,k)   ! + P(ik) P(ac)
                   + T_voov(e,j,k,a) * A_vovv(e,i,b,c) &
                   - T_voov(e,i,k,a) * A_vovv(e,j,b,c) & ! - P(ij)
                   - T_voov(e,j,i,a) * A_vovv(e,k,b,c) & ! - P(ik)
                   - T_voov(e,j,k,b) * A_vovv(e,i,a,c) & ! - P(ab)
                   - T_voov(e,j,k,c) * A_vovv(e,i,b,a) & ! - P(ac)
                   + T_voov(e,i,k,b) * A_vovv(e,j,a,c) & ! + P(ij) P(ab)
                   + T_voov(e,i,k,c) * A_vovv(e,j,b,a) & ! + P(ij) P(ac)
                   + T_voov(e,j,i,b) * A_vovv(e,k,a,c) & ! + P(ik) P(ab)
                   + T_voov(e,j,i,c) * A_vovv(e,k,b,a)   ! + P(ik) P(ac)
              enddo
            enddo
          enddo
        enddo
      enddo
      !$OMP END DO
      
      !$OMP DO collapse(3)
      do a = 1, nV
        do k = 1, nO
          do j = 1, nO
            do i = 1, nO
              do m = 1, nO
                t3_bc(i,j,k,a) = t3_bc(i,j,k,a) &
                   !+ t2(m,i,b,c) * v_ooov(j,k,m,a) &
                   !- t2(m,j,b,c) * v_ooov(i,k,m,a) & ! - P(ij)
                   !- t2(m,k,b,c) * v_ooov(j,i,m,a) & ! - P(ik)
                   !- t2(m,i,a,c) * v_ooov(j,k,m,b) & ! - P(ab)
                   !- t2(m,i,b,a) * v_ooov(j,k,m,c) & ! - P(ac)
                   !+ t2(m,j,a,c) * v_ooov(i,k,m,b) & ! + P(ij) P(ab)
                   !+ t2(m,j,b,a) * v_ooov(i,k,m,c) & ! + P(ij) P(ac)
                   !+ t2(m,k,a,c) * v_ooov(j,i,m,b) & ! + P(ik) P(ab)
                   !+ t2(m,k,b,a) * v_ooov(j,i,m,c)   ! + P(ik) P(ac)
                   + t2(m,i,b,c) * B_ooov(m,j,k,a) &
                   - t2(m,j,b,c) * B_ooov(m,i,k,a) & ! - P(ij)
                   - t2(m,k,b,c) * B_ooov(m,j,i,a) & ! - P(ik)
                   - t2(m,i,a,c) * B_ooov(m,j,k,b) & ! - P(ab)
                   - t2(m,i,b,a) * B_ooov(m,j,k,c) & ! - P(ac)
                   + t2(m,j,a,c) * B_ooov(m,i,k,b) & ! + P(ij) P(ab)
                   + t2(m,j,b,a) * B_ooov(m,i,k,c) & ! + P(ij) P(ac)
                   + t2(m,k,a,c) * B_ooov(m,j,i,b) & ! + P(ik) P(ab)
                   + t2(m,k,b,a) * B_ooov(m,j,i,c)   ! + P(ik) P(ac)
              enddo
            enddo
          enddo
        enddo
      enddo
      !$OMP END DO

      !$OMP DO
      do a = 1, nV
        delta_abc = f_v(a) + f_v(b) + f_v(c)
        do k = 1, nO
          do j = 1, nO
            do i = 1, nO
               delta = f_o(i) + f_o(j) + f_o(k) - delta_abc 
               t3_bc(i,j,k,a) = t3_bc(i,j,k,a) * (1d0 / delta)
             enddo
           enddo
         enddo
       enddo
      !$OMP END DO

      ! E_T
      !$OMP DO
      do a = 1, nV
        delta_abc = f_v(a) + f_v(b) + f_v(c)
        do k = 1, nO
          do j = 1, nO
            do i = 1, nO
              delta = f_o(i) + f_o(j) + f_o(k) - delta_abc
              e_t(a) = e_t(a) + t3_bc(i,j,k,a) * delta * t3_bc(i,j,k,a)
            enddo
          enddo
        enddo
      enddo
      !$OMP END DO nowait

      ! E_ST
      !$OMP DO
      do a = 1, nV
        do k = 1, nO
          do j = 1, nO
            do i = 1, nO
              s(i,a) = s(i,a) + v_vvoo(b,c,j,k) * t3_bc(i,j,k,a)
            enddo
          enddo
        enddo
      enddo
      !$OMP END DO nowait

      ! E_DT
      !$OMP DO
      do a = 1, nV
        do k = 1, nO
          do j = 1, nO
            do i = 1, nO
              e_dt(a) = e_dt(a) + t2(i,j,a,b) * f_ov(k,c) * t3_bc(i,j,k,a)
            enddo
          enddo
        enddo
      enddo
      !$OMP END DO
    enddo
    !$OMP MASTER
    call wall_time(tb)
    write(*,'(A1,F6.2,A5,F10.2,A2)') ' ', dble(c)/dble(nV)*100d0, '% in ', tb-ta, ' s'
    !$OMP END MASTER
  enddo
  !$OMP END PARALLEL

  do a = 2, nV
    e_t(1) = e_t(1) + e_t(a)
  enddo
  
  do a = 2, nV
    e_dt(1) = e_dt(1) + e_dt(a)
  enddo

  e_t = e_t / 36d0
  
  do a = 1, nV
    do i = 1, nO
      e_st = e_st + s(i,a) * t1(i,a)
    enddo
  enddo
  e_st = e_st * 0.25d0

  e_dt = e_dt * 0.25d0

  ! (T)
  !print*,e_t(1),e_st,e_dt(1)
  energy = e_t(1) + e_st + e_dt(1)

  deallocate(t3_bc,s)
  
end
