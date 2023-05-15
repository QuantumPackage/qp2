! Main

subroutine ccsd_par_t_space_stoch(nO,nV,t1,t2,f_o,f_v,v_vvvo,v_vvoo,v_vooo,energy)

  implicit none

  integer, intent(in)           :: nO,nV
  double precision, intent(in)  :: t1(nO,nV), f_o(nO), f_v(nV)
  double precision, intent(in)  :: t2(nO,nO,nV,nV)
  double precision, intent(in)  :: v_vvvo(nV,nV,nV,nO), v_vvoo(nV,nV,nO,nO), v_vooo(nV,nO,nO,nO)
  double precision, intent(out) :: energy

  double precision, allocatable :: W(:,:,:,:,:,:)
  double precision, allocatable :: V(:,:,:,:,:,:)
  double precision, allocatable :: W_abc(:,:,:), W_cab(:,:,:), W_bca(:,:,:)
  double precision, allocatable :: W_bac(:,:,:), W_cba(:,:,:), W_acb(:,:,:)
  double precision, allocatable :: V_abc(:,:,:), V_cab(:,:,:), V_bca(:,:,:)
  double precision, allocatable :: V_bac(:,:,:), V_cba(:,:,:), V_acb(:,:,:)
  double precision, allocatable :: X_vovv(:,:,:,:), X_ooov(:,:,:,:), X_oovv(:,:,:,:)
  double precision, allocatable :: T_voov(:,:,:,:), T_oovv(:,:,:,:)
  integer                       :: i,j,k,l,a,b,c,d
  double precision              :: e,ta,tb

  call set_multiple_levels_omp(.False.)

  allocate(X_vovv(nV,nO,nV,nV), X_ooov(nO,nO,nO,nV), X_oovv(nO,nO,nV,nV))
  allocate(T_voov(nV,nO,nO,nV),T_oovv(nO,nO,nV,nV))

  !$OMP PARALLEL &
  !$OMP SHARED(nO,nV,T_voov,T_oovv,X_vovv,X_ooov,X_oovv, &
  !$OMP t1,t2,v_vvvo,v_vooo,v_vvoo) &
  !$OMP PRIVATE(a,b,c,d,i,j,k,l) &
  !$OMP DEFAULT(NONE)

  !v_vvvo(b,a,d,i) * t2(k,j,c,d) &
  !X_vovv(d,i,b,a,i) * T_voov(d,j,c,k)

  !$OMP DO
  do a = 1, nV
    do b = 1, nV
      do i = 1, nO
        do d = 1, nV
          X_vovv(d,i,b,a) = v_vvvo(b,a,d,i)
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO nowait

  !$OMP DO
  do c = 1, nV
    do j = 1, nO
      do k = 1, nO
        do d = 1, nV
          T_voov(d,k,j,c) = t2(k,j,c,d)
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO nowait

  !v_vooo(c,j,k,l) * t2(i,l,a,b) &
  !X_ooov(l,j,k,c) * T_oovv(l,i,a,b) &

  !$OMP DO
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

  !$OMP DO
  do b = 1, nV
    do a = 1, nV
      do i = 1, nO
        do l = 1, nO
          T_oovv(l,i,a,b) = t2(i,l,a,b)
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO nowait

  !X_oovv(j,k,b,c) * T1_vo(a,i) &

  !$OMP DO
  do c = 1, nV
    do b = 1, nV
      do k = 1, nO
        do j = 1, nO
          X_oovv(j,k,b,c) = v_vvoo(b,c,j,k)
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO nowait

  !$OMP END PARALLEL

  double precision, external :: ccsd_t_task_aba
  double precision, external :: ccsd_t_task_abc

  double precision, allocatable :: memo(:), Pabc(:), waccu(:)
  logical         , allocatable :: computed(:)
  integer*2       , allocatable :: abc(:,:)
  integer*8                     :: Nabc, i8
  integer*8, allocatable :: iorder(:)
  double precision :: eocc
  double precision :: Pabc_norm, sum_w


  ! Prepare table of triplets (a,b,c)

  Nabc = (int(nV,8) * int(nV+1,8) * int(nV+2,8))/6_8 - nV
  allocate (memo(Nabc), computed(Nabc), Pabc(Nabc), waccu(0:Nabc))
  allocate (abc(4,Nabc), iorder(Nabc))

!  eocc = 3.d0/dble(nO) * sum(f_o(1:nO))
  memo(:) = 0.d0
  computed(:) = .False.
  Nabc = 0_8
  do a = 1, nV
    do b = a+1, nV
      do c = b+1, nV
        Nabc = Nabc + 1_8
!        Pabc(Nabc) = 1.d0/((f_v(a) + f_v(b) + f_v(c))*(f_v(a)*f_v(b)*f_v(c))**(1.d0/2.d0))
!        Pabc(Nabc) = 1.d0/((f_v(a) + f_v(b) + f_v(c))**2)
        Pabc(Nabc) = 1.d0/(f_v(a) + f_v(b) + f_v(c))
        abc(1,Nabc) = a
        abc(2,Nabc) = b
        abc(3,Nabc) = c
      enddo

      Nabc = Nabc + 1_8
      abc(1,Nabc) = a
      abc(2,Nabc) = b
      abc(3,Nabc) = a
!      Pabc(Nabc) = 1.d0/((f_v(a) + f_v(b) + f_v(a))*(f_v(a)*f_v(b)*f_v(a))**(1.d0/2.d0))
!      Pabc(Nabc) = 1.d0/((f_v(a) + f_v(b) + f_v(a))**2)
      Pabc(Nabc) = 1.d0/(2.d0*f_v(a) + f_v(b))

      Nabc = Nabc + 1_8
      abc(1,Nabc) = b
      abc(2,Nabc) = a
      abc(3,Nabc) = b
!      Pabc(Nabc) = 1.d0/((f_v(a) + f_v(b) + f_v(b))*(f_v(b)*f_v(a)*f_v(b))**(1.d0/2.d0))
!      Pabc(Nabc) = 1.d0/((f_v(a) + f_v(b) + f_v(b))**2)
      Pabc(Nabc) = 1.d0/(f_v(a) + 2.d0*f_v(b))
    enddo
  enddo

  do i8=1,Nabc
   iorder(i8) = i8
  enddo

  ! Sort triplets in decreasing Pabc
  call dsort_big(Pabc, iorder, Nabc)

  ! Normalize
  Pabc_norm = 0.d0
  do i8=Nabc,1,-1
    Pabc_norm = Pabc_norm + Pabc(i8)
  enddo
  Pabc_norm = 1.d0/Pabc_norm
  do i8=Nabc,1,-1
    Pabc(i8) = Pabc(i8) * Pabc_norm
  enddo

  call i8set_order_big(abc, iorder, Nabc)


  ! Cumulative distribution for sampling
  waccu(Nabc) = 0.d0
  sum_w = 0.d0
  do i8=Nabc-1,1,-1
   waccu(i8) = waccu(i8+1) - Pabc(i8)
  enddo
  waccu(:) = waccu(:) + 1.d0
  waccu(0) = 0.d0

  Pabc(:) = 1.d0/Pabc(:) * (1.d0/3.d0)

  logical :: converged
  double precision :: ET, ET2, eta, variance, average, error, sample
  integer*8 :: isample, ieta, Ncomputed
  integer*8, external :: find_sample

  allocate( W_abc(nO,nO,nO), W_cab(nO,nO,nO), W_bca(nO,nO,nO), &
            W_bac(nO,nO,nO), W_cba(nO,nO,nO), W_acb(nO,nO,nO), &
            V_abc(nO,nO,nO), V_cab(nO,nO,nO), V_bca(nO,nO,nO), &
            V_bac(nO,nO,nO), V_cba(nO,nO,nO), V_acb(nO,nO,nO) )

  converged = .False.
  ET = 0.d0
  ET2 = 0.d0
  Ncomputed = 0_8
  isample = 0_8

  average = 0.d0
  variance = 0.d0
  double precision :: t00, t01
  call wall_time(t00)
!  do ieta=1,Nabc
  do while (.not.converged)
    call random_number(eta)
!    eta = eta/dble(1000)
!    do k=0,1000-1
!    ieta = find_sample(eta+dble(k)/dble(1000),waccu,Nabc)
    ieta = find_sample(eta,waccu,Nabc)
    isample = isample+1_8

    if (.not.computed(ieta)) then
      a = abc(1,ieta)
      b = abc(2,ieta)
      c = abc(3,ieta)
      if (a/=c) then
         memo(ieta) = ccsd_t_task_abc(a,b,c,nO,nV,t1,T_oovv,T_voov,V_abc, &
                         V_acb,V_bac,V_bca,V_cab,V_cba,W_abc,W_acb,W_bac, &
                         W_bca,W_cab,W_cba,X_ooov,X_oovv,X_vovv,f_o,f_v)
      else
         memo(ieta) =  ccsd_t_task_aba(a,b,nO,nV,t1,T_oovv,T_voov,V_abc, &
                       V_acb,V_bac,V_bca,V_cab,V_cba,W_abc,W_acb,W_bac, &
                       W_bca,W_cab,W_cba,X_ooov,X_oovv,X_vovv,f_o,f_v)
      endif
      computed(ieta) = .True.
      Ncomputed += 1_8
      call wall_time(t01)
      if (t01-t00 > 1.d0) then
        t00 = t01
        print *, average, dsqrt(variance/dble(isample)), real(Ncomputed)/real(Nabc), real(isample)/real(Nabc)
      endif
!       print *, memo(ieta), Pabc(ieta), memo(ieta) * Pabc(ieta)
    endif
    sample = memo(ieta) * Pabc(ieta)
    ET = ET + sample
    ET2 = ET2 + sample*sample
    average  = ET/dble(isample)
    variance = ET2/dble(isample) - average*average
    converged = (Ncomputed >= (Nabc*90_8)/100_8) .or. (isample>=1000*Nabc)
!    enddo
  enddo
        print *, average, dsqrt(variance/dble(isample)), real(Ncomputed)/real(Nabc), real(isample)/real(Nabc)
  energy = average

!  !$OMP PARALLEL                                                     &
!      !$OMP PRIVATE(a,b,c,e)                                         &
!      !$OMP PRIVATE(W_abc, W_cab, W_bca, W_bac, W_cba, W_acb,        &
!      !$OMP         V_abc, V_cab, V_bca, V_bac, V_cba, V_acb )       &
!      !$OMP DEFAULT(SHARED)
!  allocate( W_abc(nO,nO,nO), W_cab(nO,nO,nO), W_bca(nO,nO,nO), &
!            W_bac(nO,nO,nO), W_cba(nO,nO,nO), W_acb(nO,nO,nO), &
!            V_abc(nO,nO,nO), V_cab(nO,nO,nO), V_bca(nO,nO,nO), &
!            V_bac(nO,nO,nO), V_cba(nO,nO,nO), V_acb(nO,nO,nO) )
!  e = 0d0
!  !$OMP DO SCHEDULE(dynamic)
!  do a = 1, nV
!    do b = a+1, nV
!      do c = b+1, nV
!        e = e + ccsd_t_task_abc(a,b,c,nO,nV,t1,T_oovv,T_voov,V_abc, &
!                        V_acb,V_bac,V_bca,V_cab,V_cba,W_abc,W_acb,W_bac, &
!                        W_bca,W_cab,W_cba,X_ooov,X_oovv,X_vovv,f_o,f_v)
!      enddo
!    enddo
!
!    do b = 1, nV
!      if (b == a) cycle
!      e = e + ccsd_t_task_aba(a,b,nO,nV,t1,T_oovv,T_voov,V_abc, &
!                      V_acb,V_bac,V_bca,V_cab,V_cba,W_abc,W_acb,W_bac, &
!                      W_bca,W_cab,W_cba,X_ooov,X_oovv,X_vovv,f_o,f_v)
!    enddo
!  enddo
!  !$OMP END DO NOWAIT
!
!  !$OMP CRITICAL
!  energy = energy + e
!  !$OMP END CRITICAL
!
!  deallocate(W_abc, W_cab, W_bca, W_bac, W_cba, W_acb, &
!             V_abc, V_cab, V_bca, V_bac, V_cba, V_acb )
!
!  !$OMP END PARALLEL

  deallocate(X_vovv,X_ooov,T_voov,T_oovv)
end


integer*8 function find_sample(v, w, n)
  implicit none
  BEGIN_DOC
! Finds sample v in weights w
  END_DOC
  integer*8, intent(in) :: n
  double precision, intent(in) :: v, w(0:n)
  integer*8 :: i,l,r

  l=0
  r=n

  do while(r-l > 1)
    i = shiftr(r+l,1)
    if(w(i) < v) then
      l = i
    else
      r = i
    end if
  end do
  i = r
  do r=i+1,n
    if (w(r) /= w(i)) then
      exit
    endif
  enddo
  find_sample = r-1
end function

