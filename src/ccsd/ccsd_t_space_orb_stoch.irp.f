! Main
subroutine ccsd_par_t_space_stoch(nO,nV,t1,t2,f_o,f_v,v_vvvo,v_vvoo,v_vooo,energy)

  implicit none

  integer, intent(in)           :: nO,nV
  double precision, intent(in)  :: t1(nO,nV), f_o(nO), f_v(nV)
  double precision, intent(in)  :: t2(nO,nO,nV,nV)
  double precision, intent(in)  :: v_vvvo(nV,nV,nV,nO), v_vvoo(nV,nV,nO,nO), v_vooo(nV,nO,nO,nO)
  double precision, intent(inout) :: energy

  double precision, allocatable :: X_vovv(:,:,:,:), X_ooov(:,:,:,:), X_oovv(:,:,:,:)
  double precision, allocatable :: T_voov(:,:,:,:), T_oovv(:,:,:,:)
  integer                       :: i,j,k,l,a,b,c,d
  double precision              :: e,ta,tb,eccsd

  eccsd = energy
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

  !$OMP BARRIER
  !$OMP END PARALLEL

  double precision, external :: ccsd_t_task_aba
  double precision, external :: ccsd_t_task_abc
!  logical, external :: omp_test_lock

  double precision, allocatable :: memo(:), Pabc(:), waccu(:)
  integer*8, allocatable :: sampled(:)
!  integer(omp_lock_kind), allocatable :: lock(:)
  integer*2       , allocatable :: abc(:,:)
  integer*8                     :: Nabc, i8,kiter
  integer*8, allocatable :: iorder(:)
  double precision :: eocc
  double precision :: norm
  integer :: isample


  ! Prepare table of triplets (a,b,c)

  Nabc = (int(nV,8) * int(nV+1,8) * int(nV+2,8))/6_8 - nV
  allocate (memo(Nabc), sampled(Nabc), Pabc(Nabc), waccu(0:Nabc))
  allocate (abc(4,Nabc), iorder(Nabc)) !, lock(Nabc))

!  eocc = 3.d0/dble(nO) * sum(f_o(1:nO))
  Nabc = 0_8
  do a = 1, nV
    do b = a+1, nV
      do c = b+1, nV
        Nabc = Nabc + 1_8
        Pabc(Nabc) = -1.d0/(f_v(a) + f_v(b) + f_v(c))
        abc(1,Nabc) = int(a,2)
        abc(2,Nabc) = int(b,2)
        abc(3,Nabc) = int(c,2)
      enddo

      Nabc = Nabc + 1_8
      abc(1,Nabc) = int(a,2)
      abc(2,Nabc) = int(b,2)
      abc(3,Nabc) = int(a,2)
      Pabc(Nabc) = -1.d0/(2.d0*f_v(a) + f_v(b))

      Nabc = Nabc + 1_8
      abc(1,Nabc) = int(b,2)
      abc(2,Nabc) = int(a,2)
      abc(3,Nabc) = int(b,2)
      Pabc(Nabc) = -1.d0/(f_v(a) + 2.d0*f_v(b))
    enddo
  enddo

  do i8=1,Nabc
   iorder(i8) = i8
  enddo

  ! Sort triplets in decreasing Pabc
  call dsort_big(Pabc, iorder, Nabc)

  ! Normalize
  norm = 0.d0
  do i8=Nabc,1,-1
    norm = norm + Pabc(i8)
  enddo
  norm = 1.d0/norm
  do i8=1,Nabc
    Pabc(i8) = Pabc(i8) * norm
  enddo

  call i8set_order_big(abc, iorder, Nabc)


  ! Cumulative distribution for sampling
  waccu(Nabc) = 0.d0
  do i8=Nabc-1,1,-1
   waccu(i8) = waccu(i8+1) - Pabc(i8+1)
  enddo
  waccu(:) = waccu(:) + 1.d0
  waccu(0) = 0.d0

  logical :: converged, do_comp
  double precision :: eta, variance, error, sample
  double precision :: t00, t01
  integer*8 :: ieta, Ncomputed
  integer*8, external :: binary_search

  integer :: nbuckets
  nbuckets = 100

  double precision, allocatable :: wsum(:)
  allocate(wsum(nbuckets))

  converged = .False.
  Ncomputed = 0_8

  energy = 0.d0
  variance = 0.d0
  memo(:) = 0.d0
  sampled(:) = -1_8

  integer*8 :: ileft, iright, imin
  ileft = 1_8
  iright = Nabc
  integer*8, allocatable :: bounds(:,:)

  allocate (bounds(2,nbuckets))
  do isample=1,nbuckets
    eta = 1.d0/dble(nbuckets) * dble(isample)
    ieta = binary_search(waccu,eta,Nabc)
    bounds(1,isample) = ileft
    bounds(2,isample) = ieta
    ileft = ieta+1
    wsum(isample) = sum( Pabc(bounds(1,isample):bounds(2,isample) ) )
  enddo

  Pabc(:) = 1.d0/Pabc(:)

  print '(A)', ''
  print '(A)', ' ======================= ============== =========='
  print '(A)', '        E(CCSD(T))          Error            %    '
  print '(A)', ' ======================= ============== =========='


  call wall_time(t00)
  imin = 1_8
  !$OMP PARALLEL                                                     &
      !$OMP PRIVATE(ieta,eta,a,b,c,kiter,isample)                    &
      !$OMP DEFAULT(SHARED)

  do kiter=1,Nabc

    !$OMP MASTER
    do while (imin <= Nabc)
      if (sampled(imin)>-1_8) then
        imin = imin+1
      else
        exit
      endif
    enddo

    ! Deterministic part
    if (imin < Nabc) then
      ieta=imin
      sampled(ieta) = 0_8
      a = abc(1,ieta)
      b = abc(2,ieta)
      c = abc(3,ieta)
      Ncomputed += 1_8
      !$OMP TASK DEFAULT(SHARED) FIRSTPRIVATE(a,b,c,ieta)
      if (a/=c) then
        memo(ieta) = ccsd_t_task_abc(a,b,c,nO,nV,t1,T_oovv,T_voov, &
            X_ooov,X_oovv,X_vovv,f_o,f_v) / 3.d0
      else
        memo(ieta) =  ccsd_t_task_aba(a,b,nO,nV,t1,T_oovv,T_voov,  &
            X_ooov,X_oovv,X_vovv,f_o,f_v) / 3.d0
      endif
      !$OMP END TASK
    endif

    ! Stochastic part
    call random_number(eta)
    do isample=1,nbuckets
      if (imin >= bounds(2,isample)) then
        cycle
      endif
      ieta = binary_search(waccu,(eta + dble(isample-1))/dble(nbuckets),Nabc)+1

      if (sampled(ieta) == -1_8) then
        sampled(ieta) = 0_8
        a = abc(1,ieta)
        b = abc(2,ieta)
        c = abc(3,ieta)
        Ncomputed += 1_8
        !$OMP TASK DEFAULT(SHARED) FIRSTPRIVATE(a,b,c,ieta)
        if (a/=c) then
          memo(ieta) = ccsd_t_task_abc(a,b,c,nO,nV,t1,T_oovv,T_voov, &
              X_ooov,X_oovv,X_vovv,f_o,f_v) / 3.d0
        else
          memo(ieta) =  ccsd_t_task_aba(a,b,nO,nV,t1,T_oovv,T_voov,  &
              X_ooov,X_oovv,X_vovv,f_o,f_v) / 3.d0
        endif
        !$OMP END TASK
      endif
      sampled(ieta) = sampled(ieta)+1_8

    enddo

    call wall_time(t01)
    if ((t01-t00 > 1.0d0).or.(imin >= Nabc)) then

      !$OMP TASKWAIT
      call wall_time(t01)
      t00 = t01

      double precision :: ET, ET2
      double precision :: energy_stoch, energy_det
      double precision :: scale
      double precision :: w
      double precision :: tmp
      energy_stoch = 0.d0
      energy_det   = 0.d0
      norm = 0.d0
      scale = 1.d0
      ET = 0.d0
      ET2 = 0.d0


      do isample=1,nbuckets
        if (imin >= bounds(2,isample)) then
          energy_det = energy_det + sum(memo(bounds(1,isample):bounds(2,isample)))
          scale = scale - wsum(isample)
        else
          exit
        endif
      enddo

      isample = min(isample,nbuckets)
      do ieta=bounds(1,isample), Nabc
          w = dble(max(sampled(ieta),0_8))
          tmp = w * memo(ieta) * Pabc(ieta)
          ET = ET + tmp
          ET2 = ET2 + tmp * memo(ieta) * Pabc(ieta)
          norm = norm + w
      enddo
      norm = norm/scale
      if (norm > 0.d0) then
        energy_stoch = ET / norm
        variance = ET2 / norm - energy_stoch*energy_stoch
      endif

      energy = energy_det + energy_stoch

      print '(''   '',F20.8, ''   '', ES12.4,''   '', F8.2,''  '')', eccsd+energy, dsqrt(variance/(norm-1.d0)), 100.*real(Ncomputed)/real(Nabc)
    endif
    !$OMP END MASTER
    if (imin >= Nabc) exit
  enddo

  !$OMP END PARALLEL
  print '(A)', ' ======================= ============== ========== '
  print '(A)', ''

  deallocate(X_vovv)
  deallocate(X_ooov)
  deallocate(T_voov)
  deallocate(T_oovv)
end



integer*8 function binary_search(arr, key, sze)
    implicit none
    BEGIN_DOC
! Searches the key in array arr(1:sze) between l_in and r_in, and returns its index
    END_DOC
    integer*8 :: sze, i, j, mid
    double precision :: arr(0:sze)
    double precision :: key

    if ( key < arr(1) ) then
      binary_search = 0_8
      return 
    end if

    if ( key >= arr(sze) ) then
      binary_search = sze
      return 
    end if

    i = 0_8
    j = sze + 1_8

    do while (.True.)
      mid = (i + j) / 2_8
      if ( key >= arr(mid) ) then
        i = mid
      else
        j = mid
      end if
      if (j-i <= 1_8) then
        binary_search = i
        return
      endif
    end do
end function binary_search

