use bitmasks

 BEGIN_PROVIDER [ integer(bit_kind), alphasIcfg_list , (N_int,2,N_configuration,mo_num*(mo_num))]
&BEGIN_PROVIDER [ integer, NalphaIcfg_list, (N_configuration) ]
  implicit none
  !use bitmasks
  BEGIN_DOC
  ! Documentation for alphasI
  ! Returns the associated alpha's for
  ! the input configuration Icfg.
  END_DOC

  integer                        :: idxI ! The id of the Ith CFG
  integer(bit_kind)              :: Icfg(N_int,2)
  integer                        :: NalphaIcfg
  logical,dimension(:,:),allocatable :: tableUniqueAlphas
  integer                        :: listholes(mo_num)
  integer                        :: holetype(mo_num) ! 1-> SOMO 2->DOMO
  integer                        :: nholes
  integer                        :: nvmos
  integer                        :: listvmos(mo_num)
  integer                        :: vmotype(mo_num) ! 1 -> VMO 2 -> SOMO
  integer*8                      :: Idomo
  integer*8                      :: Isomo
  integer*8                      :: Jdomo
  integer*8                      :: Jsomo
  integer*8                      :: diffSOMO
  integer*8                      :: diffDOMO
  integer*8                      :: xordiffSOMODOMO
  integer                        :: ndiffSOMO
  integer                        :: ndiffDOMO
  integer                        :: nxordiffSOMODOMO
  integer                        :: ndiffAll
  integer                        :: i,ii
  integer                        :: j,jj
  integer                        :: k,kk
  integer                        :: kstart
  integer                        :: kend
  integer                        :: Nsomo_I
  integer                        :: hole
  integer                        :: p
  integer                        :: q
  integer                        :: countalphas
  logical                        :: pqAlreadyGenQ
  logical                        :: pqExistsQ
  logical                        :: ppExistsQ
  integer*8                      :: MS

  double precision               :: t0, t1
  call wall_time(t0)

  MS = elec_alpha_num-elec_beta_num

  allocate(tableUniqueAlphas(mo_num,mo_num))
  NalphaIcfg_list = 0

  do idxI = 1, N_configuration

    Icfg  = psi_configuration(:,:,idxI)

    Isomo = iand(act_bitmask(1,1),Icfg(1,1))
    Idomo = iand(act_bitmask(1,1),Icfg(1,2))

    ! find out all pq holes possible
    nholes = 0
    ! holes in SOMO
    do ii = 1,n_act_orb
      i = list_act(ii)
      if(POPCNT(IAND(Isomo,IBSET(0_8,i-1))) .EQ. 1) then
        nholes += 1
        listholes(nholes) = i
        holetype(nholes) = 1
      endif
    end do
    ! holes in DOMO
    do ii = 1,n_act_orb
      i = list_act(ii)
      if(POPCNT(IAND(Idomo,IBSET(0_8,i-1))) .EQ. 1) then
        nholes += 1
        listholes(nholes) = i
        holetype(nholes) = 2
      endif
    end do

    ! find vmos
    listvmos = -1
    vmotype = -1
    nvmos = 0
    do ii = 1,n_act_orb
      i = list_act(ii)
      if(IAND(Idomo,(IBSET(0_8,i-1))) .EQ. 0) then
        if(IAND(Isomo,(IBSET(0_8,i-1))) .EQ. 0) then
          nvmos += 1
          listvmos(nvmos) = i
          vmotype(nvmos) = 1
        else if(POPCNT(IAND(Isomo,(IBSET(0_8,i-1)))) .EQ. 1) then
          nvmos += 1
          listvmos(nvmos) = i
          vmotype(nvmos) = 2
        end if
      end if
    end do

    tableUniqueAlphas = .FALSE.

    ! Now find the allowed (p,q) excitations
    Isomo = iand(act_bitmask(1,1),Icfg(1,1))
    Idomo = iand(act_bitmask(1,1),Icfg(1,2))
    Nsomo_I = POPCNT(Isomo)
    if(Nsomo_I .EQ. 0) then
      kstart = 1
    else
      kstart = cfg_seniority_index(max(NSOMOMin,Nsomo_I-2))
    endif
    kend = idxI-1

    do i = 1,nholes
      p = listholes(i)
      do j = 1,nvmos
        q = listvmos(j)
        if(p .EQ. q) cycle
        if(holetype(i) .EQ. 1 .AND. vmotype(j) .EQ. 1) then
          ! SOMO -> VMO
          Jsomo = IBCLR(Isomo,p-1)
          Jsomo = IBSET(Jsomo,q-1)
          Jdomo = Idomo
          kstart = max(1,cfg_seniority_index(max(NSOMOMin,Nsomo_I-2)))
          kend = idxI-1
        else if(holetype(i) .EQ. 1 .AND. vmotype(j) .EQ. 2) then
          ! SOMO -> SOMO
          Jsomo = IBCLR(Isomo,p-1)
          Jsomo = IBCLR(Jsomo,q-1)
          Jdomo = IBSET(Idomo,q-1)
          ! Check for Minimal alpha electrons (MS)
          if(POPCNT(Jsomo).ge.MS)then
            kstart = max(1,cfg_seniority_index(max(NSOMOMin,Nsomo_I-4)))
            kend = idxI-1
          else
            cycle
          endif
        else if(holetype(i) .EQ. 2 .AND. vmotype(j) .EQ. 1) then
          ! DOMO -> VMO
          Jsomo = IBSET(Isomo,p-1)
          Jsomo = IBSET(Jsomo,q-1)
          Jdomo = IBCLR(Idomo,p-1)
          kstart = cfg_seniority_index(Nsomo_I)
          kend = idxI-1
        else if(holetype(i) .EQ. 2 .AND. vmotype(j) .EQ. 2) then
          ! DOMO -> SOMO
          Jsomo = IBSET(Isomo,p-1)
          Jsomo = IBCLR(Jsomo,q-1)
          Jdomo = IBCLR(Idomo,p-1)
          Jdomo = IBSET(Jdomo,q-1)
          kstart = max(1,cfg_seniority_index(max(NSOMOMin,Nsomo_I-2)))
          kend = idxI-1
        else
          print*,"Something went wrong in obtain_associated_alphaI"
        endif
        ! Check for Minimal alpha electrons (MS)
        if(POPCNT(Jsomo).lt.MS)then
          cycle
        endif

        ! Again, we don't have to search from 1
        ! we just use seniority to find the
        ! first index with NSOMO - 2 to NSOMO + 2
        ! this is what is done in kstart, kend

        pqAlreadyGenQ = .FALSE.
        ! First check if it can be generated before
        do k = kstart, kend
          diffSOMO = IEOR(Jsomo,iand(reunion_of_act_virt_bitmask(1,1),psi_configuration(1,1,k)))
          ndiffSOMO = POPCNT(diffSOMO)
          if((ndiffSOMO .NE. 0) .AND. (ndiffSOMO .NE. 2)) cycle
          diffDOMO = IEOR(Jdomo,iand(reunion_of_act_virt_bitmask(1,1),psi_configuration(1,2,k)))
          xordiffSOMODOMO = IEOR(diffSOMO,diffDOMO)
          ndiffDOMO = POPCNT(diffDOMO)
          nxordiffSOMODOMO = POPCNT(xordiffSOMODOMO)
          nxordiffSOMODOMO += ndiffSOMO + ndiffDOMO
          !if(POPCNT(IEOR(diffSOMO,diffDOMO)) .LE. 1 .AND. ndiffDOMO .LT. 3) then
          if((ndiffSOMO+ndiffDOMO) .EQ. 0) then
            pqAlreadyGenQ = .TRUE.
            ppExistsQ = .TRUE.
            EXIT
          endif
          if((nxordiffSOMODOMO .EQ. 4) .AND. ndiffSOMO .EQ. 2) then
            pqAlreadyGenQ = .TRUE.
            EXIT
          endif
        end do

        if(pqAlreadyGenQ) cycle

        pqExistsQ = .FALSE.

        if(.NOT. pqExistsQ) then
          tableUniqueAlphas(p,q) = .TRUE.
        endif
      end do
    end do

    !print *,tableUniqueAlphas(:,:)

    ! prune list of alphas
    Isomo = Icfg(1,1)
    Idomo = Icfg(1,2)
    Jsomo = Icfg(1,1)
    Jdomo = Icfg(1,2)
    NalphaIcfg = 0
    do i = 1, nholes
      p = listholes(i)
      do j = 1, nvmos
        q = listvmos(j)
        if(p .EQ. q) cycle
        if(tableUniqueAlphas(p,q)) then
          if(holetype(i) .EQ. 1 .AND. vmotype(j) .EQ. 1) then
            ! SOMO -> VMO
            Jsomo = IBCLR(Isomo,p-1)
            Jsomo = IBSET(Jsomo,q-1)
            Jdomo = Idomo
          else if(holetype(i) .EQ. 1 .AND. vmotype(j) .EQ. 2) then
            ! SOMO -> SOMO
            Jsomo = IBCLR(Isomo,p-1)
            Jsomo = IBCLR(Jsomo,q-1)
            Jdomo = IBSET(Idomo,q-1)
            if(POPCNT(Jsomo).ge.MS)then
              kstart = max(1,cfg_seniority_index(max(NSOMOMin,Nsomo_I-4)))
              kend = idxI-1
            else
              cycle
            endif
          else if(holetype(i) .EQ. 2 .AND. vmotype(j) .EQ. 1) then
            ! DOMO -> VMO
            Jsomo = IBSET(Isomo,p-1)
            Jsomo = IBSET(Jsomo,q-1)
            Jdomo = IBCLR(Idomo,p-1)
          else if(holetype(i) .EQ. 2 .AND. vmotype(j) .EQ. 2) then
            ! DOMO -> SOMO
            Jsomo = IBSET(Isomo,p-1)
            Jsomo = IBCLR(Jsomo,q-1)
            Jdomo = IBCLR(Idomo,p-1)
            Jdomo = IBSET(Jdomo,q-1)
          else
            print*,"Something went wrong in obtain_associated_alphaI"
          endif

          ! SOMO
          !print *,i,j,"|",NalphaIcfg, Jsomo, IOR(Jdomo,ISHFT(1_8,n_core_orb)-1)
          if(POPCNT(Jsomo) .ge. NSOMOMin) then
            NalphaIcfg += 1
            alphasIcfg_list(1,1,idxI,NalphaIcfg) = Jsomo
            alphasIcfg_list(1,2,idxI,NalphaIcfg) = IOR(Jdomo,ISHFT(1_8,n_core_orb)-1)
            NalphaIcfg_list(idxI) = NalphaIcfg
          endif
        endif
      end do
    end do

    ! Check if this Icfg has been previously generated as a mono
    ppExistsQ = .False.
    Isomo = iand(reunion_of_act_virt_bitmask(1,1),Icfg(1,1))
    Idomo = iand(reunion_of_act_virt_bitmask(1,1),Icfg(1,2))
    kstart = max(1,cfg_seniority_index(max(NSOMOMin,Nsomo_I-2)))
    do k = kstart, idxI-1
      diffSOMO = IEOR(Isomo,iand(act_bitmask(1,1),psi_configuration(1,1,k)))
      ndiffSOMO = POPCNT(diffSOMO)
      if (ndiffSOMO /= 2) cycle
      diffDOMO = IEOR(Idomo,iand(act_bitmask(1,1),psi_configuration(1,2,k)))
      xordiffSOMODOMO = IEOR(diffSOMO,diffDOMO)
      ndiffDOMO = POPCNT(diffDOMO)
      nxordiffSOMODOMO = POPCNT(xordiffSOMODOMO)
      if((ndiffSOMO+ndiffDOMO+nxordiffSOMODOMO .EQ. 4)) then
        ppExistsQ = .TRUE.
        EXIT
      endif
    end do
    ! Diagonal part (pp,qq)
    if(nholes > 0 .AND. (.NOT. ppExistsQ))then
      ! SOMO
      if(POPCNT(Jsomo) .ge. NSOMOMin) then
        NalphaIcfg += 1
        alphasIcfg_list(1,1,idxI,NalphaIcfg) = Icfg(1,1)
        alphasIcfg_list(1,2,idxI,NalphaIcfg) = Icfg(1,2)
        NalphaIcfg_list(idxI) = NalphaIcfg
      endif
    endif

    NalphaIcfg = 0
  enddo ! end loop idxI
  call wall_time(t1)
  print *, 'Preparation : ', t1 - t0

END_PROVIDER

  subroutine obtain_associated_alphaI(idxI, Icfg, alphasIcfg, NalphaIcfg)
  implicit none
  use bitmasks
  BEGIN_DOC
  ! Documentation for alphasI
  ! Returns the associated alpha's for
  ! the input configuration Icfg.
  END_DOC

  integer,intent(in)                 :: idxI ! The id of the Ith CFG
  integer(bit_kind),intent(in)       :: Icfg(N_int,2)
  integer,intent(out)                :: NalphaIcfg
  integer(bit_kind),intent(out)      :: alphasIcfg(N_int,2,*)
  logical,dimension(:,:),allocatable :: tableUniqueAlphas
  integer                            :: listholes(mo_num)
  integer                            :: holetype(mo_num) ! 1-> SOMO 2->DOMO
  integer                            :: nholes
  integer                            :: nvmos
  integer                            :: listvmos(mo_num)
  integer                            :: vmotype(mo_num) ! 1 -> VMO 2 -> SOMO
  integer*8                          :: Idomo
  integer*8                          :: Isomo
  integer*8                          :: Jdomo
  integer*8                          :: Jsomo
  integer*8                          :: diffSOMO
  integer*8                          :: diffDOMO
  integer*8                          :: xordiffSOMODOMO
  integer                            :: ndiffSOMO
  integer                            :: ndiffDOMO
  integer                            :: nxordiffSOMODOMO
  integer                            :: ndiffAll
  integer                            :: i, ii
  integer                            :: j, jj
  integer                            :: k, kk
  integer                            :: kstart
  integer                            :: kend
  integer                            :: Nsomo_I
  integer                            :: hole
  integer                            :: p
  integer                            :: q
  integer                            :: countalphas
  logical                            :: pqAlreadyGenQ
  logical                            :: pqExistsQ
  logical                            :: ppExistsQ
  Isomo = iand(act_bitmask(1,1),Icfg(1,1))
  Idomo = iand(act_bitmask(1,1),Icfg(1,2))
  !print*,"Input cfg"
  !call debug_spindet(Isomo,1)
  !call debug_spindet(Idomo,1)

  ! find out all pq holes possible
  nholes = 0
  ! holes in SOMO
  do ii = 1,n_act_orb
    i = list_act(ii)
     if(POPCNT(IAND(Isomo,IBSET(0_8,i-1))) .EQ. 1) then
        nholes += 1
        listholes(nholes) = i
        holetype(nholes) = 1
     endif
  end do
  ! holes in DOMO
  do ii = 1,n_act_orb
    i = list_act(ii)
     if(POPCNT(IAND(Idomo,IBSET(0_8,i-1))) .EQ. 1) then
        nholes += 1
        listholes(nholes) = i
        holetype(nholes) = 2
     endif
  end do

  ! find vmos
  listvmos = -1
  vmotype = -1
  nvmos = 0
  do ii = 1,n_act_orb
    i = list_act(ii)
     !print *,i,IBSET(0,i-1),POPCNT(IAND(Isomo,(IBSET(0_8,i-1)))), POPCNT(IAND(Idomo,(IBSET(0_8,i-1))))
     if(POPCNT(IAND(Isomo,(IBSET(0_8,i-1)))) .EQ. 0 .AND. POPCNT(IAND(Idomo,(IBSET(0_8,i-1)))) .EQ. 0) then
        nvmos += 1
        listvmos(nvmos) = i
        vmotype(nvmos) = 1
     else if(POPCNT(IAND(Isomo,(IBSET(0_8,i-1)))) .EQ. 1 .AND. POPCNT(IAND(Idomo,(IBSET(0_8,i-1)))) .EQ. 0 ) then
        nvmos += 1
        listvmos(nvmos) = i
        vmotype(nvmos) = 2
     end if
  end do

  !print *,"Nvmo=",nvmos
  !print *,listvmos
  !print *,vmotype

  allocate(tableUniqueAlphas(mo_num,mo_num))
  tableUniqueAlphas = .FALSE.

  ! Now find the allowed (p,q) excitations
  Isomo = iand(act_bitmask(1,1),Icfg(1,1))
  Idomo = iand(act_bitmask(1,1),Icfg(1,2))
  Nsomo_I = POPCNT(Isomo)
  if(Nsomo_I .EQ. 0) then
    kstart = 1
  else
    kstart = cfg_seniority_index(max(NSOMOMin,Nsomo_I-2))
  endif
  kend = idxI-1
  !print *,"Isomo"
  !call debug_spindet(Isomo,1)
  !call debug_spindet(Idomo,1)

  !print *,"Nholes=",nholes," Nvmos=",nvmos, " idxi=",idxI
  !do i = 1,nholes
  !   print *,i,"->",listholes(i)
  !enddo
  !do i = 1,nvmos
  !   print *,i,"->",listvmos(i)
  !enddo

  do i = 1,nholes
     p = listholes(i)
     do j = 1,nvmos
        q = listvmos(j)
        if(p .EQ. q) cycle
        if(holetype(i) .EQ. 1 .AND. vmotype(j) .EQ. 1) then
           ! SOMO -> VMO
           Jsomo = IBCLR(Isomo,p-1)
           Jsomo = IBSET(Jsomo,q-1)
           Jdomo = Idomo
           kstart = max(1,cfg_seniority_index(max(NSOMOMin,Nsomo_I-2)))
           kend = idxI-1
        else if(holetype(i) .EQ. 1 .AND. vmotype(j) .EQ. 2) then
           ! SOMO -> SOMO
           Jsomo = IBCLR(Isomo,p-1)
           Jsomo = IBCLR(Jsomo,q-1)
           Jdomo = IBSET(Idomo,q-1)
           kstart = max(1,cfg_seniority_index(max(NSOMOMin,Nsomo_I-4)))
           kend = idxI-1
        else if(holetype(i) .EQ. 2 .AND. vmotype(j) .EQ. 1) then
           ! DOMO -> VMO
           Jsomo = IBSET(Isomo,p-1)
           Jsomo = IBSET(Jsomo,q-1)
           Jdomo = IBCLR(Idomo,p-1)
           kstart = cfg_seniority_index(Nsomo_I)
           kend = idxI-1
        else if(holetype(i) .EQ. 2 .AND. vmotype(j) .EQ. 2) then
           ! DOMO -> SOMO
           Jsomo = IBSET(Isomo,p-1)
           Jsomo = IBCLR(Jsomo,q-1)
           Jdomo = IBCLR(Idomo,p-1)
           Jdomo = IBSET(Jdomo,q-1)
           kstart = max(1,cfg_seniority_index(max(NSOMOMin,Nsomo_I-2)))
           kend = idxI-1
        else
           print*,"Something went wrong in obtain_associated_alphaI"
        endif

        ! Again, we don't have to search from 1
        ! we just use seniortiy to find the
        ! first index with NSOMO - 2 to NSOMO + 2
        ! this is what is done in kstart, kend

        pqAlreadyGenQ = .FALSE.
        ! First check if it can be generated before
        do k = kstart, kend
           diffSOMO = IEOR(Jsomo,iand(act_bitmask(1,1),psi_configuration(1,1,k)))
           ndiffSOMO = POPCNT(diffSOMO)
           if((ndiffSOMO .NE. 0) .AND. (ndiffSOMO .NE. 2)) cycle
           diffDOMO = IEOR(Jdomo,iand(act_bitmask(1,1),psi_configuration(1,2,k)))
           xordiffSOMODOMO = IEOR(diffSOMO,diffDOMO)
           ndiffDOMO = POPCNT(diffDOMO)
           nxordiffSOMODOMO = POPCNT(xordiffSOMODOMO)
           nxordiffSOMODOMO += ndiffSOMO + ndiffDOMO
           !if(POPCNT(IEOR(diffSOMO,diffDOMO)) .LE. 1 .AND. ndiffDOMO .LT. 3) then
           if((ndiffSOMO+ndiffDOMO) .EQ. 0) then
              pqAlreadyGenQ = .TRUE.
              ppExistsQ = .TRUE.
              EXIT
           endif
           if((nxordiffSOMODOMO .EQ. 4) .AND. ndiffSOMO .EQ. 2) then
              pqAlreadyGenQ = .TRUE.
              !EXIT
              !ppExistsQ = .TRUE.
              !print *,i,k,ndiffSOMO,ndiffDOMO
              !call debug_spindet(Jsomo,1)
              !call debug_spindet(Jdomo,1)
              !call debug_spindet(iand(reunion_of_act_virt_bitmask(1,1),psi_configuration(1,1,k)),1)
              !call debug_spindet(iand(reunion_of_act_virt_bitmask(1,1),psi_configuration(1,2,k)),1)
              EXIT
           endif
        end do

        !print *,"(,",p,",",q,")",pqAlreadyGenQ

        if(pqAlreadyGenQ) cycle

        pqExistsQ = .FALSE.
        ! now check if this exists in the selected list
        !do k = idxI+1, N_configuration
        !   diffSOMO = IEOR(OR(reunion_of_act_virt_bitmask(1,1),Jsomo),psi_configuration(1,1,k))
        !   diffDOMO = IEOR(OR(reunion_of_act_virt_bitmask(1,1),Jdomo),psi_configuration(1,2,k))
        !   ndiffSOMO = POPCNT(diffSOMO)
        !   ndiffDOMO = POPCNT(diffDOMO)
        !   if((ndiffSOMO + ndiffDOMO) .EQ. 0) then
        !      pqExistsQ = .TRUE.
        !      EXIT
        !   endif
        !end do

        if(.NOT. pqExistsQ) then
           tableUniqueAlphas(p,q) = .TRUE.
           !print *,p,q
           !call debug_spindet(Jsomo,1)
           !call debug_spindet(Jdomo,1)
        endif
     end do
  end do

  !print *,tableUniqueAlphas(:,:)

  ! prune list of alphas
  Isomo = Icfg(1,1)
  Idomo = Icfg(1,2)
  Jsomo = Icfg(1,1)
  Jdomo = Icfg(1,2)
  NalphaIcfg = 0
  do i = 1, nholes
     p = listholes(i)
     do j = 1, nvmos
        q = listvmos(j)
        if(p .EQ. q) cycle
        if(tableUniqueAlphas(p,q)) then
           if(holetype(i) .EQ. 1 .AND. vmotype(j) .EQ. 1) then
              ! SOMO -> VMO
              Jsomo = IBCLR(Isomo,p-1)
              Jsomo = IBSET(Jsomo,q-1)
              Jdomo = Idomo
           else if(holetype(i) .EQ. 1 .AND. vmotype(j) .EQ. 2) then
              ! SOMO -> SOMO
              Jsomo = IBCLR(Isomo,p-1)
              Jsomo = IBCLR(Jsomo,q-1)
              Jdomo = IBSET(Idomo,q-1)
           else if(holetype(i) .EQ. 2 .AND. vmotype(j) .EQ. 1) then
              ! DOMO -> VMO
              Jsomo = IBSET(Isomo,p-1)
              Jsomo = IBSET(Jsomo,q-1)
              Jdomo = IBCLR(Idomo,p-1)
           else if(holetype(i) .EQ. 2 .AND. vmotype(j) .EQ. 2) then
              ! DOMO -> SOMO
              Jsomo = IBSET(Isomo,p-1)
              Jsomo = IBCLR(Jsomo,q-1)
              Jdomo = IBCLR(Idomo,p-1)
              Jdomo = IBSET(Jdomo,q-1)
           else
              print*,"Something went wrong in obtain_associated_alphaI"
           endif

           ! SOMO
           NalphaIcfg += 1
           !print *,i,j,"|",NalphaIcfg
           alphasIcfg(1,1,NalphaIcfg) = Jsomo
           alphasIcfg(1,2,NalphaIcfg) = IOR(Jdomo,ISHFT(1_8,n_core_orb)-1)
           !print *,"I = ",idxI, " Na=",NalphaIcfg," - ",Jsomo, IOR(Jdomo,ISHFT(1_8,n_core_orb)-1)
        endif
     end do
  end do

  ! Check if this Icfg has been previously generated as a mono
  ppExistsQ = .False.
  Isomo = iand(act_bitmask(1,1),Icfg(1,1))
  Idomo = iand(act_bitmask(1,1),Icfg(1,2))
  do k = 1, idxI-1
     diffSOMO = IEOR(Isomo,iand(act_bitmask(1,1),psi_configuration(1,1,k)))
     diffDOMO = IEOR(Idomo,iand(act_bitmask(1,1),psi_configuration(1,2,k)))
     xordiffSOMODOMO = IEOR(diffSOMO,diffDOMO)
     ndiffSOMO = POPCNT(diffSOMO)
     ndiffDOMO = POPCNT(diffDOMO)
     nxordiffSOMODOMO = POPCNT(xordiffSOMODOMO)
     if((ndiffSOMO+ndiffDOMO+nxordiffSOMODOMO .EQ. 4) .AND. ndiffSOMO .EQ. 2) then
        ppExistsQ = .TRUE.
        EXIT
     endif
  end do
  ! Diagonal part (pp,qq)
  if(nholes > 0 .AND. (.NOT. ppExistsQ))then
     ! SOMO
     NalphaIcfg += 1
     !print *,p,q,"|",holetype(i),vmotype(j),NalphaIcfg
     !call debug_spindet(Idomo,1)
     !call debug_spindet(Jdomo,1)
     alphasIcfg(1,1,NalphaIcfg) = Icfg(1,1)
     alphasIcfg(1,2,NalphaIcfg) = Icfg(1,2)
  endif

  end subroutine

  function getNSOMO(Icfg) result(NSOMO)
    implicit none
    integer(bit_kind),intent(in)   :: Icfg(N_int,2)
    integer                        :: NSOMO
    integer                        :: i
    NSOMO = 0
    do i = 1,N_int
       NSOMO += POPCNT(Icfg(i,1))
    enddo
  end function getNSOMO
