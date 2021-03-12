use bitmasks

 BEGIN_PROVIDER [ integer(bit_kind), alphasIcfg_list , (N_int,2,N_configuration,mo_num*mo_num)]
&BEGIN_PROVIDER [ integer, NalphaIcfg_list, (N_configuration) ]
  implicit none
  !use bitmasks
  BEGIN_DOC
  ! Documentation for alphasI
  ! Returns the associated alpha's for
  ! the input configuration Icfg.
  END_DOC

  integer                            :: idxI ! The id of the Ith CFG
  integer(bit_kind)                  :: Icfg(N_int,2)
  integer                            :: NalphaIcfg
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
  integer                            :: i
  integer                            :: j
  integer                            :: k
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


  allocate(tableUniqueAlphas(mo_num,mo_num))

  do idxI = 1, N_configuration

    Icfg  = psi_configuration(:,:,idxI)

    Isomo = iand(reunion_of_act_virt_bitmask(1,1),Icfg(1,1))
    Idomo = iand(reunion_of_act_virt_bitmask(1,1),Icfg(1,2))

    ! find out all pq holes possible
    nholes = 0
    ! holes in SOMO
    do i = 1,mo_num
       if(POPCNT(IAND(Isomo,IBSET(0_8,i-1))) .EQ. 1) then
          nholes += 1
          listholes(nholes) = i
          holetype(nholes) = 1
       endif
    end do
    ! holes in DOMO
    do i = 1,mo_num
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
    do i = 1,mo_num
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

    tableUniqueAlphas = .FALSE.

    ! Now find the allowed (p,q) excitations
    Isomo = iand(reunion_of_act_virt_bitmask(1,1),Icfg(1,1))
    Idomo = iand(reunion_of_act_virt_bitmask(1,1),Icfg(1,2))
    Nsomo_I = POPCNT(Isomo)
    if(Nsomo_I .EQ. 0) then
      kstart = 1
    else
      kstart = cfg_seniority_index(max(0,Nsomo_I-2))
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
             kstart = max(1,cfg_seniority_index(max(0,Nsomo_I-2)))
             kend = idxI-1
          else if(holetype(i) .EQ. 1 .AND. vmotype(j) .EQ. 2) then
             ! SOMO -> SOMO
             Jsomo = IBCLR(Isomo,p-1)
             Jsomo = IBCLR(Jsomo,q-1)
             Jdomo = IBSET(Idomo,q-1)
             kstart = max(1,cfg_seniority_index(max(0,Nsomo_I-4)))
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
             kstart = max(1,cfg_seniority_index(max(0,Nsomo_I-2)))
             kend = idxI-1
          else
             print*,"Something went wrong in obtain_associated_alphaI"
          endif

          ! Again, we don't have to search from 1
          ! we just use seniorty to find the
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
           alphasIcfg_list(1,1,idxI,NalphaIcfg) = Jsomo
           alphasIcfg_list(1,2,idxI,NalphaIcfg) = IOR(Jdomo,ISHFT(1_8,n_core_orb)-1)
           NalphaIcfg_list(idxI) = NalphaIcfg
        endif
     end do
  end do

  ! Check if this Icfg has been previously generated as a mono
  ppExistsQ = .False.
  Isomo = iand(reunion_of_act_virt_bitmask(1,1),Icfg(1,1))
  Idomo = iand(reunion_of_act_virt_bitmask(1,1),Icfg(1,2))
  do k = 1, idxI-1
     diffSOMO = IEOR(Isomo,iand(reunion_of_act_virt_bitmask(1,1),psi_configuration(1,1,k)))
     diffDOMO = IEOR(Idomo,iand(reunion_of_act_virt_bitmask(1,1),psi_configuration(1,2,k)))
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
     alphasIcfg_list(1,1,idxI,NalphaIcfg) = Icfg(1,1)
     alphasIcfg_list(1,2,idxI,NalphaIcfg) = Icfg(1,2)
     NalphaIcfg_list(idxI) = NalphaIcfg
  endif


  enddo ! end loop idxI

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
  integer                            :: i
  integer                            :: j
  integer                            :: k
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
  Isomo = iand(reunion_of_act_virt_bitmask(1,1),Icfg(1,1))
  Idomo = iand(reunion_of_act_virt_bitmask(1,1),Icfg(1,2))
  !print*,"Input cfg"
  !call debug_spindet(Isomo,1)
  !call debug_spindet(Idomo,1)

  !print*,n_act_orb, "monum=",mo_num," n_core=",n_core_orb

  ! find out all pq holes possible
  nholes = 0
  ! holes in SOMO
  do i = 1,mo_num
     if(POPCNT(IAND(Isomo,IBSET(0_8,i-1))) .EQ. 1) then
        nholes += 1
        listholes(nholes) = i
        holetype(nholes) = 1
     endif
  end do
  ! holes in DOMO
  do i = 1,mo_num
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
  do i = 1,mo_num
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
  Isomo = iand(reunion_of_act_virt_bitmask(1,1),Icfg(1,1))
  Idomo = iand(reunion_of_act_virt_bitmask(1,1),Icfg(1,2))
  Nsomo_I = POPCNT(Isomo)
  if(Nsomo_I .EQ. 0) then
    kstart = 1
  else
    kstart = cfg_seniority_index(max(0,Nsomo_I-2))
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
           kstart = max(1,cfg_seniority_index(max(0,Nsomo_I-2)))
           kend = idxI-1
        else if(holetype(i) .EQ. 1 .AND. vmotype(j) .EQ. 2) then
           ! SOMO -> SOMO
           Jsomo = IBCLR(Isomo,p-1)
           Jsomo = IBCLR(Jsomo,q-1)
           Jdomo = IBSET(Idomo,q-1)
           kstart = max(1,cfg_seniority_index(max(0,Nsomo_I-4)))
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
           kstart = max(1,cfg_seniority_index(max(0,Nsomo_I-2)))
           kend = idxI-1
        else
           print*,"Something went wrong in obtain_associated_alphaI"
        endif

        ! Again, we don't have to search from 1
        ! we just use seniorty to find the
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
  Isomo = iand(reunion_of_act_virt_bitmask(1,1),Icfg(1,1))
  Idomo = iand(reunion_of_act_virt_bitmask(1,1),Icfg(1,2))
  do k = 1, idxI-1
     diffSOMO = IEOR(Isomo,iand(reunion_of_act_virt_bitmask(1,1),psi_configuration(1,1,k)))
     diffDOMO = IEOR(Idomo,iand(reunion_of_act_virt_bitmask(1,1),psi_configuration(1,2,k)))
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

subroutine convertOrbIdsToModelSpaceIds(Ialpha, Jcfg, p, q, extype, pmodel, qmodel)
  implicit none
  BEGIN_DOC
  ! This function converts the orbital ids
  ! in real space to those used in model space
  ! in order to identify the matrices required
  ! for the calculation of MEs.
  !
  ! The type of excitations are ordered as follows:
  ! Type 1 - SOMO -> SOMO
  ! Type 2 - DOMO -> VMO
  ! Type 3 - SOMO -> VMO
  ! Type 4 - DOMO -> SOMO
  END_DOC
  integer(bit_kind),intent(in)   :: Ialpha(N_int,2)
  integer(bit_kind),intent(in)   :: Jcfg(N_int,2)
  integer,intent(in)             :: p,q
  integer,intent(in)             :: extype
  integer,intent(out)            :: pmodel,qmodel
  !integer(bit_kind)              :: Isomo(N_int)
  !integer(bit_kind)              :: Idomo(N_int)
  !integer(bit_kind)              :: Jsomo(N_int)
  !integer(bit_kind)              :: Jdomo(N_int)
  integer*8                       :: Isomo       
  integer*8                       :: Idomo       
  integer*8                       :: Jsomo       
  integer*8                       :: Jdomo       
  integer*8                      :: mask
  integer                        :: iint, ipos
  !integer(bit_kind)              :: Isomotmp(N_int)
  !integer(bit_kind)              :: Jsomotmp(N_int)
  integer*8             :: Isomotmp
  integer*8             :: Jsomotmp
  integer                        :: pos0,pos0prev

  ! TODO Flag (print) when model space indices is > 64
  Isomo = Ialpha(1,1)
  Idomo = Ialpha(1,2)
  Jsomo = Jcfg(1,1)
  Jdomo = Jcfg(1,2)
  pos0prev = 0
  pmodel = p
  qmodel = q

  if(p .EQ. q) then
     pmodel = 1
     qmodel = 1
  else
     select case(extype)
       case (1)
          ! SOMO -> SOMO
          ! remove all domos
          !print *,"type -> SOMO -> SOMO"
          mask = ISHFT(1_8,p) - 1
          Isomotmp = IAND(Isomo,mask)
          pmodel = POPCNT(mask) - POPCNT(XOR(Isomotmp,mask))
          mask = ISHFT(1_8,q) - 1
          Isomotmp = IAND(Isomo,mask)
          qmodel = POPCNT(mask) - POPCNT(XOR(Isomotmp,mask))
       case (2)
          ! DOMO -> VMO
          ! remove all domos except one at p
          !print *,"type -> DOMO -> VMO"
          mask = ISHFT(1_8,p) - 1
          Jsomotmp = IAND(Jsomo,mask)
          pmodel = POPCNT(mask) - POPCNT(XOR(Jsomotmp,mask))
          mask = ISHFT(1_8,q) - 1
          Jsomotmp = IAND(Jsomo,mask)
          qmodel = POPCNT(mask) - POPCNT(XOR(Jsomotmp,mask))
       case (3)
          ! SOMO -> VMO
          !print *,"type -> SOMO -> VMO"
          !Isomo = IEOR(Isomo,Jsomo)
          if(p.LT.q) then
             mask = ISHFT(1_8,p) - 1
             Isomo = IAND(Isomo,mask)
             pmodel = POPCNT(mask) - POPCNT(XOR(Isomo,mask))
             mask = ISHFT(1_8,q) - 1
             Jsomo = IAND(Jsomo,mask)
             qmodel = POPCNT(mask) - POPCNT(XOR(Jsomo,mask)) + 1
          else
             mask = ISHFT(1_8,p) - 1
             Isomo = IAND(Isomo,mask)
             pmodel = POPCNT(mask) - POPCNT(XOR(Isomo,mask)) + 1
             mask = ISHFT(1_8,q) - 1
             Jsomo = IAND(Jsomo,mask)
             qmodel = POPCNT(mask) - POPCNT(XOR(Jsomo,mask))
          endif
       case (4)
          ! DOMO -> SOMO
          ! remove all domos except one at p
          !print *,"type -> DOMO -> SOMO"
          !Isomo = IEOR(Isomo,Jsomo)
          if(p.LT.q) then
             mask = ISHFT(1_8,p) - 1
             Jsomo = IAND(Jsomo,mask)
             pmodel = POPCNT(mask) - POPCNT(XOR(Jsomo,mask))
             mask = ISHFT(1_8,q) - 1
             Isomo = IAND(Isomo,mask)
             qmodel = POPCNT(mask) - POPCNT(XOR(Isomo,mask)) + 1
          else
             mask = ISHFT(1_8,p) - 1
             Jsomo = IAND(Jsomo,mask)
             pmodel = POPCNT(mask) - POPCNT(XOR(Jsomo,mask)) + 1
             mask = ISHFT(1_8,q) - 1
             Isomo = IAND(Isomo,mask)
             qmodel = POPCNT(mask) - POPCNT(XOR(Isomo,mask))
          endif
       case default
          print *,"something is wrong in convertOrbIdsToModelSpaceIds"
     end select
  endif
  !print *,p,q,"model ids=",pmodel,qmodel
end subroutine convertOrbIdsToModelSpaceIds
