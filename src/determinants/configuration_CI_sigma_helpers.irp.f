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
  integer                            :: ndiffSOMO
  integer                            :: ndiffDOMO
  integer                            :: ndiffAll
  integer                            :: i
  integer                            :: j
  integer                            :: k
  integer                            :: hole
  integer                            :: p
  integer                            :: q
  integer                            :: countalphas
  logical                            :: pqAlreadyGenQ
  logical                            :: pqExistsQ
  Isomo = iand(reunion_of_act_virt_bitmask(1,1),Icfg(1,1))
  Idomo = iand(reunion_of_act_virt_bitmask(1,1),Icfg(1,2))
  !print*,"Input cfg"
  !call debug_spindet(Isomo,1)
  !call debug_spindet(Idomo,1)

  !print*,n_act_orb, "monum=",mo_num," n_core=",n_core_orb

  ! find out all pq holes possible
  nholes = 0
  ! holes in SOMO
  do i = n_core_orb+1,n_core_orb + n_act_orb
     if(POPCNT(IAND(Isomo,IBSET(0_8,i-1))) .EQ. 1) then
        nholes += 1
        listholes(nholes) = i
        holetype(nholes) = 1
     endif
  end do
  ! holes in DOMO
  do i = n_core_orb+1,n_core_orb + n_act_orb
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
  do i = n_core_orb+1,n_core_orb + n_act_orb
     !print *,i,IBSET(0,i-1),POPCNT(IAND(Isomo,(IBSET(0,i-1)))), POPCNT(IAND(Idomo,(IBSET(0,i-1))))
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
  !print *,"Isomo"
  !call debug_spindet(Isomo,1)
  !call debug_spindet(Idomo,1)
  do i = 1,nholes
     p = listholes(i)
     do j = 1,nvmos
        q = listvmos(j)
        if(p == q) cycle
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


        pqAlreadyGenQ = .FALSE.
        ! First check if it can be generated before
        do k = 1, idxI-1
           diffSOMO = XOR(Jsomo,iand(reunion_of_act_virt_bitmask(1,1),psi_configuration(1,1,k)))
           diffDOMO = XOR(Jdomo,iand(reunion_of_act_virt_bitmask(1,1),psi_configuration(1,2,k)))
           ndiffSOMO = POPCNT(diffSOMO)
           ndiffDOMO = POPCNT(diffDOMO)
           if(POPCNT(XOR(diffSOMO,diffDOMO)) .LE. 1 .AND. ndiffDOMO .LT. 3) then
              pqAlreadyGenQ = .TRUE.
              !print *,i,k,ndiffSOMO,ndiffDOMO
              !call debug_spindet(Jsomo,1)
              !call debug_spindet(Jdomo,1)
              !call debug_spindet(iand(reunion_of_act_virt_bitmask(1,1),psi_configuration(1,1,k)),1)
              !call debug_spindet(iand(reunion_of_act_virt_bitmask(1,1),psi_configuration(1,2,k)),1)
              EXIT
           endif
        end do

        if(pqAlreadyGenQ) cycle

        pqExistsQ = .FALSE.
        ! now check if this exists in the selected list
        do k = idxI, N_configuration
           diffSOMO = XOR(OR(reunion_of_act_virt_bitmask(1,1),Jsomo),psi_configuration(1,1,k))
           diffDOMO = XOR(OR(reunion_of_act_virt_bitmask(1,1),Jdomo),psi_configuration(1,2,k))
           ndiffSOMO = POPCNT(diffSOMO)
           ndiffDOMO = POPCNT(diffDOMO)
           if((ndiffSOMO + ndiffDOMO) .EQ. 0) then
              pqExistsQ = .TRUE.
              EXIT
           endif
        end do

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

           NalphaIcfg += 1
           !print *,p,q,"|",holetype(i),vmotype(j)
           !call debug_spindet(Jsomo,1)
           !call debug_spindet(Jdomo,1)
           alphasIcfg(1,1,NalphaIcfg) = Jsomo
           alphasIcfg(1,2,NalphaIcfg) = Jdomo
        endif
     end do
  end do

  end subroutine

subroutine obtain_connected_I_foralpha(Ialpha, connectedI, nconnectedI, excitationIds, excitationTypes)
  implicit none
  use bitmasks
  BEGIN_DOC
  ! Documentation for obtain_connected_I_foralpha
  ! This function returns all those selected configurations
  ! which are connected to the input configuration
  ! Ialpha by a single excitation.
  !
  ! The type of excitations are ordered as follows:
  ! Type 1 - SOMO -> SOMO
  ! Type 2 - DOMO -> VMO
  ! Type 3 - SOMO -> VMO
  ! Type 4 - DOMO -> SOMO
  !
  ! Order of operators
  ! \alpha> = a^\dag_p a_q |I> = E_pq |I>
  END_DOC
  integer(bit_kind),intent(in)             :: Ialpha(N_int,2)
  integer(bit_kind),intent(out)            :: connectedI(N_int,2,*)
  integer,intent(out)                      :: nconnectedI
  integer,intent(out)                      :: excitationIds(2,*)
  integer,intent(out)                      :: excitationTypes(*)
  integer*8                                :: Idomo
  integer*8                                :: Isomo
  integer*8                                :: Jdomo
  integer*8                                :: Jsomo
  integer*8                                :: diffSOMO
  integer*8                                :: diffDOMO
  integer                                  :: ndiffSOMO
  integer                                  :: ndiffDOMO
  integer :: i,j,k,l,p,q,nsomoI,nsomoalpha

  nconnectedI = 0

  Isomo = Ialpha(1,1)
  Idomo = Ialpha(1,2)
  do i=1,N_configuration
     diffSOMO = XOR(Isomo,psi_configuration(1,1,i))
     diffDOMO = XOR(Idomo,psi_configuration(1,2,i))
     ndiffSOMO = POPCNT(diffSOMO)
     ndiffDOMO = POPCNT(diffDOMO)
     if(POPCNT(XOR(diffSOMO,diffDOMO)) .LE. 1 .AND. ndiffDOMO .LT. 3) then
        nconnectedI += 1
        connectedI(:,:,nconnectedI) = psi_configuration(:,:,i)
        select case(ndiffDOMO)
        case (0)
           ! SOMO -> VMO
           excitationTypes(nconnectedI) = 3
           Jsomo = XOR(Isomo, psi_configuration(1,1,i))
           p = TRAILZ(iand(Jsomo,Isomo))
           q = TRAILZ(iand(Jsomo,psi_configuration(1,1,i)))
        case (1)
           ! DOMO -> VMO
           ! or
           ! SOMO -> SOMO
           nsomoI = POPCNT(psi_configuration(1,1,i))
           nsomoalpha = POPCNT(Isomo)
           if(nsomoI .GT. nsomoalpha) then
              ! DOMO -> VMO
              excitationTypes(nconnectedI) = 2
              Jdomo = XOR(Idomo, psi_configuration(1,2,i))
              Jsomo = XOR(Jdomo,XOR(Isomo, psi_configuration(1,1,i)))
              p = TRAILZ(Jdomo)
              q = TRAILZ(Jsomo)
           else
              ! SOMO -> SOMO
              excitationTypes(nconnectedI) = 1
              Jdomo = XOR(Idomo, psi_configuration(1,2,i))
              Jsomo = XOR(Jdomo,XOR(Isomo, psi_configuration(1,1,i)))
              q = TRAILZ(Jdomo)
              p = TRAILZ(Jsomo)
           end if
        case (2)
           ! DOMO -> SOMO
           excitationTypes(nconnectedI) = 4
           Jdomo = XOR(Idomo, psi_configuration(1,2,i))
           p = TRAILZ(iand(Jdomo,Idomo))
           q = TRAILZ(iand(Jdomo,psi_configuration(1,2,i)))
        case default
           print *,"something went wront in get connectedI"
        end select
        excitationIds(1,nconnectedI)=p
        excitationIds(2,nconnectedI)=q
     endif
  end do


end subroutine obtain_connected_I_foralpha
