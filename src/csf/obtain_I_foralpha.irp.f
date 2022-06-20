subroutine obtain_connected_J_givenI(idxI, givenI, connectedI, idxs_connectedI, nconnectedI,ntotalconnectedI)
  implicit none
  use bitmasks
  BEGIN_DOC
  ! Documentation for obtain_connected_I_foralpha
  ! This function returns all those selected configurations
  ! which are connected to the input configuration
  ! givenI by a single excitation.
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
  integer          ,intent(in)             :: idxI
  integer(bit_kind),intent(in)             :: givenI(N_int,2)
  integer(bit_kind),intent(out)            :: connectedI(N_int,2,*)
  integer          ,intent(out)            :: idxs_connectedI(*)
  integer,intent(out)                      :: nconnectedI
  integer,intent(out)                      :: ntotalconnectedI
  integer*8                                :: Idomo
  integer*8                                :: Isomo
  integer*8                                :: Jdomo
  integer*8                                :: Jsomo
  integer*8                                :: IJsomo
  integer*8                                :: diffSOMO
  integer*8                                :: diffDOMO
  integer*8                                :: xordiffSOMODOMO
  integer                                  :: ndiffSOMO
  integer                                  :: ndiffDOMO
  integer                                  :: nxordiffSOMODOMO
  integer :: iii,ii,i,j,k,l,p,q,nsomoJ,nsomoalpha,starti,endi,extyp,nholes
  integer :: listholes(mo_num)
  integer :: holetype(mo_num)
  integer :: end_index
  integer :: Nsomo_I

  ! 
  ! 2 2 1 1 0 0 : 1 1 0 0 0 0
  !               0 0 1 1 0 0
  !
  ! 2 1 1 1 1 0 : 1 0 0 0 0 0
  !               0 1 1 1 1 0
  !xorS           0 1 0 0 1 0 : 2
  !xorD           0 1 0 0 0 0 : 1 
  !xorSD          0 0 0 0 1 0 : 1
  !                           -----
  !                             4
  ! 1 1 1 1 1 1 : 0 0 0 0 0 0
  !               1 1 1 1 1 1
  !               1 1 0 0 1 1 : 4
  !               1 1 0 0 0 0 : 2 
  !               0 0 0 0 1 1 : 2
  !                           -----
  !                             8
  !

  nconnectedI = 0
  ntotalconnectedI = 0
  end_index = N_configuration

  ! Since CFGs are sorted wrt to seniority
  ! we don't have to search the full CFG list
  Isomo = givenI(1,1)
  Idomo = givenI(1,2)
  Nsomo_I = POPCNT(Isomo)
  end_index = min(N_configuration,cfg_seniority_index(min(Nsomo_I+6,elec_num))-1)
  if(end_index .LT. 0) end_index= N_configuration
  !end_index = N_configuration
  !print *,"Start and End = ",idxI, end_index


  p = 0
  q = 0
  do i=idxI,end_index
     !if(.True.) then
     !  nconnectedI += 1
     !  connectedI(:,:,nconnectedI) = psi_configuration(:,:,i)
     !  idxs_connectedI(nconnectedI)=i
     !  cycle
     !endif
     Isomo = givenI(1,1)
     Idomo = givenI(1,2)
     Jsomo = psi_configuration(1,1,i)
     Jdomo = psi_configuration(1,2,i)
     diffSOMO = IEOR(Isomo,Jsomo)
     ndiffSOMO = POPCNT(diffSOMO)
     diffDOMO = IEOR(Idomo,Jdomo)
     xordiffSOMODOMO = IEOR(diffSOMO,diffDOMO)
     ndiffDOMO = POPCNT(diffDOMO)
     nxordiffSOMODOMO = POPCNT(xordiffSOMODOMO)
     nxordiffSOMODOMO += ndiffSOMO + ndiffDOMO 
     if((nxordiffSOMODOMO .EQ. 4) .AND. ndiffSOMO .EQ. 2) then
       !-------
       ! MONO |
       !-------
        nconnectedI += 1
        connectedI(:,:,nconnectedI) = psi_configuration(:,:,i)
        idxs_connectedI(nconnectedI)=i
        ntotalconnectedI += max(1,(psi_config_data(i,2)-psi_config_data(i,1)+1))
     else if((nxordiffSOMODOMO .EQ. 8) .AND. ndiffSOMO .EQ. 4) then
       !----------------------------
       ! DOMO -> VMO + DOMO -> VMO |
       !----------------------------
       nconnectedI += 1
       connectedI(:,:,nconnectedI) = psi_configuration(:,:,i)
       idxs_connectedI(nconnectedI)=i
        ntotalconnectedI += max(1,(psi_config_data(i,2)-psi_config_data(i,1)+1))
     else if((nxordiffSOMODOMO .EQ. 6) .AND. ndiffSOMO .EQ. 2) then
       !----------------------------
       ! DOUBLE
       !----------------------------
       nconnectedI += 1
       connectedI(:,:,nconnectedI) = psi_configuration(:,:,i)
       idxs_connectedI(nconnectedI)=i
        ntotalconnectedI += max(1,(psi_config_data(i,2)-psi_config_data(i,1)+1))
     else if((nxordiffSOMODOMO .EQ. 2) .AND. ndiffSOMO .EQ. 3) then
       !-----------------
       ! DOUBLE
       !-----------------
       nconnectedI += 1
       connectedI(:,:,nconnectedI) = psi_configuration(:,:,i)
       idxs_connectedI(nconnectedI)=i
        ntotalconnectedI += max(1,(psi_config_data(i,2)-psi_config_data(i,1)+1))
     else if((nxordiffSOMODOMO .EQ. 4) .AND. ndiffSOMO .EQ. 0) then
       !-----------------
       ! DOUBLE
       !-----------------
       nconnectedI += 1
       connectedI(:,:,nconnectedI) = psi_configuration(:,:,i)
       idxs_connectedI(nconnectedI)=i
        ntotalconnectedI += max(1,(psi_config_data(i,2)-psi_config_data(i,1)+1))
     else if((ndiffSOMO + ndiffDOMO) .EQ. 0) then
       !--------
       ! I = I |
       !--------
       nconnectedI += 1
       connectedI(:,:,nconnectedI) = psi_configuration(:,:,i)
       idxs_connectedI(nconnectedI)= i
        ! find out all pq holes possible
        nholes = 0
        ! holes in SOMO
        Isomo = psi_configuration(1,1,i)
        Idomo = psi_configuration(1,2,i)
        do iii = 1,n_act_orb
          ii = list_act(iii)
           if(POPCNT(IAND(Isomo,IBSET(0_8,ii-1))) .EQ. 1) then
              nholes += 1
              listholes(nholes) = ii
              holetype(nholes) = 1
           endif
        end do
        ! holes in DOMO
        do iii = 1,n_act_orb
          ii = list_act(iii)
           if(POPCNT(IAND(Idomo,IBSET(0_8,ii-1))) .EQ. 1) then
              nholes += 1
              listholes(nholes) = ii
              holetype(nholes) = 2
           endif
        end do
        ntotalconnectedI += max(1,(psi_config_data(i,2)-psi_config_data(i,1)+1)*nholes)
     endif
  end do

end subroutine obtain_connected_J_givenI

subroutine obtain_connected_I_foralpha(idxI, Ialpha, connectedI, idxs_connectedI, nconnectedI, excitationIds, excitationTypes, diagfactors)
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
  integer          ,intent(in)             :: idxI
  integer(bit_kind),intent(in)             :: Ialpha(N_int,2)
  integer(bit_kind),intent(out)            :: connectedI(N_int,2,*)
  integer          ,intent(out)            :: idxs_connectedI(*)
  integer,intent(out)                      :: nconnectedI
  integer,intent(out)                      :: excitationIds(2,*)
  integer,intent(out)                      :: excitationTypes(*)
  real*8 ,intent(out)                      :: diagfactors(*)
  integer*8                                :: Idomo
  integer*8                                :: Isomo
  integer*8                                :: Jdomo
  integer*8                                :: Jsomo
  integer*8                                :: IJsomo
  integer*8                                :: diffSOMO
  integer*8                                :: diffDOMO
  integer*8                                :: xordiffSOMODOMO
  integer                                  :: ndiffSOMO
  integer                                  :: ndiffDOMO
  integer                        :: nxordiffSOMODOMO
  integer                        :: iii,ii,i,j,k,l,p,q,nsomoJ,nsomoalpha,starti,endi,extyp,nholes
  integer                        :: listholes(mo_num)
  integer                        :: holetype(mo_num)
  integer                        :: end_index
  integer                        :: Nsomo_alpha
  integer*8                      :: MS
  MS = elec_alpha_num-elec_beta_num

  nconnectedI = 0
  end_index = N_configuration

  ! Since CFGs are sorted wrt to seniority
  ! we don't have to search the full CFG list
  Isomo = Ialpha(1,1)
  Idomo = Ialpha(1,2)
  Nsomo_alpha = POPCNT(Isomo)
  end_index = min(N_configuration,cfg_seniority_index(min(Nsomo_alpha+4,elec_num))-1)
  if(end_index .LT. 0) end_index= N_configuration
  end_index = N_configuration


  p = 0
  q = 0
  if (N_int > 1) stop 'obtain_connected_i_foralpha : N_int > 1'
  do i=idxI,end_index
     Isomo = Ialpha(1,1)
     Idomo = Ialpha(1,2)
     Jsomo = psi_configuration(1,1,i)
     Jdomo = psi_configuration(1,2,i)
     ! Check for Minimal alpha electrons (MS)
     if(POPCNT(Isomo).lt.MS)then
       cycle
     endif
     diffSOMO = IEOR(Isomo,Jsomo)
     ndiffSOMO = POPCNT(diffSOMO)
     !if(idxI.eq.1)then
     !  print *," \t idxI=",i," diffS=",ndiffSOMO," popJs=", POPCNT(Jsomo)," popIs=",POPCNT(Isomo)
     !endif
     diffDOMO = IEOR(Idomo,Jdomo)
     xordiffSOMODOMO = IEOR(diffSOMO,diffDOMO)
     ndiffDOMO = POPCNT(diffDOMO)
     nxordiffSOMODOMO = POPCNT(xordiffSOMODOMO)
     nxordiffSOMODOMO += ndiffSOMO + ndiffDOMO 
     if((nxordiffSOMODOMO .EQ. 4) .AND. ndiffSOMO .EQ. 2) then
        select case(ndiffDOMO)
        case (0)
           ! SOMO -> VMO
           !print *,"obt SOMO -> VMO"
           extyp = 3
           IJsomo = IEOR(Isomo, Jsomo)
!IRP_IF WITHOUT_TRAILZ
!           p = (popcnt(ieor( IAND(Isomo,IJsomo) , IAND(Isomo,IJsomo) -1))-1) + 1
!IRP_ELSE
           p = TRAILZ(IAND(Isomo,IJsomo)) + 1
!IRP_ENDIF
           IJsomo = IBCLR(IJsomo,p-1)
!IRP_IF WITHOUT_TRAILZ
!           q = (popcnt(ieor(IJsomo,IJsomo-1))-1) + 1
!IRP_ELSE
           q = TRAILZ(IJsomo) + 1
!IRP_ENDIF
        case (1)
           ! DOMO -> VMO
           ! or
           ! SOMO -> SOMO
           nsomoJ = POPCNT(Jsomo)
           nsomoalpha = POPCNT(Isomo)
           if(nsomoJ .GT. nsomoalpha) then
              ! DOMO -> VMO
              !print *,"obt DOMO -> VMO"
              extyp = 2
!IRP_IF WITHOUT_TRAILZ
!              p = (popcnt(ieor( IEOR(Idomo,Jdomo),IEOR(Idomo,Jdomo) -1))-1) + 1
!IRP_ELSE
              p = TRAILZ(IEOR(Idomo,Jdomo)) + 1
!IRP_ENDIF
              Isomo = IEOR(Isomo, Jsomo)
              Isomo = IBCLR(Isomo,p-1)
!IRP_IF WITHOUT_TRAILZ
!              q = (popcnt(ieor(Isomo,Isomo-1))-1) + 1
!IRP_ELSE
              q = TRAILZ(Isomo) + 1
!IRP_ENDIF
           else
              ! SOMO -> SOMO
              !print *,"obt SOMO -> SOMO"
              extyp = 1
!IRP_IF WITHOUT_TRAILZ
!              q = (popcnt(ieor( IEOR(Idomo,Jdomo), IEOR(Idomo,Jdomo)-1))-1) + 1
!IRP_ELSE
              q = TRAILZ(IEOR(Idomo,Jdomo)) + 1
!IRP_ENDIF
              Isomo = IEOR(Isomo, Jsomo)
              Isomo = IBCLR(Isomo,q-1)
!IRP_IF WITHOUT_TRAILZ
!              p = (popcnt(ieor(Isomo,Isomo-1))-1) + 1
!IRP_ELSE
              p = TRAILZ(Isomo) + 1
!IRP_ENDIF
              ! Check for Minimal alpha electrons (MS)
              !if(POPCNT(Isomo).lt.MS)then
              !  cycle
              !endif
           end if
        case (2)
           ! DOMO -> SOMO
           !print *,"obt DOMO -> SOMO"
           extyp = 4
           IJsomo = IEOR(Isomo, Jsomo)
!IRP_IF WITHOUT_TRAILZ
!           p = (popcnt(ieor( IAND(Jsomo,IJsomo), IAND(Jsomo,IJsomo)-1))-1) + 1
!IRP_ELSE
           p = TRAILZ(IAND(Jsomo,IJsomo)) + 1
!IRP_ENDIF
           IJsomo = IBCLR(IJsomo,p-1)
!IRP_IF WITHOUT_TRAILZ
!           q = (popcnt(ieor( IJsomo , IJsomo -1))-1) + 1
!IRP_ELSE
           q = TRAILZ(IJsomo) + 1
!IRP_ENDIF
        case default
           print *,"something went wront in get connectedI"
        end select
        starti = psi_config_data(i,1)
        endi   = psi_config_data(i,2)
        nconnectedI += 1
        do k=1,N_int
          connectedI(k,1,nconnectedI) = psi_configuration(k,1,i)
          connectedI(k,2,nconnectedI) = psi_configuration(k,2,i)
        enddo
        idxs_connectedI(nconnectedI)=starti
        excitationIds(1,nconnectedI)=p
        excitationIds(2,nconnectedI)=q
        excitationTypes(nconnectedI) = extyp
        diagfactors(nconnectedI) = 1.0d0
     else if((ndiffSOMO + ndiffDOMO) .EQ. 0) then
        ! find out all pq holes possible
        nholes = 0
        ! holes in SOMO
        Isomo = psi_configuration(1,1,i)
        Idomo = psi_configuration(1,2,i)
        do iii = 1,n_act_orb
          ii = list_act(iii)
           if(POPCNT(IAND(Isomo,IBSET(0_8,ii-1))) .EQ. 1) then
              nholes += 1
              listholes(nholes) = ii
              holetype(nholes) = 1
           endif
        end do
        ! holes in DOMO
        do iii = 1,n_act_orb
          ii = list_act(iii)
           if(POPCNT(IAND(Idomo,IBSET(0_8,ii-1))) .EQ. 1) then
              nholes += 1
              listholes(nholes) = ii
              holetype(nholes) = 2
           endif
        end do

        do k=1,nholes
           p = listholes(k)
           q = p
           extyp = 1
           if(holetype(k) .EQ. 1) then
              starti = psi_config_data(i,1)
              endi   = psi_config_data(i,2)
              nconnectedI += 1
              connectedI(:,:,nconnectedI) = psi_configuration(:,:,i)
              idxs_connectedI(nconnectedI)=starti
              excitationIds(1,nconnectedI)=p
              excitationIds(2,nconnectedI)=q
              excitationTypes(nconnectedI) = extyp
              diagfactors(nconnectedI) = 1.0d0
           else
              starti = psi_config_data(i,1)
              endi   = psi_config_data(i,2)
              nconnectedI += 1
              connectedI(:,:,nconnectedI) = psi_configuration(:,:,i)
              idxs_connectedI(nconnectedI)=starti
              excitationIds(1,nconnectedI)=p
              excitationIds(2,nconnectedI)=q
              excitationTypes(nconnectedI) = extyp
              diagfactors(nconnectedI) = 2.0d0
           endif
        enddo
     endif
  end do

end subroutine obtain_connected_I_foralpha
