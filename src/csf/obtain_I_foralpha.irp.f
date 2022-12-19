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
  integer :: listall(N_int*bit_kind_size), nelall

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
  Nsomo_I = 0
  do i=1,N_int
    Isomo = givenI(i,1)
    Idomo = givenI(i,2)
    Nsomo_I += POPCNT(Isomo)
  end do
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

     ndiffSOMO = 0
     ndiffDOMO = 0
     nxordiffSOMODOMO = 0
     do ii=1,N_int
       Isomo = givenI(ii,1)
       Idomo = givenI(ii,2)
       Jsomo = psi_configuration(ii,1,i)
       Jdomo = psi_configuration(ii,2,i)
       diffSOMO = IEOR(Isomo,Jsomo)
       ndiffSOMO += POPCNT(diffSOMO)
       diffDOMO = IEOR(Idomo,Jdomo)
       xordiffSOMODOMO = IEOR(diffSOMO,diffDOMO)
       ndiffDOMO += POPCNT(diffDOMO)
       nxordiffSOMODOMO += POPCNT(xordiffSOMODOMO)
       nxordiffSOMODOMO += POPCNT(diffSOMO) + POPCNT(diffDOMO)
     end do

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
        !Isomo = psi_configuration(1,1,i)
        !Idomo = psi_configuration(1,2,i)
        !do iii = 1,n_act_orb
        !  ii = list_act(iii)
        !   if(POPCNT(IAND(Isomo,IBSET(0_8,ii-1))) .EQ. 1) then
        !      nholes += 1
        !      listholes(nholes) = ii
        !      holetype(nholes) = 1
        !   endif
        !end do
       
        call bitstring_to_list(psi_configuration(1,1,i),listall,nelall,N_int)

        do iii=1,nelall
          nholes += 1
          listholes(nholes) = listall(iii)
          holetype(nholes) = 1
        end do

        ! holes in DOMO
        !do iii = 1,n_act_orb
        !  ii = list_act(iii)
        !   if(POPCNT(IAND(Idomo,IBSET(0_8,ii-1))) .EQ. 1) then
        !      nholes += 1
        !      listholes(nholes) = ii
        !      holetype(nholes) = 2
        !   endif
        !end do
        
        call bitstring_to_list(psi_configuration(1,2,i),listall,nelall,N_int)

        do iii=1,nelall
          if(listall(iii) .gt. n_core_orb)then
            nholes += 1
            listholes(nholes) = listall(iii)
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
  integer(bit_kind)                        :: Jcfg(N_int,2)
  integer(bit_kind)                        :: Icfg(N_int,2)
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
  integer                        :: end_index, ishift
  integer                        :: Nsomo_alpha, pp,qq, nperm, iint, ipos
  integer*8                      :: MS
  integer                        :: exc(0:2,2,2), tz, m, n, high, low
  integer                        :: listall(N_int*bit_kind_size), nelall
  integer                        :: nconnectedExtradiag, nconnectedDiag
  integer(bit_kind)              :: hole, particle, tmp
  MS = elec_alpha_num-elec_beta_num

  nconnectedExtradiag=0
  nconnectedDiag=0
  nconnectedI = 0
  end_index = N_configuration

  ! Since CFGs are sorted wrt to seniority
  ! we don't have to search the full CFG list
  !Isomo = Ialpha(1,1)
  !Idomo = Ialpha(1,2)
  !Nsomo_alpha = POPCNT(Isomo)
  Icfg = Ialpha
  Nsomo_alpha = 0
  !print *," Ialpha="
  do ii=1,N_int
    Isomo = Ialpha(ii,1)
    Idomo = Ialpha(ii,2)
    Nsomo_alpha += POPCNT(Isomo)
    !print *,Isomo, Idomo, "Nsomo=",Nsomo_alpha
  end do
  end_index = min(N_configuration,cfg_seniority_index(min(Nsomo_alpha+4,elec_num))-1)
  if(end_index .LT. 0 .OR. end_index .lt. idxI) end_index= N_configuration
  end_index = N_configuration


  p = 0
  q = 0
  !if (N_int > 1) stop 'obtain_connected_i_foralpha : N_int > 1'
  do i=idxI,end_index
     ! Check for Minimal alpha electrons (MS)
     if(Nsomo_alpha .lt. MS)then
       cycle
     endif

     ndiffSOMO = 0
     ndiffDOMO = 0
     nxordiffSOMODOMO = 0
     nsomoJ=0
     nsomoalpha=0
     do ii=1,N_int
       Isomo = Ialpha(ii,1)
       Idomo = Ialpha(ii,2)
       Jsomo = psi_configuration(ii,1,i)
       Jdomo = psi_configuration(ii,2,i)
       nsomoJ += POPCNT(Jsomo)
       nsomoalpha += POPCNT(Isomo)
       diffSOMO = IEOR(Isomo,Jsomo)
       ndiffSOMO += POPCNT(diffSOMO)
       diffDOMO = IEOR(Idomo,Jdomo)
       xordiffSOMODOMO = IEOR(diffSOMO,diffDOMO)
       ndiffDOMO += POPCNT(diffDOMO)
       nxordiffSOMODOMO += POPCNT(xordiffSOMODOMO)
       nxordiffSOMODOMO += POPCNT(diffSOMO) + POPCNT(diffDOMO)
     end do
     !if(idxI.eq.218)then
     !  print *,"I=",idxI,"Nsomo_alpha=",Nsomo_alpha,"nxordiffSOMODOMO(4)=",nxordiffSOMODOMO, " ndiffSOMO(2)=",ndiffSOMO, " ndiffDOMO=",ndiffDOMO
     !endif
     !Jcfg = psi_configuration(:,:,i)
     !print *,"nxordiffSOMODOMO(4)=",nxordiffSOMODOMO, " ndiffSOMO(2)=",ndiffSOMO

     if((nxordiffSOMODOMO .EQ. 4) .AND. ndiffSOMO .EQ. 2) then
        select case(ndiffDOMO)
        case (0)
           ! SOMO -> VMO
           !print *,"obt SOMO -> VMO"
           extyp = 3
           !if(N_int .eq. 1) then
           !  IJsomo = IEOR(Isomo, Jsomo)
           !  p = TRAILZ(IAND(Isomo,IJsomo)) + 1
           !  IJsomo = IBCLR(IJsomo,p-1)
           !  q = TRAILZ(IJsomo) + 1
           !  !print *," p=",p," q=",q
           !  !call get_single_excitation_cfg(Jcfg, Icfg, p, q, N_int)
           !else
             ! Find p
             do ii=1,N_int
               Isomo = Ialpha(ii,1)
               Jsomo = psi_configuration(ii,1,i)
               IJsomo = IEOR(Isomo, Jsomo)
               if(popcnt(IAND(Isomo,IJsomo)) > 0)then
                 p = TRAILZ(IAND(Isomo,IJsomo)) + 1 + (ii-1) * bit_kind_size
                 EXIT
               endif
             end do
             ! Find q
             do ii=1,N_int
               Isomo = Ialpha(ii,1)
               Jsomo = psi_configuration(ii,1,i)
               IJsomo = IEOR(Isomo, Jsomo)
               iint = shiftr(p-1,bit_kind_shift) + 1
               ipos = p-shiftl((iint-1),bit_kind_shift)
               if(iint .eq. ii)then
                 IJsomo = IBCLR(IJsomo,ipos-1)
               endif
               if(popcnt(IJsomo) > 0)then
                 q = TRAILZ(IJsomo) + 1 + (ii-1) * bit_kind_size
                 EXIT
               endif
             enddo
           !endif
           !assert ( p == pp)
           !assert ( q == qq)
           !print *," 1--- p=",p," q=",q
        case (1)
           ! DOMO -> VMO
           ! or
           ! SOMO -> SOMO
           if(nsomoJ .GT. nsomoalpha) then
              ! DOMO -> VMO
              !print *,"obt DOMO -> VMO"
              extyp = 2
              !if(N_int.eq.1)then
              !  p = TRAILZ(IEOR(Idomo,Jdomo)) + 1
              !  Isomo = IEOR(Isomo, Jsomo)
              !  Isomo = IBCLR(Isomo,p-1)
              !  q = TRAILZ(Isomo) + 1
              !else

                ! Find p
                do ii=1,N_int
                  Isomo = Ialpha(ii,1)
                  Jsomo = psi_configuration(ii,1,i)
                  Idomo = Ialpha(ii,2)
                  Jdomo = psi_configuration(ii,2,i)
                  if(popcnt(IEOR(Idomo,Jdomo)) > 0)then
                    p = TRAILZ(IEOR(Idomo,Jdomo)) + 1 + (ii-1) * bit_kind_size
                    EXIT
                  endif
                end do
                ! Find q
                do ii=1,N_int
                  Isomo = Ialpha(ii,1)
                  Jsomo = psi_configuration(ii,1,i)
                  IJsomo = IEOR(Isomo, Jsomo)
                  iint = shiftr(p-1,bit_kind_shift) + 1
                  ipos = p-shiftl((iint-1),bit_kind_shift)
                  if(iint .eq. ii)then
                    IJsomo = IBCLR(IJsomo,ipos-1)
                  endif
                  if(popcnt(IJsomo) > 0)then
                    q = TRAILZ(IJsomo) + 1 + (ii-1) * bit_kind_size
                    EXIT
                  endif
                end do
              !endif
           !assert ( p == pp)
           !assert ( q == qq)
           else
              ! SOMO -> SOMO
              !print *,"obt SOMO -> SOMO"
              extyp = 1
              !if(N_int.eq.1)then
              !  q = TRAILZ(IEOR(Idomo,Jdomo)) + 1
              !  Isomo = IEOR(Isomo, Jsomo)
              !  Isomo = IBCLR(Isomo,q-1)
              !  p = TRAILZ(Isomo) + 1
              !  ! Check for Minimal alpha electrons (MS)
              !  !if(POPCNT(Isomo).lt.MS)then
              !  !  cycle
              !  !endif
              !else
                ! Find p
                !print *,"Ialpha somo=",Ialpha(1,1), Ialpha(2,1)," Ialpha domo=",Ialpha(1,2), Ialpha(2,2)
                !print *,"J somo=",psi_configuration(1,1,i), psi_configuration(2,1,i)," J domo=",psi_configuration(1,2,i),&
                !psi_configuration(2,2,i)
                do ii=1,N_int
                  Isomo = Ialpha(ii,1)
                  Jsomo = psi_configuration(ii,1,i)
                  Idomo = Ialpha(ii,2)
                  Jdomo = psi_configuration(ii,2,i)
                  if(popcnt(IEOR(Idomo,Jdomo)) > 0)then
                    q = TRAILZ(IEOR(Idomo,Jdomo)) + 1 + (ii-1) * bit_kind_size
                    EXIT
                  endif
                enddo
                ! Find q
                do ii=1,N_int
                  Isomo = Ialpha(ii,1)
                  Jsomo = psi_configuration(ii,1,i)
                  IJsomo = IEOR(Isomo, Jsomo)
                  iint = shiftr(q-1,bit_kind_shift) + 1
                  ipos = q-shiftl((iint-1),bit_kind_shift)
                  if(iint .eq. ii)then
                    IJsomo = IBCLR(IJsomo,ipos-1)
                  endif
                  !print *,"ii=",ii," Isomo=",Isomo
                  if(popcnt(IJsomo) > 0)then
                    p = TRAILZ(IJsomo) + 1 + (ii-1) * bit_kind_size
                    EXIT
                  endif
                enddo
              !endif
           !assert ( p == pp)
           !assert ( q == qq)
           endif
           !print *," 2--- p=",p," q=",q
        case (2)
           ! DOMO -> SOMO
           !print *,"obt DOMO -> SOMO"
           extyp = 4
           !if(N_int.eq.1)then
           !  IJsomo = IEOR(Isomo, Jsomo)
           !  p = TRAILZ(IAND(Jsomo,IJsomo)) + 1
           !  IJsomo = IBCLR(IJsomo,p-1)
           !  q = TRAILZ(IJsomo) + 1
           !else
             ! Find p
             do ii=1,N_int
               Isomo = Ialpha(ii,1)
               Jsomo = psi_configuration(ii,1,i)
               Idomo = Ialpha(ii,2)
               Jdomo = psi_configuration(ii,2,i)
               IJsomo = IEOR(Isomo, Jsomo)
               if(popcnt(IAND(Jsomo,IJsomo)) > 0)then
                 p = TRAILZ(IAND(Jsomo,IJsomo)) + 1 + (ii-1) * bit_kind_size
                 EXIT
               endif
             enddo
             ! Find q
             do ii=1,N_int
               Isomo = Ialpha(ii,1)
               Jsomo = psi_configuration(ii,1,i)
               IJsomo = IEOR(Isomo, Jsomo)
               iint = shiftr(p-1,bit_kind_shift) + 1
               ipos = p-shiftl((iint-1),bit_kind_shift)
               if(iint .eq. ii)then
                 IJsomo = IBCLR(IJsomo,ipos-1)
               endif
               if(popcnt(IJsomo) > 0)then
                 q = TRAILZ(IJsomo) + 1 + (ii-1) * bit_kind_size
                 EXIT
               endif
             enddo
           !endif
           !assert ( p == pp)
           !assert ( q == qq)
           !print *," 3--- p=",p," q=",q
        case default
           print *,"something went wront in get connectedI"
        end select
        starti = psi_config_data(i,1)
        endi   = psi_config_data(i,2)
        nconnectedExtradiag+=1
        nconnectedI += 1
        do ii=1,N_int
          connectedI(ii,1,nconnectedI) = psi_configuration(ii,1,i)
          connectedI(ii,2,nconnectedI) = psi_configuration(ii,2,i)
        enddo
        idxs_connectedI(nconnectedI)=starti
        excitationIds(1,nconnectedI)=p
        excitationIds(2,nconnectedI)=q
        excitationTypes(nconnectedI) = extyp
        diagfactors(nconnectedI) = 1.0d0
     else if((ndiffSOMO + ndiffDOMO) .EQ. 0) then
        ! find out all pq holes possible
        !print *,"I = ",i
        !print *,"I somo= ",psi_configuration(1,1,i), " domo=", psi_configuration(1,2,i)
        !print *,"alp somo= ",Ialpha(1,1), " domo=", Ialpha(1,2)
        nholes = 0
        ! holes in SOMO
        !Isomo = psi_configuration(1,1,i)
        !Idomo = psi_configuration(1,2,i)
        !do iii = 1,n_act_orb
        !  ii = list_act(iii)
        !   if(POPCNT(IAND(Isomo,IBSET(0_8,ii-1))) .EQ. 1) then
        !      nholes += 1
        !      listholes(nholes) = ii
        !      holetype(nholes) = 1
        !   endif
        !end do
        call bitstring_to_list(psi_configuration(1,1,i),listall,nelall,N_int)

        do iii=1,nelall
          nholes += 1
          listholes(nholes) = listall(iii)
          holetype(nholes) = 1
        end do

        ! holes in DOMO
        !do iii = 1,n_act_orb
        !  ii = list_act(iii)
        !   if(POPCNT(IAND(Idomo,IBSET(0_8,ii-1))) .EQ. 1) then
        !      nholes += 1
        !      listholes(nholes) = ii
        !      holetype(nholes) = 2
        !   endif
        !end do
        nelall=0
        listall=0
        call bitstring_to_list(psi_configuration(1,2,i),listall,nelall,N_int)

        do iii=1,nelall
          if(listall(iii) .gt. n_core_orb)then
            nholes += 1
            listholes(nholes) = listall(iii)
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
              nconnectedDiag+=1
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
              nconnectedDiag+=1
              nconnectedI += 1
              connectedI(:,:,nconnectedI) = psi_configuration(:,:,i)
              idxs_connectedI(nconnectedI)=starti
              excitationIds(1,nconnectedI)=p
              excitationIds(2,nconnectedI)=q
              excitationTypes(nconnectedI) = extyp
              diagfactors(nconnectedI) = 2.0d0
           endif
           !print *,excitationIds(1,nconnectedI), excitationIds(2,nconnectedI)
        enddo
     endif
  end do
  !print *,"nconnectedExtradiag=",nconnectedExtradiag," nconnectedDiad=",nconnectedDiag

end subroutine obtain_connected_I_foralpha
