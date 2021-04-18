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
  integer(bit_kind),intent(in)         :: Ialpha(N_int,2)
  integer(bit_kind),intent(in)         :: Jcfg(N_int,2)
  integer,intent(in)                   :: p,q
  integer,intent(in)                   :: extype
  integer,intent(out)                  :: pmodel,qmodel
  integer*8                            :: Isomo
  integer*8                            :: Idomo
  integer*8                            :: Jsomo
  integer*8                            :: Jdomo
  integer*8                            :: mask
  integer*8                            :: Isomotmp
  integer*8                            :: Jsomotmp
  integer                              :: pos0,pos0prev

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
     !print *,"input pq=",p,q,"extype=",extype
     !call debug_spindet(Isomo,1)
     !call debug_spindet(Idomo,1)
     !call debug_spindet(Jsomo,1)
     !call debug_spindet(Jdomo,1)
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
          mask = ISHFT(1_8,p) - 1
          Isomo = IAND(Isomo,mask)
          pmodel = POPCNT(mask) - POPCNT(XOR(Isomo,mask))
          mask = ISHFT(1_8,q) - 1
          Jsomo = IAND(Jsomo,mask)
          qmodel = POPCNT(mask) - POPCNT(XOR(Jsomo,mask))
       case (4)
          ! DOMO -> SOMO
          ! remove all domos except one at p
          !print *,"type -> DOMO -> SOMO"
          !Isomo = IEOR(Isomo,Jsomo)
          mask = ISHFT(1_8,p) - 1
          Jsomo = IAND(Jsomo,mask)
          pmodel = POPCNT(mask) - POPCNT(XOR(Jsomo,mask))
          mask = ISHFT(1_8,q) - 1
          Isomo = IAND(Isomo,mask)
          qmodel = POPCNT(mask) - POPCNT(XOR(Isomo,mask))
       case default
          print *,"something is wrong in convertOrbIdsToModelSpaceIds"
     end select
  endif
     !print *,p,q,"model ids=",pmodel,qmodel
end subroutine convertOrbIdsToModelSpaceIds
