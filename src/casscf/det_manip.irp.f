use bitmasks

subroutine do_signed_mono_excitation(key1,key2,nu,ihole,ipart,       &
      ispin,phase,ierr)
  BEGIN_DOC
  ! we create the mono-excitation, and determine, if possible,
  ! the phase and the number in the list of determinants
  END_DOC
  implicit none
  integer(bit_kind)              :: key1(N_int,2),key2(N_int,2)
  integer(bit_kind), allocatable :: keytmp(:,:)
  integer                        :: exc(0:2,2,2),ihole,ipart,ierr,nu,ispin
  real*8                         :: phase
  logical                        :: found
  allocate(keytmp(N_int,2))
  
  nu=-1
  phase=1.D0
  ierr=0
  call det_copy(key1,key2,N_int)
  !        write(6,*) ' key2 before excitation ',ihole,' -> ',ipart,' spin = ',ispin
  !        call print_det(key2,N_int)
  call do_single_excitation(key2,ihole,ipart,ispin,ierr)
  !        write(6,*) ' key2 after ',ihole,' -> ',ipart,' spin = ',ispin
  !        call print_det(key2,N_int)
  !        write(6,*) ' excitation ',ihole,' -> ',ipart,' gives ierr = ',ierr
  if (ierr.eq.1) then
    ! excitation is possible
    ! get the phase
    call get_single_excitation(key1,key2,exc,phase,N_int)
    ! get the number in the list
    found=.false.
    nu=0

    !TODO BOTTLENECK
    do while (.not.found)
      nu+=1
      if (nu.gt.N_det) then
        ! the determinant is possible, but not in the list
        found=.true.
        nu=-1
      else
        call det_extract(keytmp,nu,N_int)
        integer                        :: i,ii
        found=.true.
        do ii=1,2
          do i=1,N_int
            if (keytmp(i,ii).ne.key2(i,ii)) then
              found=.false.
            end if
          end do
        end do
      end if
    end do
  end if
  !
  ! we found the new string, the phase, and possibly the number in the list
  !
end subroutine do_signed_mono_excitation

subroutine det_extract(key,nu,Nint)
  BEGIN_DOC
  ! extract a determinant from the list of determinants
  END_DOC
  implicit none
  integer                        :: ispin,i,nu,Nint
  integer(bit_kind)              :: key(Nint,2)
  do ispin=1,2
    do i=1,Nint
      key(i,ispin)=psi_det(i,ispin,nu)
    end do
  end do
end subroutine det_extract

subroutine det_copy(key1,key2,Nint)
  use bitmasks ! you need to include the bitmasks_module.f90 features
  BEGIN_DOC
  ! copy a determinant from key1 to key2
  END_DOC
  implicit none
  integer                        :: ispin,i,Nint
  integer(bit_kind)              :: key1(Nint,2),key2(Nint,2)
  do ispin=1,2
    do i=1,Nint
      key2(i,ispin)=key1(i,ispin)
    end do
  end do
end subroutine det_copy

subroutine do_spinfree_mono_excitation(key_in,key_out1,key_out2      &
      ,nu1,nu2,ihole,ipart,phase1,phase2,ierr,jerr)
  BEGIN_DOC
  ! we create the spin-free mono-excitation E_pq=(a^+_p a_q + a^+_P a_Q)
  ! we may create two determinants as result
  !
  END_DOC
  implicit none
  integer(bit_kind)              :: key_in(N_int,2),key_out1(N_int,2)
  integer(bit_kind)              :: key_out2(N_int,2)
  integer                        :: ihole,ipart,ierr,jerr,nu1,nu2
  integer                        :: ispin
  real*8                         :: phase1,phase2
  
  !        write(6,*) ' applying E_',ipart,ihole,' on determinant '
  !        call print_det(key_in,N_int)
  
  ! spin alpha
  ispin=1
  call do_signed_mono_excitation(key_in,key_out1,nu1,ihole           &
      ,ipart,ispin,phase1,ierr)
  !        if (ierr.eq.1) then
  !         write(6,*) ' 1 result is ',nu1,phase1
  !         call print_det(key_out1,N_int)
  !        end if
  ! spin beta
  ispin=2
  call do_signed_mono_excitation(key_in,key_out2,nu2,ihole           &
      ,ipart,ispin,phase2,jerr)
  !        if (jerr.eq.1) then
  !         write(6,*) ' 2 result is ',nu2,phase2
  !         call print_det(key_out2,N_int)
  !        end if
  
end subroutine do_spinfree_mono_excitation

