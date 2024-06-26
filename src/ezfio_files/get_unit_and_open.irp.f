
integer function getUnitAndOpen(f,mode)
  implicit none

  BEGIN_DOC
!  :f:
!     file name
!
!  :mode:
!     'R' : READ, UNFORMATTED
!     'W' : WRITE, UNFORMATTED
!     'r' : READ, FORMATTED
!     'w' : WRITE, FORMATTED
!     'a' : APPEND, FORMATTED
!     'x' : READ/WRITE, FORMATTED
!
  END_DOC

  character*(*)      :: f
  character*(256)    :: new_f
  integer            :: iunit
  logical            :: is_open, exists
  character          :: mode

  is_open = .True.
  iunit = 20
  new_f = f
  do while (is_open)
    inquire(unit=iunit,opened=is_open)
    if (.not.is_open) then
      getUnitAndOpen = iunit
    endif
    iunit = iunit+1
  enddo
  if (mode.eq.'r') then
    inquire(file=f,exist=exists)
    if (.not.exists) then
      open(unit=getUnitAndOpen,file=f,status='NEW',action='WRITE',form='FORMATTED')
      close(unit=getUnitAndOpen)
    endif
    open(unit=getUnitAndOpen,file=f,status='OLD',action='READ',form='FORMATTED')
  else if (mode.eq.'R') then
    inquire(file=f,exist=exists)
    if (.not.exists) then
      open(unit=getUnitAndOpen,file=f,status='NEW',action='WRITE',form='UNFORMATTED')
      close(unit=getUnitAndOpen)
    endif
    open(unit=getUnitAndOpen,file=f,status='OLD',action='READ',form='UNFORMATTED')
  else if (mode.eq.'W') then
    open(unit=getUnitAndOpen,file=new_f,status='UNKNOWN',action='READWRITE',form='UNFORMATTED')
  else if (mode.eq.'A') then
    open(unit=getUnitAndOpen,file=new_f,status='UNKNOWN',action='READWRITE',position='APPEND',form='UNFORMATTED')
  else if (mode.eq.'w') then
    open(unit=getUnitAndOpen,file=new_f,status='UNKNOWN',action='READWRITE',form='FORMATTED')
  else if (mode.eq.'a') then
    open(unit=getUnitAndOpen,file=new_f,status='UNKNOWN',action='READWRITE',position='APPEND',form='FORMATTED')
  else if (mode.eq.'x') then
    open(unit=getUnitAndOpen,file=new_f,form='FORMATTED')
  endif
end function getUnitAndOpen

