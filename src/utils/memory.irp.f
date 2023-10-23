BEGIN_PROVIDER [ integer, qp_max_mem ]
 implicit none
 BEGIN_DOC
 ! Maximum memory in Gb
 END_DOC
 character*(128) :: env
 integer, external :: get_total_available_memory

 qp_max_mem = get_total_available_memory()
 call write_int(6,qp_max_mem,'Total available memory (GB)')
 call getenv('QP_MAXMEM',env)
 if (trim(env) /= '') then
    call lock_io()
    read(env,*) qp_max_mem
    call unlock_io()
 endif
 call write_int(6,qp_max_mem,'Target maximum memory (GB)')

END_PROVIDER

subroutine resident_memory(value)
  use c_functions
  implicit none
  BEGIN_DOC
! Returns the current used memory in gigabytes used by the current process.
  END_DOC
  integer :: iunit
  integer, external :: getUnitAndOpen
  character*(32) :: key
  double precision, intent(out) :: value

  call lock_io()
  call usleep(10)

  value = 0.d0
IRP_IF MACOS
IRP_ELSE
  iunit = getUnitAndOpen('/proc/self/status','r')
  do
    read(iunit,*,err=10,end=20) key, value
    if (trim(key) == 'VmRSS:') then
      exit
    endif
    10 continue
  end do
  20 continue
  close(iunit)
IRP_ENDIF
  value = value / (1024.d0*1024.d0)
  call unlock_io()
end function

subroutine total_memory(value)
  implicit none
  BEGIN_DOC
! Returns the current used memory in gigabytes used by the current process.
  END_DOC
  integer :: iunit
  integer, external :: getUnitAndOpen
  character*(32) :: key
  double precision, intent(out) :: value

  call lock_io()
  value = 0.d0
IRP_IF MACOS
IRP_ELSE
  iunit = getUnitAndOpen('/proc/self/status','r')
  do
    read(iunit,*,err=10,end=20) key, value
    if (trim(key) == 'VmSize:') then
      exit
    endif
    10 continue
  end do
  20 continue
  close(iunit)
IRP_ENDIF
  value = value / (1024.d0*1024.d0)
  call unlock_io()
end function

double precision function memory_of_double(n)
  implicit none
  BEGIN_DOC
! Computes the memory required for n double precision elements in gigabytes.
  END_DOC
  integer, intent(in) :: n
  double precision, parameter :: f = 8.d0 / (1024.d0*1024.d0*1024.d0)
  memory_of_double = dble(n) * f
end function

double precision function memory_of_int(n)
  implicit none
  BEGIN_DOC
! Computes the memory required for n double precision elements in gigabytes.
  END_DOC
  integer, intent(in) :: n
  double precision, parameter :: f = 4.d0 / (1024.d0*1024.d0*1024.d0)
  memory_of_int = dble(n) * f
end function

subroutine check_mem(rss_in,routine)
  implicit none
  BEGIN_DOC
! Checks if n gigabytes can be allocated. If not, exit the run.
  END_DOC
  double precision, intent(in) :: rss_in
  character*(*) :: routine
  double precision :: mem
  call total_memory(mem)
  mem += rss_in
  if (mem > qp_max_mem) then
    call print_memory_usage()
    print *,  'Not enough memory: aborting in ', routine
    print *,  mem, ' GB required'
    stop -1
  endif
end

subroutine print_memory_usage()
  implicit none
  BEGIN_DOC
! Prints the memory usage in the output
  END_DOC
  double precision :: rss, mem
  call resident_memory(rss)
  call total_memory(mem)

  write(*,'(A,F14.3,A,F14.3,A)') &
    '.. >>>>> [ RES  MEM : ', rss , &
        ' GB ] [ VIRT MEM : ', mem, ' GB ] <<<<< ..'
end

integer function get_total_available_memory() result(res)
  implicit none
  BEGIN_DOC
! Returns the total available memory on the current machine
  END_DOC

  character(len=128) :: line
  integer :: status
  integer :: iunit
  integer*8, parameter :: KB = 1024
  integer*8, parameter :: GiB = 1024**3
  integer, external :: getUnitAndOpen

  iunit = getUnitAndOpen('/proc/meminfo','r')

  res = 512
  do
      read(iunit, '(A)', END=10) line
      if (line(1:10) == "MemTotal: ") then
        read(line(11:), *, ERR=20) res
        res = int((res*KB) / GiB,4)
        exit
 20     continue
      end if
  end do
 10 continue
  close(iunit)

end function get_total_available_memory


