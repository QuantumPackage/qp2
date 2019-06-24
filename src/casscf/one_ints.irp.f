! -*- F90 -*-
BEGIN_PROVIDER [real*8, one_ints, (mo_num,mo_num)]
      implicit none 
      integer :: i,j,kk
      logical :: lread
      real*8 :: rdum
      do i=1,mo_num
       do j=1,mo_num
        one_ints(i,j)=0.D0
       end do 
      end do
      open(unit=12,file='onetrf.tmp',status='old',form='formatted')
      lread=.true.
      do while (lread)
       read(12,*,iostat=kk) i,j,rdum
       if (kk.ne.0) then
        lread=.false.
       else
        one_ints(i,j)=rdum
        one_ints(j,i)=rdum
       end if
      end do
      close(12)
      write(6,*) ' read MCSCF natural one-electron integrals '
END_PROVIDER

