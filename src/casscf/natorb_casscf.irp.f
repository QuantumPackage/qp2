! -*- F90 -*-
BEGIN_PROVIDER [real*8, occnum, (mo_num)]
    implicit none
    integer :: i,kk,j
    logical :: lread
    real*8 :: rdum
    do i=1,mo_num
     occnum(i)=0.D0
    end do
    do i=1,n_core_orb
     occnum(list_core(i))=2.D0
    end do

    open(unit=12,file='D0tu.dat',form='formatted',status='old')
    lread=.true.
    do while (lread)
     read(12,*,iostat=kk) i,j,rdum
     if (kk.ne.0) then
      lread=.false.
     else
      if (i.eq.j) then
       occnum(list_act(i))=rdum
      else
       write(6,*) ' WARNING - no natural orbitals !'
       write(6,*) i,j,rdum
      end if
     end if
    end do
    close(12)
    write(6,*) ' read occupation numbers '
    do i=1,mo_num
     write(6,*) i,occnum(i)
    end do

END_PROVIDER

BEGIN_PROVIDER [real*8, P0tuvx_no, (n_act_orb,n_act_orb,n_act_orb,n_act_orb)]
      implicit none
      integer :: i,j,k,l,kk
      real*8 :: rdum
      logical :: lread

      do i=1,n_act_orb
       do j=1,n_act_orb
        do k=1,n_act_orb
         do l=1,n_act_orb
          P0tuvx_no(l,k,j,i)=0.D0
         end do
        end do
       end do
      end do

      open(unit=12,file='P0tuvx.dat',form='formatted',status='old')
      lread=.true.
      do while (lread)
       read(12,*,iostat=kk) i,j,k,l,rdum
       if (kk.ne.0) then
        lread=.false.
       else
        P0tuvx_no(i,j,k,l)=rdum
       end if
      end do
      close(12)
      write(6,*) ' read the 2-particle density matrix '
END_PROVIDER
