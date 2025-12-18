program read_t2_ccsd
 implicit none
 io_amplitudes = "Read"
 touch io_amplitudes
 call routine

end

subroutine routine
 implicit none
 double precision, allocatable :: t2(:,:,:,:),t1(:,:)
 integer, allocatable :: index_t2(:,:),index_t1(:,:)
 integer, allocatable :: iorder_t2(:),iorder_t1(:)
 double precision, allocatable :: abs_t2(:), abs_t1(:)
 integer :: nO,nV,ntot_t2,ntot_t1
 integer :: i_tmp,i,j,a,b,index_large
 nO = cc_nOa
 nV = cc_nVa
 print*,'nO,nV = ',nO,nV
 ntot_t2 = nO*nO*nV*nV
 ntot_t1 = nO*nV
 ! reading T1 and T2 amplitudes 
 allocate(t2(nO, nO, nV, nV),t1(nO, nV))
 call read_t2(nO,nV,t2)
 call read_t1(nO,nV,t1)
 ! Sorting T1 and T2 amplitudes 
 allocate(index_t2(4,ntot_t2),index_t1(2,ntot_t1) )
 allocate(abs_t2(ntot_t2),abs_t1(ntot_t1))
 allocate(iorder_t2(ntot_t2),iorder_t1(ntot_t1))
 i_tmp = 0
 do b = 1, nV
   do a = 1, nV
     do j = 1, nO
       do i = 1, nO
        i_tmp += 1
        iorder_t2(i_tmp) = i_tmp
        abs_t2(i_tmp) = -dabs(t2(i,j,a,b))
!        print*,'abs_t2 = ',abs_t2(i_tmp)
        index_t2(1,i_tmp) = i
        index_t2(2,i_tmp) = j
        index_t2(3,i_tmp) = a
        index_t2(4,i_tmp) = b
       enddo
     enddo
   enddo
 enddo
 if(i_tmp .ne.ntot_t2)then
  print*,'pb !!'
  print*,'i_tmp .ne.ntot_t2',i_tmp,ntot_t2
 endif
 call dsort(abs_t2, iorder_t2, ntot_t2) 
 print*,'Two largest T2 amplitudes '
 index_large = iorder_t2(1)
 i = index_t2(1,index_large) 
 j = index_t2(2,index_large) 
 a = index_t2(3,index_large) 
 b = index_t2(4,index_large) 
 print*,'----'
 print*,'i,j,a,b',i+n_core_orb,j+n_core_orb,a,b
 print*,'T2     ',t2(i,j,a,b)
 print*,'----'
 index_large = iorder_t2(2)
 i = index_t2(1,index_large) 
 j = index_t2(2,index_large) 
 a = index_t2(3,index_large) 
 b = index_t2(4,index_large) 
 print*,'i,j,a,b',i+n_core_orb,j+n_core_orb,a,b
 print*,'T2     ',t2(i,j,a,b)
 print*,'----'
 i_tmp = 0
 do a = 1, nV
   do i = 1, nO
    i_tmp += 1
    iorder_t1(i_tmp) = i_tmp
    abs_t1(i_tmp) = -dabs(t1(i,a))
    index_t1(1,i_tmp) = i
    index_t1(2,i_tmp) = a
    enddo
  enddo
 if(i_tmp .ne.ntot_t1)then
  print*,'pb !!'
  print*,'i_tmp .ne.ntot_t1',i_tmp,ntot_t1
 endif
 call dsort(abs_t1, iorder_t1, ntot_t1) 
 print*,'Two largest T1 amplitudes'
 index_large = iorder_t1(1)
 i = index_t1(1,index_large) 
 a = index_t1(2,index_large) 
 print*,'----'
 print*,'i,a',i+n_core_orb,a
 print*,'t1     ',t1(i,a)
 print*,'----'
 index_large = iorder_t1(2)
 i = index_t1(1,index_large) 
 a = index_t1(2,index_large) 
 print*,'----'
 print*,'i,a',i+n_core_orb,a
 print*,'t1     ',t1(i,a)
 print*,'----'


end
