! -*- F90 -*- 
 BEGIN_PROVIDER [real*8, Fipq, (mo_num,mo_num) ]
&BEGIN_PROVIDER [real*8, Fapq, (mo_num,mo_num) ]
BEGIN_DOC
! the inactive and the active Fock matrices, in molecular 
! orbitals
! we create them in MOs, quite expensive
!
! for an implementation in AOs we need first the natural orbitals 
! for forming an active density matrix in AOs
!
END_DOC
       implicit none
       double precision, allocatable :: integrals_array1(:,:)
       double precision, allocatable :: integrals_array2(:,:)
       integer :: p,q,k,kk,t,tt,u,uu
       allocate(integrals_array1(mo_num,mo_num))
       allocate(integrals_array2(mo_num,mo_num))

       do p=1,mo_num
        do q=1,mo_num
         Fipq(p,q)=one_ints(p,q)
         Fapq(p,q)=0.D0
        end do
       end do

! the inactive Fock matrix
       do k=1,n_core_orb
        kk=list_core(k)
        do p=1,mo_num
         do q=1,mo_num
          Fipq(p,q)+=2.D0*bielec_pqxx(p,q,k,k) -bielec_pxxq(p,k,k,q)
         end do
        end do
       end do

! the active Fock matrix, D0tu is diagonal
       do t=1,n_act_orb
        tt=list_act(t)
        do p=1,mo_num
         do q=1,mo_num
          Fapq(p,q)+=occnum(tt)  &
              *(bielec_pqxx(p,q,tt,tt)-0.5D0*bielec_pxxq(p,tt,tt,q))
         end do
        end do
       end do

if (bavard) then
integer :: i
         write(6,*)
         write(6,*) ' the effective Fock matrix over MOs'
         write(6,*)

         write(6,*)
         write(6,*) ' the diagonal of the inactive effective Fock matrix '
         write(6,'(5(i3,F12.5))') (i,Fipq(i,i),i=1,mo_num)
         write(6,*)
         write(6,*)
         write(6,*) ' the diagonal of the active Fock matrix '
         write(6,'(5(i3,F12.5))') (i,Fapq(i,i),i=1,mo_num)
         write(6,*)
end if


END_PROVIDER


