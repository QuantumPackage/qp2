BEGIN_PROVIDER [real*8, Fipq, (mo_num,mo_num) ]
   BEGIN_DOC
   ! the inactive Fock matrix, in molecular orbitals
   END_DOC
   implicit none
   integer                        :: p,q,k,kk,t,tt,u,uu
   
   do q=1,mo_num
     do p=1,mo_num
       Fipq(p,q)=one_ints_no(p,q)
     end do
   end do
   
   ! the inactive Fock matrix
   do k=1,n_core_orb
     kk=list_core(k)
     do q=1,mo_num
       do p=1,mo_num
         Fipq(p,q)+=2.D0*bielec_pqxx_no(p,q,k,k) -bielec_pxxq_no(p,k,k,q)
       end do
     end do
   end do
   
   if (bavard) then
     integer                        :: i
     write(6,*)
     write(6,*) ' the diagonal of the inactive effective Fock matrix '
     write(6,'(5(i3,F12.5))') (i,Fipq(i,i),i=1,mo_num)
     write(6,*)
   end if
   
   
END_PROVIDER
 
 
BEGIN_PROVIDER [real*8, Fapq, (mo_num,mo_num) ]
   BEGIN_DOC
   ! the active active Fock matrix, in molecular orbitals
   ! we create them in MOs, quite expensive
   !
   ! for an implementation in AOs we need first the natural orbitals
   ! for forming an active density matrix in AOs
   !
   END_DOC
   implicit none
   integer                        :: p,q,k,kk,t,tt,u,uu
   
   Fapq = 0.d0
   
   ! the active Fock matrix, D0tu is diagonal
   do t=1,n_act_orb
     tt=list_act(t)
     do q=1,mo_num
       do p=1,mo_num
         Fapq(p,q)+=occnum(tt)                                       &
             *(bielec_pqxx_no(p,q,tt,tt)-0.5D0*bielec_pxxq_no(p,tt,tt,q))
       end do
     end do
   end do
   
   if (bavard) then
     integer                        :: i
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
 
 
