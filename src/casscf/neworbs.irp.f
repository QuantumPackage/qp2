! -*- F90 -*-
BEGIN_PROVIDER [real*8, SXmatrix, (nMonoEx+1,nMonoEx+1)]
      implicit none
      integer :: i,j
      do i=1,nMonoEx+1
       do j=1,nMonoEx+1
        SXmatrix(i,j)=0.D0
       end do
      end do

      do i=1,nMonoEx
       SXmatrix(1,i+1)=gradvec2(i)
       SXmatrix(1+i,1)=gradvec2(i)
      end do

      do i=1,nMonoEx
       do j=1,nMonoEx
        SXmatrix(i+1,j+1)=hessmat2(i,j)
        SXmatrix(j+1,i+1)=hessmat2(i,j)
       end do
      end do

   if (bavard) then
      do i=2,nMonoEx+1
       write(6,*) ' diagonal of the Hessian : ',i,hessmat2(i,i)
      end do
   end if
      

END_PROVIDER

 BEGIN_PROVIDER [real*8, SXeigenvec, (nMonoEx+1,nMonoEx+1)]
&BEGIN_PROVIDER [real*8, SXeigenval, (nMonoEx+1)]
 END_PROVIDER

 BEGIN_PROVIDER [real*8, SXvector, (nMonoEx+1)]
&BEGIN_PROVIDER [real*8, energy_improvement]
      implicit none
      integer :: ierr,matz,i
      real*8 :: c0

      call lapack_diag(SXeigenval,SXeigenvec,SXmatrix,nMonoEx+1,nMonoEx+1)
      write(6,*) ' SXdiag : lowest 5 eigenvalues '
      write(6,*) ' 1 - ',SXeigenval(1),SXeigenvec(1,1)
      write(6,*) ' 2 - ',SXeigenval(2),SXeigenvec(1,2)
      write(6,*) ' 3 - ',SXeigenval(3),SXeigenvec(1,3)
      write(6,*) ' 4 - ',SXeigenval(4),SXeigenvec(1,4)
      write(6,*) ' 5 - ',SXeigenval(5),SXeigenvec(1,5)
      write(6,*) 
      write(6,*) ' SXdiag : lowest eigenvalue = ',SXeigenval(1)
      energy_improvement = SXeigenval(1)

integer :: best_vector
real*8 :: best_overlap
      best_overlap=0.D0
      do i=1,nMonoEx+1
       if (SXeigenval(i).lt.0.D0) then
        if (abs(SXeigenvec(1,i)).gt.best_overlap) then
         best_overlap=abs(SXeigenvec(1,i))
         best_vector=i
        end if
       end if
      end do

      write(6,*) ' SXdiag : eigenvalue for best overlap with ' 
      write(6,*) '  previous orbitals = ',SXeigenval(best_vector)
      energy_improvement = SXeigenval(best_vector)
      
      c0=SXeigenvec(1,best_vector)
      write(6,*) ' weight of the 1st element ',c0
      do i=1,nMonoEx+1
       SXvector(i)=SXeigenvec(i,best_vector)/c0
!      write(6,*) ' component No ',i,' : ',SXvector(i)
      end do

END_PROVIDER


BEGIN_PROVIDER [real*8, NewOrbs, (ao_num,mo_num) ]
      implicit none
      integer :: i,j,ialph

! form the exponential of the Orbital rotations
      call get_orbrotmat
! form the new orbitals
      do i=1,ao_num
       do j=1,mo_num
        NewOrbs(i,j)=0.D0
       end do
      end do

      do ialph=1,ao_num
       do i=1,mo_num
        wrkline(i)=mo_coef(ialph,i)
       end do
       do i=1,mo_num
        do j=1,mo_num
         NewOrbs(ialph,i)+=Umat(i,j)*wrkline(j)
        end do
       end do
      end do

END_PROVIDER

 BEGIN_PROVIDER [real*8, Tpotmat, (mo_num,mo_num) ]
&BEGIN_PROVIDER [real*8, Umat, (mo_num,mo_num) ]
&BEGIN_PROVIDER [real*8, wrkline, (mo_num) ]
&BEGIN_PROVIDER [real*8, Tmat, (mo_num,mo_num) ]
END_PROVIDER

      subroutine get_orbrotmat
      implicit none
      integer :: i,j,indx,k,iter,t,a,ii,tt,aa
      real*8 :: sum
      logical :: converged


! the orbital rotation matrix T
      do i=1,mo_num
       do j=1,mo_num
        Tmat(i,j)=0.D0
        Umat(i,j)=0.D0
        Tpotmat(i,j)=0.D0
       end do
       Tpotmat(i,i)=1.D0
      end do

      indx=1
      do i=1,n_core_orb
       ii=list_core(i)
       do t=1,n_act_orb
        tt=list_act(t)
        indx+=1
        Tmat(ii,tt)= SXvector(indx)
        Tmat(tt,ii)=-SXvector(indx)
       end do
      end do
      do i=1,n_core_orb
       ii=list_core(i)
       do a=1,n_virt_orb
        aa=list_virt(a)
        indx+=1
        Tmat(ii,aa)= SXvector(indx)
        Tmat(aa,ii)=-SXvector(indx)
       end do
      end do
      do t=1,n_act_orb
       tt=list_act(t)
       do a=1,n_virt_orb
        aa=list_virt(a)
        indx+=1
        Tmat(tt,aa)= SXvector(indx)
        Tmat(aa,tt)=-SXvector(indx)
       end do
      end do

      write(6,*) ' the T matrix '
      do indx=1,nMonoEx
       i=excit(1,indx)
       j=excit(2,indx)
!       if (abs(Tmat(i,j)).gt.1.D0) then
!        write(6,*) ' setting matrix element ',i,j,' of ',Tmat(i,j),' to ' &
!             , sign(1.D0,Tmat(i,j))
!        Tmat(i,j)=sign(1.D0,Tmat(i,j))
!        Tmat(j,i)=-Tmat(i,j)
!       end if
        if (abs(Tmat(i,j)).gt.1.D-9) write(6,9901) i,j,excit_class(indx),Tmat(i,j)
  9901 format('   ',i4,' -> ',i4,' (',A3,') : ',E14.6)
      end do

      write(6,*) 
      write(6,*) ' forming the matrix exponential '
      write(6,*) 
! form the exponential
      iter=0
      converged=.false.
      do while (.not.converged)
       iter+=1
! add the next term
       do i=1,mo_num
        do j=1,mo_num
         Umat(i,j)+=Tpotmat(i,j)
        end do
       end do
! next power of T, we multiply Tpotmat with Tmat/iter
       do i=1,mo_num
        do j=1,mo_num
         wrkline(j)=Tpotmat(i,j)/dble(iter)
         Tpotmat(i,j)=0.D0
        end do
        do j=1,mo_num
         do k=1,mo_num
          Tpotmat(i,j)+=wrkline(k)*Tmat(k,j)
         end do
        end do
       end do
! Convergence test
       sum=0.D0
       do i=1,mo_num
        do j=1,mo_num
         sum+=abs(Tpotmat(i,j))
        end do
       end do
       write(6,*) ' Iteration No ',iter,' Sum = ',sum
       if (sum.lt.1.D-6) then
        converged=.true.
       end if
       if (iter.ge.NItExpMax) then
        stop ' no convergence '
       end if
      end do
      write(6,*)
      write(6,*) ' Converged ! '
      write(6,*)

      end subroutine get_orbrotmat

BEGIN_PROVIDER [integer, NItExpMax]
   NItExpMax=100
END_PROVIDER


