BEGIN_PROVIDER [real*8, SXmatrix, (nMonoEx+1,nMonoEx+1)]
  implicit none
  BEGIN_DOC
  ! Single-excitation matrix
  END_DOC
  
  integer                        :: i,j
  
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
  implicit none
  BEGIN_DOC
  ! Eigenvectors/eigenvalues of the single-excitation matrix
  END_DOC
  call lapack_diag(SXeigenval,SXeigenvec,SXmatrix,nMonoEx+1,nMonoEx+1)
END_PROVIDER

 BEGIN_PROVIDER [real*8, SXvector, (nMonoEx+1)]
&BEGIN_PROVIDER [real*8, energy_improvement]
  implicit none
  BEGIN_DOC
  ! Best eigenvector of the single-excitation matrix
  END_DOC
  integer                        :: ierr,matz,i
  real*8                         :: c0
  
  write(6,*) ' SXdiag : lowest 5 eigenvalues '
  write(6,*) ' 1 - ',SXeigenval(1),SXeigenvec(1,1)
  write(6,*) ' 2 - ',SXeigenval(2),SXeigenvec(1,2)
  write(6,*) ' 3 - ',SXeigenval(3),SXeigenvec(1,3)
  write(6,*) ' 4 - ',SXeigenval(4),SXeigenvec(1,4)
  write(6,*) ' 5 - ',SXeigenval(5),SXeigenvec(1,5)
  write(6,*)
  write(6,*) ' SXdiag : lowest eigenvalue = ',SXeigenval(1)
  energy_improvement = SXeigenval(1)
  
  integer                        :: best_vector
  real*8                         :: best_overlap
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
  BEGIN_DOC
  ! Updated orbitals
  END_DOC
  integer                        :: i,j,ialph
  
  call dgemm('N','T', ao_num,mo_num,mo_num,1.d0,                     &
      NatOrbsFCI, size(NatOrbsFCI,1),                                &
      Umat, size(Umat,1), 0.d0,                                      &
      NewOrbs, size(NewOrbs,1))
  
END_PROVIDER

BEGIN_PROVIDER [real*8, Umat, (mo_num,mo_num) ]
  implicit none
  BEGIN_DOC
  ! Orbital rotation matrix
  END_DOC
  integer                        :: i,j,indx,k,iter,t,a,ii,tt,aa
  logical                        :: converged
  
  real*8 :: Tpotmat (mo_num,mo_num), Tpotmat2 (mo_num,mo_num) 
  real*8 :: Tmat(mo_num,mo_num) 
  real*8 :: f
  
  ! the orbital rotation matrix T
  Tmat(:,:)=0.D0
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
  
  ! Form the exponential

  Tpotmat(:,:)=0.D0
  Umat(:,:)   =0.D0
  do i=1,mo_num
    Tpotmat(i,i)=1.D0
    Umat(i,i)   =1.d0
  end do
  iter=0
  converged=.false.
  do while (.not.converged)
    iter+=1
    f = 1.d0 / dble(iter)
    Tpotmat2(:,:) = Tpotmat(:,:) * f
    call dgemm('N','N', mo_num,mo_num,mo_num,1.d0,                   &
        Tpotmat2, size(Tpotmat2,1),                                  &
        Tmat, size(Tmat,1), 0.d0,                                    &
        Tpotmat, size(Tpotmat,1))
    Umat(:,:) = Umat(:,:) + Tpotmat(:,:)
    
    converged = ( sum(abs(Tpotmat(:,:))) < 1.d-6).or.(iter>30)
  end do
END_PROVIDER
  


