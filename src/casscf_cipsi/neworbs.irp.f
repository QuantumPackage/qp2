 BEGIN_PROVIDER [real*8, SXmatrix, (nMonoEx+1,nMonoEx+1)]
&BEGIN_PROVIDER [integer, n_guess_sx_mat ]
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
  if(diag_hess_cas)then  
   do i = 1, nMonoEx
    SXmatrix(i+1,i+1) = hessdiag(i)
   enddo
  else
   do i=1,nMonoEx
     do j=1,nMonoEx
       SXmatrix(i+1,j+1)=hessmat(i,j)
       SXmatrix(j+1,i+1)=hessmat(i,j)
     end do
   end do
  endif
  
  do i = 1, nMonoEx
   SXmatrix(i+1,i+1) += level_shift_casscf
  enddo
  n_guess_sx_mat = 1
  do i = 1, nMonoEx
   if(SXmatrix(i+1,i+1).lt.0.d0 )then
    n_guess_sx_mat += 1
   endif
  enddo
  if (bavard) then
    do i=2,nMonoEx
      write(6,*) ' diagonal of the Hessian : ',i,hessmat(i,i)
    end do
  end if
  
END_PROVIDER

 BEGIN_PROVIDER [real*8, SXeigenvec, (nMonoEx+1,nMonoEx+1)]
&BEGIN_PROVIDER [real*8, SXeigenval, (nMonoEx+1)]
  implicit none
  BEGIN_DOC
  ! Eigenvectors/eigenvalues of the single-excitation matrix
  END_DOC
  if(nMonoEx+1.gt.n_det_max_full)then
   if(bavard)then
    print*,'Using the Davidson algorithm to diagonalize the SXmatrix'
   endif
   double precision, allocatable :: u_in(:,:),energies(:)
   allocate(u_in(nMonoEx+1,n_states_diag),energies(n_guess_sx_mat))
   call davidson_diag_sx_mat(n_guess_sx_mat, u_in, energies)
   integer :: i,j
   SXeigenvec = 0.d0
   SXeigenval = 0.d0
   do i = 1, n_guess_sx_mat
    SXeigenval(i) = energies(i)
    do j = 1, nMonoEx+1
     SXeigenvec(j,i) = u_in(j,i)
    enddo
   enddo
  else
   if(bavard)then
    print*,'Diagonalize the SXmatrix with Jacobi'
   endif
   call lapack_diag(SXeigenval,SXeigenvec,SXmatrix,nMonoEx+1,nMonoEx+1)
  endif
  if (bavard) then
    write(6,*) ' SXdiag : lowest eigenvalues '
    write(6,*) ' 1 - ',SXeigenval(1),SXeigenvec(1,1)
    if(n_guess_sx_mat.gt.0)then
     write(6,*) ' 2 - ',SXeigenval(2),SXeigenvec(1,2)
     write(6,*) ' 3 - ',SXeigenval(3),SXeigenvec(1,3)
     write(6,*) ' 4 - ',SXeigenval(4),SXeigenvec(1,4)
     write(6,*) ' 5 - ',SXeigenval(5),SXeigenvec(1,5)
    endif
    write(6,*)
    write(6,*) ' SXdiag : lowest eigenvalue = ',SXeigenval(1)
  endif
END_PROVIDER

 BEGIN_PROVIDER [real*8, energy_improvement]
 implicit none
 if(state_following_casscf)then
  energy_improvement = SXeigenval(best_vector_ovrlp_casscf)
 else 
  energy_improvement = SXeigenval(1)
 endif
 END_PROVIDER 



 BEGIN_PROVIDER [ integer, best_vector_ovrlp_casscf ]
&BEGIN_PROVIDER [ double precision,  best_overlap_casscf ]
  implicit none
  integer :: i
  double precision :: c0
  best_overlap_casscf = 0.D0
  best_vector_ovrlp_casscf = -1000
  do i=1,nMonoEx+1
    if (SXeigenval(i).lt.0.D0) then
      if (dabs(SXeigenvec(1,i)).gt.best_overlap_casscf) then
        best_overlap_casscf=dabs(SXeigenvec(1,i))
        best_vector_ovrlp_casscf = i
      end if
    end if
  end do
  if(best_vector_ovrlp_casscf.lt.0)then 
   best_vector_ovrlp_casscf = minloc(SXeigenval,nMonoEx+1) 
  endif
  c0=SXeigenvec(1,best_vector_ovrlp_casscf)
  if (bavard) then
    write(6,*) ' SXdiag : eigenvalue for best overlap with '
    write(6,*) ' previous orbitals = ',SXeigenval(best_vector_ovrlp_casscf)
    write(6,*) ' weight of the 1st element ',c0
  endif
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, SXvector, (nMonoEx+1)]
  implicit none
  BEGIN_DOC
  ! Best eigenvector of the single-excitation matrix
  END_DOC
  integer           :: i
  double precision  :: c0
  c0=SXeigenvec(1,best_vector_ovrlp_casscf)
  do i=1,nMonoEx+1
    SXvector(i)=SXeigenvec(i,best_vector_ovrlp_casscf)/c0
  end do
 END_PROVIDER


BEGIN_PROVIDER [double precision, NewOrbs, (ao_num,mo_num) ]
  implicit none
  BEGIN_DOC
  ! Updated orbitals
  END_DOC
  integer                        :: i,j,ialph
  
  if(state_following_casscf)then
   print*,'Using the state following casscf '
   call dgemm('N','T', ao_num,mo_num,mo_num,1.d0,                     &
       NatOrbsFCI, size(NatOrbsFCI,1),                                &
       Umat, size(Umat,1), 0.d0,                                      &
       NewOrbs, size(NewOrbs,1))

    level_shift_casscf *= 0.5D0
    level_shift_casscf = max(level_shift_casscf,0.002d0)
   !touch level_shift_casscf
  else
   if(best_vector_ovrlp_casscf.ne.1.and.n_orb_swap.ne.0)then
     print*,'Taking the lowest root for the CASSCF'
     print*,'!!! SWAPPING MOS !!!!!!'
     level_shift_casscf *= 2.D0
    level_shift_casscf = min(level_shift_casscf,0.5d0)
     print*,'level_shift_casscf = ',level_shift_casscf
     NewOrbs = switch_mo_coef
    !mo_coef = switch_mo_coef
    !soft_touch mo_coef
    !call save_mos_no_occ
    !stop
   else 
    level_shift_casscf *= 0.5D0
    level_shift_casscf = max(level_shift_casscf,0.002d0)
   !touch level_shift_casscf
    call dgemm('N','T', ao_num,mo_num,mo_num,1.d0,                     &
        NatOrbsFCI, size(NatOrbsFCI,1),                                &
        Umat, size(Umat,1), 0.d0,                                      &
        NewOrbs, size(NewOrbs,1))
   endif
  endif
  
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
  do i=1,n_core_inact_orb
    ii=list_core_inact(i)
    do t=1,n_act_orb
      tt=list_act(t)
      indx+=1
      Tmat(ii,tt)= SXvector(indx)
      Tmat(tt,ii)=-SXvector(indx)
    end do
  end do
  do i=1,n_core_inact_orb
    ii=list_core_inact(i)
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
  


