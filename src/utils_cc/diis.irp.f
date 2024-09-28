! Code

subroutine diis_cc(all_err,all_t,sze,m,iter,t)

  implicit none

  BEGIN_DOC
  ! DIIS. Take the error vectors and the amplitudes of the previous
  ! iterations to compute the new amplitudes
  END_DOC
  
  ! {err_i}_{i=1}^{m_it} -> B -> c
  ! {t_i}_{i=1}^{m_it}, c, {err_i}_{i=1}^{m_it} -> t_{m_it+1}

  integer, intent(in)             :: m,iter,sze
  double precision, intent(in)    :: all_err(sze,m)
  double precision, intent(in)    :: all_t(sze,m)
  
  double precision, intent(out)   :: t(sze)
  
  double precision, allocatable   :: B(:,:), c(:), zero(:)
  integer                         :: m_iter
  integer                         :: i,j,k
  integer                         :: info
  integer, allocatable            :: ipiv(:)
  double precision                :: accu
  
  m_iter = min(m,iter)
  !print*,'m_iter',m_iter
  allocate(B(m_iter+1,m_iter+1), c(m_iter), zero(m_iter+1))
  allocate(ipiv(m+1))

  ! B(i,j) =  < err(iter-m_iter+j),err(iter-m_iter+i) > ! iter-m_iter will be zero for us
  B = 0d0
  !$OMP PARALLEL &
  !$OMP SHARED(B,m,m_iter,sze,all_err) &
  !$OMP PRIVATE(i,j,k,accu) &
  !$OMP DEFAULT(NONE)
  do j = 1, m_iter
    do i = 1, m_iter
      accu = 0d0
      !$OMP DO
      do k = 1, sze
        ! the errors of the ith iteration are in all_err(:,m+1-i)
        accu = accu + all_err(k,m+1-i) * all_err(k,m+1-j)
      enddo
      !$OMP END DO NOWAIT
      !$OMP CRITICAL
      B(i,j) = B(i,j) + accu
      !$OMP END CRITICAL
    enddo
  enddo
  !$OMP END PARALLEL
  
  do i = 1, m_iter
    B(i,m_iter+1) = -1.d0
  enddo
  do j = 1, m_iter
    B(m_iter+1,j) = -1.d0
  enddo
  ! Debug
  !print*,'B'
  !do i = 1, m_iter+1
  !  write(*,'(100(F10.6))') B(i,:)
  !enddo

  ! (0 0 .... 0 -1)
  zero = 0d0
  zero(m_iter+1) = -1d0

  ! Solve B.c = zero
  call dgesv(m_iter+1, 1, B, size(B,1), ipiv, zero, size(zero,1), info)
  if (info /= 0) then
    print*,'DIIS error in dgesv:', info
    call abort
  endif
  ! c corresponds to the m_iter first solutions
  c = zero(1:m_iter)
  ! Debug
  !print*,'c',c
  !print*,'all_t' 
  !do i = 1, m
  !  write(*,'(100(F10.6))') all_t(:,i)
  !enddo
  !print*,'all_err' 
  !do i = 1, m
  !  write(*,'(100(F10.6))') all_err(:,i)
  !enddo

  ! update T
  !$OMP PARALLEL &
  !$OMP SHARED(t,c,m,all_err,all_t,sze,m_iter) &
  !$OMP PRIVATE(i,j,accu) &
  !$OMP DEFAULT(NONE)
  !$OMP DO
  do i = 1, sze
    t(i) = 0d0
  enddo
  !$OMP END DO
  do i = 1, m_iter
    !$OMP DO
    do j = 1, sze
      t(j) = t(j) + c(i) * (all_t(j,m+1-i) + all_err(j,m+1-i))
    enddo
    !$OMP END DO
  enddo
  !$OMP END PARALLEL

  !print*,'new t',t

  deallocate(ipiv,B,c,zero)

end

! Update all err

subroutine update_all_err(err,all_err,sze,m,iter)

  implicit none

  BEGIN_DOC
  ! Shift all the err vectors of the previous iterations to add the new one
  ! The last err vector is placed in the last position and all the others are
  ! moved toward the first one.
  END_DOC

  integer, intent(in)             :: m, iter, sze
  double precision, intent(in)    :: err(sze)
  double precision, intent(inout) :: all_err(sze,m)
  integer                         :: i,j
  integer                         :: m_iter

  m_iter = min(m,iter)

  ! Shift
  !$OMP PARALLEL &
  !$OMP SHARED(m,all_err,err,sze) &
  !$OMP PRIVATE(i,j) &
  !$OMP DEFAULT(NONE)
  do i = 1, m-1
    !$OMP DO
    do j = 1, sze
      all_err(j,i) = all_err(j,i+1)
    enddo
    !$OMP END DO
  enddo
  
  ! Debug
  !print*,'shift err'
  !do i = 1, m
  !  print*,i, all_err(:,i)
  !enddo

  ! New
  !$OMP DO
  do i = 1, sze
    all_err(i,m) = err(i)
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  ! Debug
  !print*,'Updated err'
  !do i = 1, m
  !  print*,i, all_err(:,i)
  !enddo

end

! Update all t

subroutine update_all_t(t,all_t,sze,m,iter)

  implicit none

  BEGIN_DOC
  ! Shift all the t vectors of the previous iterations to add the new one
  ! The last t vector is placed in the last position and all the others are
  ! moved toward the first one.
  END_DOC

  integer, intent(in)             :: m, iter, sze
  double precision, intent(in)    :: t(sze)
  double precision, intent(inout) :: all_t(sze,m)
  integer                         :: i,j
  integer                         :: m_iter

  m_iter = min(m,iter)

  ! Shift
  !$OMP PARALLEL &
  !$OMP SHARED(m,all_t,t,sze) &
  !$OMP PRIVATE(i,j) &
  !$OMP DEFAULT(NONE)
  do i = 1, m-1
    !$OMP DO
    do j = 1, sze
      all_t(j,i) = all_t(j,i+1)
    enddo
    !$OMP END DO
  enddo

  ! New
  !$OMP DO
  do i = 1, sze
    all_t(i,m) = t(i)
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  ! Debug
  !print*,'Updated t'
  !do i = 1, m
  !  print*,i, all_t(:,i)
  !enddo

end

! Err1

subroutine compute_err1(nO,nV,f_o,f_v,r1,err1)

  implicit none

  BEGIN_DOC
  ! Compute the error vector for the t1
  END_DOC

  integer, intent(in)           :: nO, nV
  double precision, intent(in)  :: f_o(nO), f_v(nV), r1(nO,nV)
  
  double precision, intent(out) :: err1(nO,nV)

  integer                       :: i,a
  
  !$OMP PARALLEL &
  !$OMP SHARED(err1,r1,f_o,f_v,nO,nV,cc_level_shift) &
  !$OMP PRIVATE(i,a) &
  !$OMP DEFAULT(NONE)
  !$OMP DO
  do a = 1, nV
    do i = 1, nO
      err1(i,a) = - r1(i,a) / (f_o(i) - f_v(a) - cc_level_shift)
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

end

! Err2

subroutine compute_err2(nO,nV,f_o,f_v,r2,err2)

  implicit none

  BEGIN_DOC
  ! Compute the error vector for the t2
  END_DOC

  integer, intent(in)           :: nO, nV
  double precision, intent(in)  :: f_o(nO), f_v(nV), r2(nO,nO,nV,nV)
  
  double precision, intent(out) :: err2(nO,nO,nV,nV)

  integer                       :: i,j,a,b

  !$OMP PARALLEL &
  !$OMP SHARED(err2,r2,f_o,f_v,nO,nV,cc_level_shift) &
  !$OMP PRIVATE(i,j,a,b) &
  !$OMP DEFAULT(NONE)
  !$OMP DO collapse(3)
  do b = 1, nV
    do a = 1, nV
      do j = 1, nO
        do i = 1, nO       
          err2(i,j,a,b) = - r2(i,j,a,b) / (f_o(i) + f_o(j) - f_v(a) - f_v(b) - cc_level_shift)
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

end

! Update t

subroutine update_t_ccsd(nO,nV,nb_iter,f_o,f_v,r1,r2,t1,t2,all_err1,all_err2,all_t1,all_t2)

  implicit none

  integer, intent(in)             :: nO,nV,nb_iter
  double precision, intent(in)    :: f_o(nO), f_v(nV)
  double precision, intent(in)    :: r1(nO,nV), r2(nO,nO,nV,nV)
  
  double precision, intent(inout) :: t1(nO,nV), t2(nO,nO,nV,nV)
  double precision, intent(inout) :: all_err1(nO*nV, cc_diis_depth), all_err2(nO*nO*nV*nV, cc_diis_depth)
  double precision, intent(inout) :: all_t1(nO*nV, cc_diis_depth), all_t2(nO*nO*nV*nV, cc_diis_depth)

  double precision, allocatable   :: err1(:,:), err2(:,:,:,:)
  double precision, allocatable   :: tmp_err1(:), tmp_err2(:)
  double precision, allocatable   :: tmp_t1(:), tmp_t2(:)
  
  if (cc_update_method == 'diis') then

    allocate(err1(nO,nV), err2(nO,nO,nV,nV))
    allocate(tmp_err1(nO*nV), tmp_err2(nO*nO*nV*nV))
    allocate(tmp_t1(nO*nV), tmp_t2(nO*nO*nV*nV))

    ! DIIS T1, it is not always good since the t1 can be small
    ! That's why there is a call to update the t1 in the standard way
    ! T1 error tensor
    !call compute_err1(nO,nV,f_o,f_v,r1,err1)
    ! Transfo errors and parameters in vectors
    !tmp_err1 = reshape(err1,(/nO*nV/))
    !tmp_t1   = reshape(t1  ,(/nO*nV/))
    ! Add the error and parameter vectors with those of the previous iterations
    !call update_all_err(tmp_err1,all_err1,nO*nV,cc_diis_depth,nb_iter+1)
    !call update_all_t  (tmp_t1  ,all_t1  ,nO*nV,cc_diis_depth,nb_iter+1)
    ! Diis and reshape T as a tensor
    !call diis_cc(all_err1,all_t1,nO*nV,cc_diis_depth,nb_iter+1,tmp_t1)
    !t1 = reshape(tmp_t1  ,(/nO,nV/))
    call update_t1(nO,nV,f_o,f_v,r1,t1)

    ! DIIS T2
    ! T2 error tensor
    call compute_err2(nO,nV,f_o,f_v,r2,err2)
    ! Transfo errors and parameters in vectors
    tmp_err2 = reshape(err2,(/nO*nO*nV*nV/))
    tmp_t2   = reshape(t2  ,(/nO*nO*nV*nV/))
    ! Add the error and parameter vectors with those of the previous iterations
    call update_all_err(tmp_err2,all_err2,nO*nO*nV*nV,cc_diis_depth,nb_iter+1)
    call update_all_t  (tmp_t2  ,all_t2  ,nO*nO*nV*nV,cc_diis_depth,nb_iter+1)
    ! Diis and reshape T as a tensor
    call diis_cc(all_err2,all_t2,nO*nO*nV*nV,cc_diis_depth,nb_iter+1,tmp_t2)
    t2 = reshape(tmp_t2  ,(/nO,nO,nV,nV/))

    deallocate(tmp_t1,tmp_t2,tmp_err1,tmp_err2,err1,err2)

  ! Standard update as T = T - Delta
  elseif (cc_update_method == 'none') then
     
    call update_t1(nO,nV,f_o,f_v,r1,t1)
    call update_t2(nO,nV,f_o,f_v,r2,t2)
    
  else
    print*,'Unkonw cc_method_method: '//cc_update_method
  endif
  
end

! Update t v2

subroutine update_t_ccsd_diis(nO,nV,nb_iter,f_o,f_v,r1,r2,t1,t2,all_err1,all_err2,all_t1,all_t2)

  implicit none

  integer, intent(in)             :: nO,nV,nb_iter
  double precision, intent(in)    :: f_o(nO), f_v(nV)
  double precision, intent(in)    :: r1(nO,nV), r2(nO,nO,nV,nV)
  
  double precision, intent(inout) :: t1(nO,nV), t2(nO,nO,nV,nV)
  double precision, intent(inout) :: all_err1(nO*nV, cc_diis_depth), all_err2(nO*nO*nV*nV, cc_diis_depth)
  double precision, intent(inout) :: all_t1(nO*nV, cc_diis_depth), all_t2(nO*nO*nV*nV, cc_diis_depth)

  double precision, allocatable   :: all_t(:,:), all_err(:,:), tmp_t(:)
  double precision, allocatable   :: err1(:,:), err2(:,:,:,:)
  double precision, allocatable   :: tmp_err1(:), tmp_err2(:)
  double precision, allocatable   :: tmp_t1(:), tmp_t2(:)

  integer                         :: i,j
  
  ! Allocate
  allocate(all_err(nO*nV+nO*nO*nV*nV,cc_diis_depth), all_t(nO*nV+nO*nO*nV*nV,cc_diis_depth))
  allocate(tmp_t(nO*nV+nO*nO*nV*nV))
  allocate(err1(nO,nV), err2(nO,nO,nV,nV))
  allocate(tmp_err1(nO*nV), tmp_err2(nO*nO*nV*nV))
  allocate(tmp_t1(nO*nV), tmp_t2(nO*nO*nV*nV))

  ! Compute the errors and reshape them as vector
  call compute_err1(nO,nV,f_o,f_v,r1,err1)
  call compute_err2(nO,nV,f_o,f_v,r2,err2)
  tmp_err1 = reshape(err1,(/nO*nV/))
  tmp_err2 = reshape(err2,(/nO*nO*nV*nV/))
  tmp_t1   = reshape(t1  ,(/nO*nV/))
  tmp_t2   = reshape(t2  ,(/nO*nO*nV*nV/))
  
  ! Update the errors and parameters for the diis
  call update_all_err(tmp_err1,all_err1,nO*nV,cc_diis_depth,nb_iter+1)
  call update_all_t  (tmp_t1  ,all_t1  ,nO*nV,cc_diis_depth,nb_iter+1)
  call update_all_err(tmp_err2,all_err2,nO*nO*nV*nV,cc_diis_depth,nb_iter+1)
  call update_all_t  (tmp_t2  ,all_t2  ,nO*nO*nV*nV,cc_diis_depth,nb_iter+1)

  ! Gather the different parameters and errors
  !$OMP PARALLEL &
  !$OMP SHARED(nO,nV,all_err,all_err1,all_err2,cc_diis_depth,&
  !$OMP all_t,all_t1,all_t2) &
  !$OMP PRIVATE(i,j) &
  !$OMP DEFAULT(NONE)
  do j = 1, cc_diis_depth
    !$OMP DO 
    do i = 1, nO*nV
      all_err(i,j) = all_err1(i,j)
    enddo
    !$OMP END DO NOWAIT
  enddo
  do j = 1, cc_diis_depth
    !$OMP DO
    do i = 1, nO*nO*nV*nV
      all_err(i+nO*nV,j) = all_err2(i,j)
    enddo
    !$OMP END DO NOWAIT
  enddo
  do j = 1, cc_diis_depth
    !$OMP DO 
    do i = 1, nO*nV
      all_t(i,j) = all_t1(i,j)
    enddo
    !$OMP END DO NOWAIT
  enddo
  do j = 1, cc_diis_depth
    !$OMP DO 
    do i = 1, nO*nO*nV*nV
      all_t(i+nO*nV,j) = all_t2(i,j)
    enddo
    !$OMP END DO
  enddo
  !$OMP END PARALLEL
  
  ! Diis
  call diis_cc(all_err,all_t,nO*nV+nO*nO*nV*nV,cc_diis_depth,nb_iter+1,tmp_t)

  ! Split the resulting vector
  !$OMP PARALLEL &
  !$OMP SHARED(nO,nV,tmp_t,tmp_t1,tmp_t2) &
  !$OMP PRIVATE(i) &
  !$OMP DEFAULT(NONE)
  !$OMP DO
  do i = 1, nO*nV
    tmp_t1(i) = tmp_t(i)
  enddo
  !$OMP END DO NOWAIT
  !$OMP DO
  do i = 1, nO*nO*nV*nV
    tmp_t2(i) = tmp_t(i+nO*nV) 
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  ! Reshape as tensors
  t1 = reshape(tmp_t1 ,(/nO,nV/))
  t2 = reshape(tmp_t2 ,(/nO,nO,nV,nV/))

  ! Deallocate
  deallocate(tmp_t1,tmp_t2,tmp_err1,tmp_err2,err1,err2,all_t,all_err)

end

! Update t v3

subroutine update_t_ccsd_diis_v3(nO,nV,nb_iter,f_o,f_v,r1,r2,t1,t2,all_err,all_t)

  implicit none

  integer, intent(in)             :: nO,nV,nb_iter
  double precision, intent(in)    :: f_o(nO), f_v(nV)
  double precision, intent(in)    :: r1(nO,nV), r2(nO,nO,nV,nV)
  
  double precision, intent(inout) :: t1(nO*nV), t2(nO*nO*nV*nV)
  double precision, intent(inout) :: all_err(nO*nV+nO*nO*nV*nV, cc_diis_depth)
  double precision, intent(inout) :: all_t(nO*nV+nO*nO*nV*nV, cc_diis_depth)

  double precision, allocatable   :: tmp(:)

  integer                         :: i,j
  
  ! Allocate
  allocate(tmp(nO*nV+nO*nO*nV*nV))

  ! Compute the errors
  call compute_err1(nO,nV,f_o,f_v,r1,tmp(1:nO*nV))
  call compute_err2(nO,nV,f_o,f_v,r2,tmp(nO*nV+1:nO*nV+nO*nO*nV*nV))
  
  ! Update the errors and parameters for the diis
  call update_all_err(tmp,all_err,nO*nV+nO*nO*nV*nV,cc_diis_depth,nb_iter+1)

  !$OMP PARALLEL &
  !$OMP SHARED(nO,nV,tmp,t1,t2) &
  !$OMP PRIVATE(i) &
  !$OMP DEFAULT(NONE)
  !$OMP DO
  do i = 1, nO*nV
    tmp(i) = t1(i)
  enddo
  !$OMP END DO
  !$OMP DO
  do i = 1, nO*nO*nV*nV
    tmp(i+nO*nV) = t2(i)
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
   
  call update_all_t(tmp,all_t,nO*nV+nO*nO*nV*nV,cc_diis_depth,nb_iter+1)

  ! Diis
  call diis_cc(all_err,all_t,nO*nV+nO*nO*nV*nV,cc_diis_depth,nb_iter+1,tmp)

  ! Split the resulting vector
  !$OMP PARALLEL &
  !$OMP SHARED(nO,nV,tmp,t1,t2) &
  !$OMP PRIVATE(i) &
  !$OMP DEFAULT(NONE)
  !$OMP DO
  do i = 1, nO*nV
    t1(i) = tmp(i)
  enddo
  !$OMP END DO
  !$OMP DO
  do i = 1, nO*nO*nV*nV
    t2(i) = tmp(i+nO*nV) 
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  ! Deallocate
  deallocate(tmp)

end
