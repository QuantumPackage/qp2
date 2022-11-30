program truncate_wf
 implicit none
 BEGIN_DOC
! Truncate the wave function
 END_DOC
 read_wf = .True.
 if (s2_eig) then
  call routine_s2
 else
  call routine
 endif
end

subroutine routine
 implicit none
 integer :: ndet_max
 print*, 'Max number of determinants ?'
 read(5,*) ndet_max
 integer(bit_kind), allocatable :: psi_det_tmp(:,:,:)
 double precision, allocatable :: psi_coef_tmp(:,:)
 allocate(psi_det_tmp(N_int,2,ndet_max),psi_coef_tmp(ndet_max, N_states))
 
 integer :: i,j
 double precision :: accu(N_states)
 accu = 0.d0
 do i = 1, ndet_max
  do j = 1, N_int
   psi_det_tmp(j,1,i) = psi_det_sorted(j,1,i)
   psi_det_tmp(j,2,i) = psi_det_sorted(j,2,i)
  enddo
  do j = 1, N_states
   psi_coef_tmp(i,j) = psi_coef_sorted(i,j)
   accu(j) += psi_coef_tmp(i,j) **2
  enddo
 enddo
 do j = 1, N_states
  accu(j) = 1.d0/dsqrt(accu(j))
 enddo
 do j = 1, N_states
  do i = 1, ndet_max
   psi_coef_tmp(i,j) = psi_coef_tmp(i,j) * accu(j)
  enddo
 enddo

 call save_wavefunction_general(ndet_max,N_states,psi_det_tmp,size(psi_coef_tmp,1),psi_coef_tmp)
 
end

subroutine routine_s2
 implicit none
 integer :: ndet_max
 double precision :: wmin
 integer(bit_kind), allocatable :: psi_det_tmp(:,:,:)
 double precision, allocatable :: psi_coef_tmp(:,:)
 integer :: i,j,k
 double precision :: accu(N_states)

 integer :: weights(0:16), ix
 double precision :: x

 weights(:) = 0
 do i=1,N_det
   x = -dlog(1.d-32+sum(weight_configuration(det_to_configuration(i),:)))/dlog(10.d0)
   ix = min(int(x), 16)
   weights(ix) += 1
 enddo

 print *,  'Histogram of the weights of the CFG'
 do i=0,15
   print *, '  10^{-', i, '}   ', weights(i) 
 end do
 print *, '< 10^{-', 15, '}   ', weights(16) 


 print*, 'Min weight of the configuration?'
 read(5,*) wmin

 ndet_max = 0
 do i=1,N_det
   if (maxval(weight_configuration( det_to_configuration(i),:)) < wmin) cycle
   ndet_max = ndet_max+1
 enddo

 allocate(psi_det_tmp(N_int,2,ndet_max),psi_coef_tmp(ndet_max, N_states))
 
 accu = 0.d0
 k=0
 do i = 1, N_det
  if (maxval(weight_configuration( det_to_configuration(i),:)) < wmin) cycle
  k = k+1
  do j = 1, N_int
   psi_det_tmp(j,1,k) = psi_det(j,1,i)
   psi_det_tmp(j,2,k) = psi_det(j,2,i)
  enddo
  do j = 1, N_states
   psi_coef_tmp(k,j) = psi_coef(i,j)
   accu(j) += psi_coef_tmp(k,j)**2
  enddo
 enddo
 do j = 1, N_states
  accu(j) = 1.d0/dsqrt(accu(j))
 enddo
 do j = 1, N_states
  do i = 1, ndet_max
   psi_coef_tmp(i,j) = psi_coef_tmp(i,j) * accu(j)
  enddo
 enddo

 call save_wavefunction_general(ndet_max,N_states,psi_det_tmp,size(psi_coef_tmp,1),psi_coef_tmp)
 
end
