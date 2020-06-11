program print_h_omp_debug
 implicit none
 read_wf = .True.
 touch read_wf
 call routine_omp

end

subroutine routine_omp
 use bitmasks
  implicit none
  integer     :: h_size
  complex*16, allocatable :: u_tmp(:,:), s_tmp(:,:),v_tmp(:,:)
  integer :: i,n_st
  h_size=N_det
  BEGIN_DOC
  ! |vec2> = H|vec1>
  !
  ! TODO: implement
  ! maybe reuse parts of H_S2_u_0_nstates_{openmp,zmq}?
  END_DOC
  n_st=min(1000,h_size)
  allocate(u_tmp(n_st,h_size),s_tmp(n_st,h_size),v_tmp(n_st,h_size))

  u_tmp=(0.d0,0.d0)
  v_tmp=(0.d0,0.d0)
  s_tmp=(0.d0,0.d0)

  do i=1,n_st
    u_tmp(i,i)=(1.d0,0.d0)
  enddo

  call h_s2_u_0_nstates_openmp_complex(v_tmp,s_tmp,u_tmp,n_st,h_size)
  do i = 1, n_st
    v_tmp(i,i) += nuclear_repulsion
  enddo
 do i = 1, n_st
  write(*,'(I3,X,A3,2000(E24.15))')i,' | ',v_tmp(i,:)
 enddo
  deallocate(u_tmp,v_tmp,s_tmp)
end
