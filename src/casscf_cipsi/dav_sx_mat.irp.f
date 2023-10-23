

subroutine davidson_diag_sx_mat(N_st, u_in, energies)
 implicit none
 integer, intent(in) :: N_st
 double precision, intent(out) :: u_in(nMonoEx+1,n_states_diag), energies(N_st)
 integer :: i,j,N_st_tmp, dim_in, sze, N_st_diag_in
 integer, allocatable :: list_guess(:)
 double precision, allocatable :: H_jj(:)
 logical    :: converged                                                                                                
 N_st_diag_in = n_states_diag
 provide SXmatrix
 sze = nMonoEx+1
 dim_in = sze
 allocate(H_jj(sze), list_guess(sze))
 H_jj(1) = 0.d0
 N_st_tmp = 1
 list_guess(1) = 1
 do j = 2, nMonoEx+1
  H_jj(j) = SXmatrix(j,j)
  if(H_jj(j).lt.0.d0)then 
   list_guess(N_st_tmp) = j
   N_st_tmp += 1
  endif
 enddo
 if(N_st_tmp .ne. N_st)then
  print*,'Pb in davidson_diag_sx_mat'
  print*,'N_st_tmp .ne. N_st'
  print*,N_st_tmp, N_st
  stop
 endif
 print*,'Number of possibly interesting states = ',N_st
 print*,'Corresponding diagonal elements of the SX matrix '
 u_in = 0.d0
 do i = 1, min(N_st, N_st_diag_in)
! do i = 1, N_st
  j = list_guess(i) 
  print*,'i,j',i,j
  print*,'SX(i,i) = ',H_jj(j)
  u_in(j,i) = 1.d0
 enddo
 call davidson_general(u_in,H_jj,energies,dim_in,sze,N_st,N_st_diag_in,converged,SXmatrix)
 print*,'energies = ',energies

end
