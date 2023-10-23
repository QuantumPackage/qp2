 BEGIN_PROVIDER [ double precision, lccsd_coef, (N_det, N_states)]
&BEGIN_PROVIDER [ double precision, lccsd_energies, (N_states)]
 implicit none
 double precision, allocatable :: Dress_jj(:), H_jj(:), u_in(:,:)
 double precision :: ebefore, eafter, ecorr, thresh
 integer :: i,it,degree
 logical :: converged
 external H_u_0_nstates_openmp
 allocate(Dress_jj(N_det),H_jj(N_det),u_in(N_det,N_states_diag))
 thresh = 1.d-6
 converged = .False.
 Dress_jj = 0.d0
 u_in = 0.d0
 it = 0
 ! initial guess
 do i = 1, N_states_diag
  u_in(i,i) = 1.d0
 enddo
 do i = 1,N_det
  call i_H_j(psi_det(1,1,i),psi_det(1,1,i),N_int,H_jj(i))
 enddo
 ebefore = H_jj(1)
 do while (.not.converged)
  it += 1
  print*,'N_det = ',N_det
  call davidson_general_ext_rout_diag_dressed(u_in,H_jj,Dress_jj,lccsd_energies,& 
                   N_det,N_states,N_states_diag,converged,H_u_0_nstates_openmp)
  ecorr = lccsd_energies(1) - H_jj(1) 
  print*,'---------------------'
  print*,'it = ',it
  print*,'ecorr = ',ecorr
  Dress_jj(1) = 0.d0
  do i = 2, N_det
    if(ecorr + H_jj(i) .lt. H_jj(1))then
     print*,'Warning, some dets are not dressed: ' 
     call get_excitation_degree(ref_bitmask,psi_det(1,1,i),degree,N_int)
     print*,'degree, Delta E, coef', degree, H_jj(i)-H_jj(1), u_in(i,1)/u_in(1,1)
    else
     Dress_jj(i) = ecorr
    endif
  enddo
  eafter = lccsd_energies(1)
  converged = (dabs(eafter - ebefore).lt.thresh)
  ebefore = eafter
 enddo
 do i = 1, N_states
  lccsd_coef(1:N_det,i) = u_in(1:N_det,i)
 enddo

END_PROVIDER 
