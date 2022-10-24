BEGIN_PROVIDER [ double precision, extrapolated_energy, (N_iter,N_states) ]
 implicit none
 BEGIN_DOC
 ! Extrapolated energy, using E_var = f(PT2) where PT2=0
 END_DOC
! integer :: i
 extrapolated_energy = 0.D0
END_PROVIDER 

 subroutine get_extrapolated_energy(Niter,ept2,pt1,extrap_energy)
 implicit none
 integer, intent(in)  :: Niter
 double precision, intent(in) :: ept2(Niter),pt1(Niter),extrap_energy(Niter)
 call extrapolate_data(Niter,ept2,pt1,extrap_energy)
 end

subroutine save_iterations(e_, pt2_,n_)
  implicit none
  BEGIN_DOC
! Update the energy in the EZFIO file.
  END_DOC
  integer, intent(in) :: n_
  double precision, intent(in) :: e_(N_states), pt2_(N_states)
  integer :: i

  if (N_iter == 101) then
    do i=2,N_iter-1
      energy_iterations(1:N_states,N_iter-1) = energy_iterations(1:N_states,N_iter)
      pt2_iterations(1:N_states,N_iter-1) = pt2_iterations(1:N_states,N_iter) 
    enddo
    N_iter = N_iter-1
    TOUCH N_iter
  endif

  energy_iterations(1:N_states,N_iter) = e_(1:N_states)
     pt2_iterations(1:N_states,N_iter) = pt2_(1:N_states)
  n_det_iterations(N_iter) = n_
  call ezfio_set_iterations_N_iter(N_iter)
  call ezfio_set_iterations_energy_iterations(energy_iterations)
  call ezfio_set_iterations_pt2_iterations(pt2_iterations)
  call ezfio_set_iterations_n_det_iterations(n_det_iterations)
end

