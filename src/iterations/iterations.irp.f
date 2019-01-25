BEGIN_PROVIDER [ double precision, extrapolated_energy, (N_iter,N_states) ]
 implicit none
 BEGIN_DOC
 ! Extrapolated energy, using E_var = f(PT2) where PT2=0
 END_DOC
 integer :: i
 do i=1,min(N_states,N_det)
   call extrapolate_data(N_iter,                               &
       energy_iterations(i,1:N_iter),                          &
          pt2_iterations(i,1:N_iter),                          &
       extrapolated_energy(1:N_iter,i))
 enddo
END_PROVIDER


subroutine save_iterations(e_, pt2_,n_)
  implicit none
  BEGIN_DOC
! Update the energy in the EZFIO file.
  END_DOC
  integer, intent(in) :: n_
  double precision, intent(in) :: e_(N_states), pt2_(N_states)

  if (N_iter > 100) then
    return
  endif

  energy_iterations(1:N_states,N_iter) = e_(1:N_states)
     pt2_iterations(1:N_states,N_iter) = pt2_(1:N_states)
  n_det_iterations(N_iter) = n_
  call ezfio_set_iterations_N_iter(N_iter)
  call ezfio_set_iterations_energy_iterations(energy_iterations)
  call ezfio_set_iterations_pt2_iterations(pt2_iterations)
  call ezfio_set_iterations_n_det_iterations(n_det_iterations)
end

