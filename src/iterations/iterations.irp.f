BEGIN_PROVIDER [ integer, N_iter  ]
  implicit none
  BEGIN_DOC
! Number of CIPSI iterations
  END_DOC

  N_iter = 0
END_PROVIDER

BEGIN_PROVIDER [ integer, N_iter_max ]
 implicit none
 BEGIN_DOC
 ! Max number of iterations for extrapolations
 END_DOC
 N_iter_max = 8
END_PROVIDER

 BEGIN_PROVIDER [ double precision, energy_iterations , (n_states,N_iter_max) ]
&BEGIN_PROVIDER [ double precision, pt2_iterations , (n_states,N_iter_max) ]
&BEGIN_PROVIDER [ double precision, extrapolated_energy, (N_iter_max,N_states) ]
  implicit none
  BEGIN_DOC
! The energy at each iteration for the extrapolations
  END_DOC

   energy_iterations = 0.d0
   pt2_iterations = 0.d0
   extrapolated_energy = 0.d0
END_PROVIDER

subroutine increment_n_iter(e, pt2_data)
  use selection_types
  implicit none
  BEGIN_DOC
! Does what is necessary to increment n_iter
  END_DOC
  double precision, intent(in) :: e(*)
  type(pt2_type), intent(in)   :: pt2_data
  integer :: k, i

  if (N_det < N_states) return

  if (N_iter < N_iter_max) then
    N_iter += 1
  else
    do k=2,N_iter
      energy_iterations(1:N_states,k-1) = energy_iterations(1:N_states,k)
      pt2_iterations(1:N_states,k-1) = pt2_iterations(1:N_states,k)
    enddo
  endif
  energy_iterations(1:N_states,N_iter) = e(1:N_states)
  pt2_iterations(1:N_states,N_iter) = pt2_data % rpt2(1:N_states)

  if (N_iter < 2) then
    extrapolated_energy(1,:) = energy_iterations(:,1) + pt2_iterations(:,1)
    extrapolated_energy(2,:) = energy_iterations(:,2) + pt2_iterations(:,2)
  else
    do i=1,N_states
      call extrapolate_data(N_iter,                               &
          energy_iterations(i,1:N_iter),                          &
             pt2_iterations(i,1:N_iter),                          &
          extrapolated_energy(1:N_iter,i))
    enddo
  endif
end
