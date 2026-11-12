
! ---

BEGIN_PROVIDER [ integer, n_states_diag  ]
  implicit none
  BEGIN_DOC
! Number of states to consider during the Davdison diagonalization
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then

    call ezfio_has_davidson_keywords_n_states_diag(has)
    if (has) then
      call ezfio_get_davidson_keywords_n_states_diag(n_states_diag)
    else
      print *, 'davidson_keywords/n_states_diag not found in EZFIO file'
      stop 1
    endif
    n_states_diag = max(2,N_states * N_states_diag)
  endif

END_PROVIDER

! ---
