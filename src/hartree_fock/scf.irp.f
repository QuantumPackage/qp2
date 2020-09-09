program scf
  BEGIN_DOC
! 
! The :ref:`scf` program performs *Restricted* Hartree-Fock
! calculations (the spatial part of the |MOs| is common for alpha and beta
! spinorbitals).
! 
! It performs the following actions:
! 
! #. Compute/Read all the one- and two-electron integrals, and store them
!    in memory
! #. Check in the |EZFIO| database if there is a set of |MOs|.
!    If there is, it will read them as initial guess. Otherwise, it will
!    create a guess.
! #. Perform the |SCF| iterations
! 
! For the keywords related to the |SCF| procedure, see the ``scf_utils``
! directory where you will find all options.
! 
! At each iteration, the |MOs| are saved in the |EZFIO| database. Hence,
! if the calculation crashes for any unexpected reason, the calculation
! can be restarted by running again the |SCF| with the same |EZFIO|
! database.
! 
! To start again a fresh |SCF| calculation, the |MOs| can be reset by
! running the :ref:`qp_reset` command.
! 
! The `DIIS`_ algorithm is implemented, as well as the `level-shifting`_ 
! method. If the |SCF| does not converge, try again with a higher value of
! :option:`level_shift`.
! 
! .. _DIIS: https://en.wikipedia.org/w/index.php?title=DIIS
! .. _level-shifting: https://doi.org/10.1002/qua.560070407
!
  END_DOC
  call create_guess
  call orthonormalize_mos
  call run
end

subroutine create_guess
  implicit none
  BEGIN_DOC
!   Create a MO guess if no MOs are present in the EZFIO directory
  END_DOC
  logical                        :: exists
  PROVIDE ezfio_filename
  if (is_complex) then
!    call ezfio_has_mo_basis_mo_coef_complex(exists)
    call ezfio_has_mo_basis_mo_coef_kpts(exists)
  else
    call ezfio_has_mo_basis_mo_coef(exists)
  endif
  if (.not.exists) then
    if (mo_guess_type == "HCore") then
      if (is_complex) then
        !mo_coef_complex = ao_ortho_lowdin_coef_complex
        mo_coef_kpts = ao_ortho_lowdin_coef_kpts
        TOUCH mo_coef_kpts
        mo_label = 'Guess'
        !call mo_as_eigvectors_of_mo_matrix_complex(mo_one_e_integrals_kpts,     &
        call mo_as_eigvectors_of_mo_matrix_kpts(mo_one_e_integrals_kpts,     &
            size(mo_one_e_integrals_kpts,1),                            &
            size(mo_one_e_integrals_kpts,2),                            &
            size(mo_one_e_integrals_kpts,3),                            &
            mo_label,1,.false.)
        SOFT_TOUCH mo_coef_kpts mo_label
      else
        mo_coef = ao_ortho_lowdin_coef
        TOUCH mo_coef
        mo_label = 'Guess'
        call mo_as_eigvectors_of_mo_matrix(mo_one_e_integrals,     &
            size(mo_one_e_integrals,1),                            &
            size(mo_one_e_integrals,2),                            &
            mo_label,1,.false.)
        SOFT_TOUCH mo_coef mo_label
      endif
    else if (mo_guess_type == "Huckel") then
      if (is_complex) then
        !call huckel_guess_complex
        call huckel_guess_kpts
      else
        call huckel_guess
      endif
    else
      print *,  'Unrecognized MO guess type : '//mo_guess_type
      stop 1
    endif
  endif
end

subroutine run

  BEGIN_DOC
!   Run SCF calculation
  END_DOC

  use bitmasks
  implicit none

  integer                        :: i_it, i, j, k

  mo_label = "Orthonormalized"
  if (is_complex) then
    !call roothaan_hall_scf_complex
    call roothaan_hall_scf_kpts
  else
    call roothaan_hall_scf
  endif
  call ezfio_set_hartree_fock_energy(SCF_energy)
  print*,'hf 1e,2e,total energy'
  print*,hf_one_electron_energy
  print*,hf_two_electron_energy
  print*,hf_energy

end


