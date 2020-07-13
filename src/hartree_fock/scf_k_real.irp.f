program scf_k_real
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
  call create_guess_k_real
  call orthonormalize_mos_k_real
  call run_k_real
end

subroutine create_guess_k_real
  implicit none
  BEGIN_DOC
!   Create a MO guess if no MOs are present in the EZFIO directory
  END_DOC
  logical                        :: exists
  PROVIDE ezfio_filename
  call ezfio_has_mo_basis_mo_coef_kpts(exists)
  if (.not.exists) then
    if (mo_guess_type == "HCore") then
        !mo_coef_complex = ao_ortho_lowdin_coef_complex
        mo_coef_kpts = ao_ortho_lowdin_coef_kpts_real
        TOUCH mo_coef_kpts
        mo_label = 'Guess'
        !call mo_as_eigvectors_of_mo_matrix_complex(mo_one_e_integrals_kpts,     &
        call mo_as_eigvectors_of_mo_matrix_kpts_real(mo_one_e_integrals_kpts_real,     &
            size(mo_one_e_integrals_kpts_real,1),                            &
            size(mo_one_e_integrals_kpts_real,2),                            &
            size(mo_one_e_integrals_kpts_real,3),                            &
            mo_label,1,.false.)
        SOFT_TOUCH mo_coef_kpts mo_label
    else if (mo_guess_type == "Huckel") then
        call huckel_guess_kpts_real
    else
      print *,  'Unrecognized MO guess type : '//mo_guess_type
      stop 1
    endif
  endif
end

subroutine run_k_real

  BEGIN_DOC
!   Run SCF calculation
  END_DOC

  use bitmasks
  implicit none

  integer                        :: i_it, i, j, k

  mo_label = "Orthonormalized"
  call roothaan_hall_scf_kpts_real
  call ezfio_set_hartree_fock_energy(SCF_energy)
  print*,'hf 1e,2e,total energy'
  print*,hf_one_electron_energy
  print*,hf_two_electron_energy
  print*,hf_energy

end


