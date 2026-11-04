program ks_scf
  BEGIN_DOC
! Produce `Kohn_Sham` MO orbital
! output: mo_basis.mo_num mo_basis.mo_label mo_basis.ao_md5 mo_basis.mo_coef mo_basis.mo_occ
! output: kohn_sham.energy
! optional: mo_basis.mo_coef
  END_DOC

  io_mo_one_e_integrals = "None"
  touch io_mo_one_e_integrals
  io_ao_one_e_integrals = "None"
  touch io_ao_one_e_integrals
  density_for_dft ="KS"
  touch density_for_dft
  print*, '**************************'
  print*, 'mu_erf_dft = ',mu_erf_dft
  print*, '**************************'
  call check_coherence_functional
  call create_guess
  call orthonormalize_mos
  call run
end

subroutine check_coherence_functional
  BEGIN_DOC
  ! ks_scf start-up check: stop with a diagnostic if
  ! exchange_functional or correlation_functional is "None" or
  ! contains the substring "short_range".
  ! TODO(doc): intent unknown to the author (author, 2026-10-31).
  ! The code literally stops when a functional is "None" or
  ! contains the substring "short_range"; in-tree functionals are
  ! named sr_*, so the substring branch is dead for them, and the
  ! diagnostic text is copy-pasted from
  ! kohn_sham_rs/rs_ks_scf.irp.f (there the test is "sr").
  END_DOC
 implicit none
 integer :: ifound_x,ifound_c
 if(exchange_functional.eq."None")then
  ifound_x = 1
 else
  ifound_x = index(exchange_functional,"short_range")
 endif

 if(correlation_functional.eq."None")then
  ifound_c = 1
 else
  ifound_c = index(correlation_functional,"short_range")
 endif
 print*,ifound_x,ifound_c
 if(ifound_x .ne.0 .or. ifound_c .ne. 0)then
  print*,'YOU ARE USING THE RANGE SEPARATED KS PROGRAM BUT YOUR INPUT KEYWORD FOR '
  print*,'exchange_functional is ',exchange_functional
  print*,'correlation_functional is ',correlation_functional
  print*,'CHANGE THE exchange_functional and correlation_functional keywords to range separated functionals'
  print*,'or switch to the KS_SCF program that uses regular functionals'
  stop
 endif

end



subroutine create_guess
  implicit none
  BEGIN_DOC
!   Create a MO guess if no MOs are present in the EZFIO directory
  END_DOC
  logical                        :: exists
  PROVIDE ezfio_filename
  call ezfio_has_mo_basis_mo_coef(exists)
  if (.not.exists) then
    if (mo_guess_type == "HCore") then
      mo_coef = ao_ortho_lowdin_coef
      TOUCH mo_coef
      mo_label = 'Guess'
      call mo_as_eigvectors_of_mo_matrix(mo_one_e_integrals,size(mo_one_e_integrals,1),size(mo_one_e_integrals,2),mo_label,.false.)
      SOFT_TOUCH mo_coef mo_label
    else if (mo_guess_type == "Huckel") then
      call huckel_guess
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

  double precision               :: EHF

  EHF = KS_energy

  mo_label = "Orthonormalized"

! Choose SCF algorithm

  write(json_unit,*) '"scf" : ['
  call Roothaan_Hall_SCF
  write(json_unit,*) ']'

  call json_close

  call pg_print_mo_energies(Fock_matrix_diag_mo, mo_num)
  call pg_report_mo_symmetry()
end


