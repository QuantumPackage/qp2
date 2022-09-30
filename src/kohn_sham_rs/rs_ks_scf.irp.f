program rs_ks_scf
  BEGIN_DOC
! Produce `Range_separated_Kohn_Sham` MO orbital
! output: mo_basis.mo_num mo_basis.mo_label mo_basis.ao_md5 mo_basis.mo_coef mo_basis.mo_occ
! output: kohn_sham.energy
! optional: mo_basis.mo_coef
  END_DOC

  io_mo_one_e_integrals = "None"
  touch io_mo_one_e_integrals
  io_ao_one_e_integrals = "None"
  touch io_ao_one_e_integrals

  read_wf = .False.
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
 implicit none
 integer :: ifound_x,ifound_c
 if(exchange_functional.eq."None")then
  ifound_x = 1
 else
  ifound_x = index(exchange_functional,"sr")
 endif

 if(correlation_functional.eq."None")then
  ifound_c = 1
 else
  ifound_c = index(correlation_functional,"sr")
 endif
 print*,ifound_x,ifound_c
 if(ifound_x .eq.0 .or. ifound_c .eq. 0)then
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
    print*,'Creating a guess for the MOs'
    print*,'mo_guess_type = ',mo_guess_type
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

  EHF = RS_KS_energy

  mo_label = "Orthonormalized"

  level_shift += 1.d0
  touch level_shift
  call Roothaan_Hall_SCF
  call ezfio_set_kohn_sham_rs_energy(SCF_energy)

 write(*, '(A22,X,F16.10)') 'one_e_energy = ',one_e_energy
 write(*, '(A22,X,F16.10)') 'two_e_energy = ',two_e_energy
 write(*, '(A22,X,F16.10)') 'e_exchange_dft      = ',e_exchange_dft
 write(*, '(A22,X,F16.10)') 'e_correlation_dft   = ',e_correlation_dft
 write(*, '(A22,X,F16.10)') 'Fock_matrix_energy  = ',Fock_matrix_energy


end


