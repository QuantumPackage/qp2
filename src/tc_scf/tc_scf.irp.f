program tc_scf

  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC

  implicit none

  print *, 'starting ...'

  my_grid_becke  = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
!  my_n_pt_r_grid = 10 ! small grid for quick debug
!  my_n_pt_a_grid = 26 ! small grid for quick debug
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  !call create_guess
  !call orthonormalize_mos

  call routine_scf()

end

! ---

subroutine create_guess

  BEGIN_DOC
  !   Create a MO guess if no MOs are present in the EZFIO directory
  END_DOC

  implicit none
  logical :: exists

  PROVIDE ezfio_filename
  call ezfio_has_mo_basis_mo_coef(exists)

  if (.not.exists) then
    mo_label = 'Guess'
    if (mo_guess_type == "HCore") then
      mo_coef = ao_ortho_lowdin_coef
      call restore_symmetry(ao_num, mo_num, mo_coef, size(mo_coef, 1), 1.d-10)
      TOUCH mo_coef
      call mo_as_eigvectors_of_mo_matrix(mo_one_e_integrals,     &
          size(mo_one_e_integrals,1),                            &
          size(mo_one_e_integrals,2),                            &
          mo_label,1,.false.)
      call restore_symmetry(ao_num, mo_num, mo_coef, size(mo_coef,1), 1.d-10)
      SOFT_TOUCH mo_coef
    else if (mo_guess_type == "Huckel") then
      call huckel_guess
    else
      print *,  'Unrecognized MO guess type : '//mo_guess_type
      stop 1
    endif
    SOFT_TOUCH mo_label
  endif

end subroutine create_guess

! ---

subroutine routine_scf()

  implicit none
  integer                       :: i, j, it
  double precision              :: e_save, e_delta, rho_delta
  double precision, allocatable :: rho_old(:,:), rho_new(:,:)

  allocate(rho_old(ao_num,ao_num), rho_new(ao_num,ao_num))

  it = 0
  print*,'iteration = ', it

  !print*,'grad_hermit = ', grad_hermit
  print*,'***'
  print*,'TC HF total energy = ', TC_HF_energy
  print*,'TC HF 1 e   energy = ', TC_HF_one_electron_energy
  print*,'TC HF 2 e   energy = ', TC_HF_two_e_energy
  if(.not. bi_ortho)then
  print*,'TC HF 3 body       = ', diag_three_elem_hf
  endif
  print*,'***'
  e_delta = 10.d0
  e_save  = 0.d0 !TC_HF_energy
  rho_delta = 10.d0


  if(bi_ortho)then

   mo_l_coef = fock_tc_leigvec_ao
   mo_r_coef = fock_tc_reigvec_ao
   rho_old   = TCSCF_bi_ort_dm_ao
   call ezfio_set_bi_ortho_mos_mo_l_coef(mo_l_coef)
   call ezfio_set_bi_ortho_mos_mo_r_coef(mo_r_coef)
   TOUCH mo_l_coef mo_r_coef


  else

   print*,'grad_hermit = ',grad_hermit
   call save_good_hermit_tc_eigvectors
   TOUCH mo_coef 
   call save_mos

  endif

  ! ---

  if(bi_ortho) then

    !do while( it .lt. n_it_tcscf_max .and. (e_delta .gt. dsqrt(thresh_tcscf)) )
    !do while( it .lt. n_it_tcscf_max .and. (e_delta .gt. thresh_tcscf) )
    do while( it .lt. n_it_tcscf_max .and. (rho_delta .gt. thresh_tcscf) )

      it += 1
      print*,'iteration = ', it
      print*,'***'
      print*,'TC HF total energy = ', TC_HF_energy
      print*,'TC HF 1 e   energy = ', TC_HF_one_electron_energy
      print*,'TC HF 2 non hermit = ', TC_HF_two_e_energy
      print*,'***'
      e_delta = dabs( TC_HF_energy - e_save )
      print*, 'it, delta E = ', it, e_delta
      e_save    = TC_HF_energy
      mo_l_coef = fock_tc_leigvec_ao
      mo_r_coef = fock_tc_reigvec_ao

      rho_new   = TCSCF_bi_ort_dm_ao
      !print*, rho_new
      rho_delta = 0.d0
      do i = 1, ao_num 
        do j = 1, ao_num 
          rho_delta += dabs(rho_new(j,i) - rho_old(j,i))
        enddo
      enddo
      print*, ' rho_delta =', rho_delta
      rho_old = rho_new

      call ezfio_set_bi_ortho_mos_mo_l_coef(mo_l_coef)
      call ezfio_set_bi_ortho_mos_mo_r_coef(mo_r_coef)
      TOUCH mo_l_coef mo_r_coef

      call ezfio_set_tc_scf_bitc_energy(TC_HF_energy)

    enddo

  else
   do while( (grad_hermit.gt.dsqrt(thresh_tcscf)) .and. it .lt. n_it_tcscf_max )
      print*,'grad_hermit = ',grad_hermit
      it += 1
      print*,'iteration = ', it
      print*,'***'
      print*,'TC HF total energy = ', TC_HF_energy
      print*,'TC HF 1 e   energy = ', TC_HF_one_electron_energy
      print*,'TC HF 2 e   energy = ', TC_HF_two_e_energy
      print*,'TC HF 3 body       = ', diag_three_elem_hf
      print*,'***'
      call save_good_hermit_tc_eigvectors
      TOUCH mo_coef 
      call save_mos

    enddo

  endif

  print*,'Energy converged !'
  print*,'Diagonal Fock elements '
  do i = 1, mo_num
   print*,i,Fock_matrix_tc_mo_tot(i,i)
  enddo

  deallocate(rho_old, rho_new)

end subroutine routine_scf

! ---

