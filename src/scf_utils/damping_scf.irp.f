subroutine damping_SCF
  implicit none
  double precision               :: E
  double precision, allocatable  :: D_alpha(:,:), D_beta(:,:)
  double precision               :: E_new
  double precision, allocatable  :: D_new_alpha(:,:), D_new_beta(:,:), F_new(:,:)
  double precision, allocatable  :: delta_alpha(:,:), delta_beta(:,:)
  double precision               :: lambda, E_half, a, b, delta_D, delta_E, E_min

  integer                        :: i,j,k
  logical                        :: saving
  character                      :: save_char

  allocate(                                                          &
      D_alpha( ao_num, ao_num ),                               &
      D_beta( ao_num, ao_num ),                                &
      F_new( ao_num, ao_num ),                                 &
      D_new_alpha( ao_num, ao_num ),                           &
      D_new_beta( ao_num, ao_num ),                            &
      delta_alpha( ao_num, ao_num ),                           &
      delta_beta( ao_num, ao_num ))

  do j=1,ao_num
    do i=1,ao_num
      D_alpha(i,j) = SCF_density_matrix_ao_alpha(i,j)
      D_beta (i,j) = SCF_density_matrix_ao_beta (i,j)
    enddo
  enddo


  call write_time(6)

  write(6,'(A4,1X,A16, 1X, A16, 1X, A16, 1X, A4 )')  &
    '====','================','================','================', '===='
  write(6,'(A4,1X,A16, 1X, A16, 1X, A16, 1X, A4 )')  &
    '  N ', 'Energy  ', 'Energy diff  ', 'Density diff  ', 'Save'
  write(6,'(A4,1X,A16, 1X, A16, 1X, A16, 1X, A4 )')  &
    '====','================','================','================', '===='

  E = SCF_energy + 1.d0
  E_min = SCF_energy
  delta_D = 0.d0
  do k=1,n_it_scf_max

    delta_E = SCF_energy - E
    E = SCF_energy

    if ( (delta_E < 0.d0).and.(dabs(delta_E) < thresh_scf) ) then
      exit
    endif

    saving = E < E_min
    if (saving) then
      call save_mos
      save_char = 'X'
      E_min = E
    else
      save_char = ' '
    endif

    write(6,'(I4,1X,F16.10, 1X, F16.10, 1X, F16.10, 3X, A )')  &
      k, E, delta_E, delta_D, save_char

    if(frozen_orb_scf)then
     call initialize_mo_coef_begin_iteration
    endif
    D_alpha = SCF_density_matrix_ao_alpha
    D_beta  = SCF_density_matrix_ao_beta
    mo_coef = eigenvectors_fock_matrix_mo
    if(frozen_orb_scf)then
     call reorder_core_orb
     call initialize_mo_coef_begin_iteration
    endif
    TOUCH mo_coef

    D_new_alpha = SCF_density_matrix_ao_alpha
    D_new_beta  = SCF_density_matrix_ao_beta
    F_new = Fock_matrix_ao
    E_new = SCF_energy

    delta_alpha = D_new_alpha - D_alpha
    delta_beta  = D_new_beta  - D_beta

    lambda = .5d0
    E_half = 0.d0
    do while (E_half > E)
      SCF_density_matrix_ao_alpha = D_alpha + lambda * delta_alpha
      SCF_density_matrix_ao_beta  = D_beta  + lambda * delta_beta
      TOUCH SCF_density_matrix_ao_alpha SCF_density_matrix_ao_beta
      mo_coef = eigenvectors_fock_matrix_mo
      if(frozen_orb_scf)then
       call reorder_core_orb
       call initialize_mo_coef_begin_iteration
      endif
      TOUCH mo_coef
      E_half = SCF_energy
      if ((E_half > E).and.(E_new < E)) then
        lambda = 1.d0
        exit
      else if ((E_half > E).and.(lambda > 5.d-4)) then
        lambda = 0.5d0 * lambda
        E_new = E_half
      else
        exit
      endif
    enddo

    a = (E_new + E - 2.d0*E_half)*2.d0
    b = -E_new - 3.d0*E + 4.d0*E_half
    lambda = -lambda*b/(a+1.d-16)
    D_alpha = (1.d0-lambda) * D_alpha + lambda * D_new_alpha
    D_beta  = (1.d0-lambda) * D_beta  + lambda * D_new_beta
    delta_E = SCF_energy - E
    do j=1,ao_num
      do i=1,ao_num
        delta_D = delta_D + &
        (D_alpha(i,j) - SCF_density_matrix_ao_alpha(i,j))*(D_alpha(i,j) - SCF_density_matrix_ao_alpha(i,j)) + &
        (D_beta (i,j) - SCF_density_matrix_ao_beta (i,j))*(D_beta (i,j) - SCF_density_matrix_ao_beta (i,j))
      enddo
    enddo
    delta_D = dsqrt(delta_D/dble(ao_num)**2)
    SCF_density_matrix_ao_alpha = D_alpha
    SCF_density_matrix_ao_beta  = D_beta
    TOUCH SCF_density_matrix_ao_alpha SCF_density_matrix_ao_beta
    mo_coef = eigenvectors_fock_matrix_mo
    if(frozen_orb_scf)then
     call reorder_core_orb
     call initialize_mo_coef_begin_iteration
    endif
    TOUCH mo_coef

  enddo
  write(6,'(A4,1X,A16, 1X, A16, 1X, A16, 1X, A4 )')  '====','================','================','================', '===='
  write(6,*)

  if(.not.frozen_orb_scf)then
   call mo_as_eigvectors_of_mo_matrix(Fock_matrix_mo,size(Fock_matrix_mo,1),size(Fock_matrix_mo,2),mo_label,1,.true.)
  endif

  call write_double(6, E_min, 'Hartree-Fock energy')
  call ezfio_set_hartree_fock_energy(E_min)

  call write_time(6)

  deallocate(D_alpha,D_beta,F_new,D_new_alpha,D_new_beta,delta_alpha,delta_beta)
end
