
! ---

BEGIN_PROVIDER [double precision, j1e_val, (n_points_final_grid)]

  implicit none
  integer          :: ipoint, i, j, p
  double precision :: x, y, z, dx, dy, dz, d2
  double precision :: a, c, tmp
  double precision :: time0, time1

  PROVIDE j1e_type

  call wall_time(time0)
  print*, ' providing j1e_val ...'

  if(j1e_type .eq. "None") then

    j1e_val = 0.d0

  elseif(j1e_type .eq. "Gauss") then

    ! \sum_{A} \sum_p c_{p_A} \exp(-\alpha_{p_A} (r - R_A)^2)

    PROVIDE j1e_size j1e_coef j1e_expo

    do ipoint = 1, n_points_final_grid

      x = final_grid_points(1,ipoint)
      y = final_grid_points(2,ipoint)
      z = final_grid_points(3,ipoint)

      tmp = 0.d0
      do j = 1, nucl_num

        dx = x - nucl_coord(j,1)
        dy = y - nucl_coord(j,2)
        dz = z - nucl_coord(j,3)
        d2 = dx*dx + dy*dy + dz*dz

        do p = 1, j1e_size

          c = j1e_coef(p,j)
          a = j1e_expo(p,j)

          tmp = tmp + c * dexp(-a*d2)
        enddo
      enddo

      j1e_val(ipoint) = tmp
    enddo

  else

    print *, ' Error in j1e_val: Unknown j1e_type = ', j1e_type
    stop

  endif

  call wall_time(time1)
  print*, ' Wall time for j1e_val (min) = ', (time1 - time0) / 60.d0
  call print_memory_usage()

END_PROVIDER

! ---

 BEGIN_PROVIDER [double precision, j1e_gradx, (n_points_final_grid)]
&BEGIN_PROVIDER [double precision, j1e_grady, (n_points_final_grid)]
&BEGIN_PROVIDER [double precision, j1e_gradz, (n_points_final_grid)]

  implicit none
  integer                       :: ipoint, i, j, ij, p
  integer                       :: ierr
  logical                       :: exists
  double precision              :: x, y, z, dx, dy, dz, d2
  double precision              :: a, c, g, tmp_x, tmp_y, tmp_z
  double precision              :: cx, cy, cz
  double precision              :: time0, time1
  double precision, allocatable :: Pa(:,:), Pb(:,:), Pt(:,:)
  double precision, allocatable :: coef_fit(:), coef_fit2(:,:), coef_fit3(:,:)

  PROVIDE j1e_type

  call wall_time(time0)
  print*, ' providing j1e_grad ...'

  if(j1e_type .eq. "None") then

    j1e_gradx = 0.d0
    j1e_grady = 0.d0
    j1e_gradz = 0.d0

  elseif(j1e_type .eq. "Gauss") then

    ! - \sum_{A} (r - R_A) \sum_p c_{p_A} \exp(-\alpha_{p_A} (r - R_A)^2)

    PROVIDE j1e_size j1e_coef j1e_expo

    do ipoint = 1, n_points_final_grid

      x = final_grid_points(1,ipoint)
      y = final_grid_points(2,ipoint)
      z = final_grid_points(3,ipoint)

      tmp_x = 0.d0
      tmp_y = 0.d0
      tmp_z = 0.d0
      do j = 1, nucl_num

        dx = x - nucl_coord(j,1)
        dy = y - nucl_coord(j,2)
        dz = z - nucl_coord(j,3)
        d2 = dx*dx + dy*dy + dz*dz

        do p = 1, j1e_size

          c = j1e_coef(p,j)
          a = j1e_expo(p,j)
          g = c * a * dexp(-a*d2)

          tmp_x = tmp_x + g * dx
          tmp_y = tmp_y + g * dy
          tmp_z = tmp_z + g * dz
        enddo
      enddo

      j1e_gradx(ipoint) = -2.d0 * tmp_x
      j1e_grady(ipoint) = -2.d0 * tmp_y
      j1e_gradz(ipoint) = -2.d0 * tmp_z
    enddo

  elseif(j1e_type .eq. "Charge_Harmonizer") then

    ! -[(N-1)/2N] x \sum_{\mu,\nu} P_{\mu,\nu} \int dr2 [\grad_r1 J_2e(r1,r2)] \phi_\mu(r2) \phi_nu(r2) 

    PROVIDE elec_alpha_num elec_beta_num elec_num
    PROVIDE mo_coef
    PROVIDE int2_grad1_u2e_ao

    allocate(Pa(ao_num,ao_num), Pb(ao_num,ao_num), Pt(ao_num,ao_num))

    call dgemm( 'N', 'T', ao_num, ao_num, elec_alpha_num, 1.d0       &
              , mo_coef, size(mo_coef, 1), mo_coef, size(mo_coef, 1) &
              , 0.d0, Pa, size(Pa, 1))

    if(elec_alpha_num .eq. elec_beta_num) then
      Pb = Pa
    else
      call dgemm( 'N', 'T', ao_num, ao_num, elec_beta_num, 1.d0        &
                , mo_coef, size(mo_coef, 1), mo_coef, size(mo_coef, 1) &
                , 0.d0, Pb, size(Pb, 1))
    endif
    Pt = Pa + Pb

    g = -0.5d0 * (dble(elec_num) - 1.d0) / dble(elec_num)

    call dgemv("T", ao_num*ao_num, n_points_final_grid, g, int2_grad1_u2e_ao(1,1,1,1), ao_num*ao_num, Pt, 1, 0.d0, j1e_gradx, 1)
    call dgemv("T", ao_num*ao_num, n_points_final_grid, g, int2_grad1_u2e_ao(1,1,1,2), ao_num*ao_num, Pt, 1, 0.d0, j1e_grady, 1)
    call dgemv("T", ao_num*ao_num, n_points_final_grid, g, int2_grad1_u2e_ao(1,1,1,3), ao_num*ao_num, Pt, 1, 0.d0, j1e_gradz, 1)

    FREE int2_grad1_u2e_ao

    deallocate(Pa, Pb, Pt)

!  elseif(j1e_type .eq. "Charge_Harmonizer_AO") then
!
!    ! \grad_1 \sum_{\eta} C_{\eta} \chi_{\eta}
!    ! where 
!    !       \chi_{\eta} are the AOs
!    !       C_{\eta} are fitted to mimic (j1e_type .eq. "Charge_Harmonizer")
!    !
!    ! The - sign is in the parameters C_{\eta}
!
!    PROVIDE aos_grad_in_r_array
!
!    allocate(coef_fit(ao_num))
!
!    if(mpi_master) then
!      call ezfio_has_jastrow_j1e_coef_ao(exists)
!    endif 
!    IRP_IF MPI_DEBUG
!      print *,  irp_here, mpi_rank
!      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!    IRP_ENDIF
!    IRP_IF MPI
!      include 'mpif.h'
!      call MPI_BCAST(coef_fit, ao_num, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
!      if (ierr /= MPI_SUCCESS) then
!        stop 'Unable to read j1e_coef_ao with MPI'
!      endif
!    IRP_ENDIF
!    if(exists) then
!      if(mpi_master) then
!        write(6,'(A)') '.. >>>>> [ IO READ: j1e_coef_ao ] <<<<< ..'
!        call ezfio_get_jastrow_j1e_coef_ao(coef_fit)
!        IRP_IF MPI
!          call MPI_BCAST(coef_fit, ao_num, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
!          if (ierr /= MPI_SUCCESS) then
!            stop 'Unable to read j1e_coef_ao with MPI'
!          endif
!        IRP_ENDIF
!      endif
!    else
!
!      call get_j1e_coef_fit_ao(ao_num, coef_fit)
!      call ezfio_set_jastrow_j1e_coef_ao(coef_fit)
!
!    endif
!
!    !$OMP PARALLEL                               &
!    !$OMP DEFAULT (NONE)                         &
!    !$OMP PRIVATE (i, ipoint, c)                 &
!    !$OMP SHARED (n_points_final_grid, ao_num,   &
!    !$OMP         aos_grad_in_r_array, coef_fit, &
!    !$OMP         j1e_gradx, j1e_grady, j1e_gradz)
!    !$OMP DO SCHEDULE (static)
!    do ipoint = 1, n_points_final_grid
!
!      j1e_gradx(ipoint) = 0.d0
!      j1e_grady(ipoint) = 0.d0
!      j1e_gradz(ipoint) = 0.d0
!      do i = 1, ao_num
!        c = coef_fit(i)
!        j1e_gradx(ipoint) = j1e_gradx(ipoint) + c * aos_grad_in_r_array(i,ipoint,1)
!        j1e_grady(ipoint) = j1e_grady(ipoint) + c * aos_grad_in_r_array(i,ipoint,2)
!        j1e_gradz(ipoint) = j1e_gradz(ipoint) + c * aos_grad_in_r_array(i,ipoint,3)
!      enddo
!    enddo
!    !$OMP END DO
!    !$OMP END PARALLEL
!
!    deallocate(coef_fit)

  elseif(j1e_type .eq. "Charge_Harmonizer_AO") then

    ! \grad_1 \sum_{\eta,\beta} C_{\eta,\beta} \chi_{\eta} \chi_{\beta}
    ! where 
    !       \chi_{\eta} are the AOs
    !       C_{\eta,\beta} are fitted to mimic (j1e_type .eq. "Charge_Harmonizer")
    !
    ! The - sign is in the parameters C_{\eta,\beta}

    PROVIDE aos_grad_in_r_array

    allocate(coef_fit2(ao_num,ao_num))

    if(mpi_master) then
      call ezfio_has_jastrow_j1e_coef_ao2(exists)
    endif 
    IRP_IF MPI_DEBUG
      print *,  irp_here, mpi_rank
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    IRP_ENDIF
    IRP_IF MPI
      include 'mpif.h'
      call MPI_BCAST(coef_fit2, (ao_num*ao_num), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) then
        stop 'Unable to read j1e_coef_ao2 with MPI'
      endif
    IRP_ENDIF
    if(exists) then
      if(mpi_master) then
        write(6,'(A)') '.. >>>>> [ IO READ: j1e_coef_ao2 ] <<<<< ..'
        call ezfio_get_jastrow_j1e_coef_ao2(coef_fit2)
        IRP_IF MPI
          call MPI_BCAST(coef_fit2, (ao_num*ao_num), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
          if (ierr /= MPI_SUCCESS) then
            stop 'Unable to read j1e_coef_ao2 with MPI'
          endif
        IRP_ENDIF
      endif
    else

      call get_j1e_coef_fit_ao2(ao_num, coef_fit2)
      call ezfio_set_jastrow_j1e_coef_ao2(coef_fit2)

    endif

    !$OMP PARALLEL                                &
    !$OMP DEFAULT (NONE)                          &
    !$OMP PRIVATE (i, j, ipoint, c)               &
    !$OMP SHARED (n_points_final_grid, ao_num,    &
    !$OMP         aos_grad_in_r_array, coef_fit2, &
    !$OMP         aos_in_r_array, j1e_gradx, j1e_grady, j1e_gradz)
    !$OMP DO SCHEDULE (static)
    do ipoint = 1, n_points_final_grid

      j1e_gradx(ipoint) = 0.d0
      j1e_grady(ipoint) = 0.d0
      j1e_gradz(ipoint) = 0.d0

      do i = 1, ao_num
        do j = 1, ao_num
          c = coef_fit2(j,i)

          j1e_gradx(ipoint) += c * (aos_in_r_array(i,ipoint) * aos_grad_in_r_array(j,ipoint,1) + aos_grad_in_r_array(i,ipoint,1) * aos_in_r_array(j,ipoint))
          j1e_grady(ipoint) += c * (aos_in_r_array(i,ipoint) * aos_grad_in_r_array(j,ipoint,2) + aos_grad_in_r_array(i,ipoint,2) * aos_in_r_array(j,ipoint))
          j1e_gradz(ipoint) += c * (aos_in_r_array(i,ipoint) * aos_grad_in_r_array(j,ipoint,3) + aos_grad_in_r_array(i,ipoint,3) * aos_in_r_array(j,ipoint))
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    deallocate(coef_fit2)

!  elseif(j1e_type .eq. "Charge_Harmonizer_AO3") then
!
!    ! \sum_{\eta} \vec{C}_{\eta} \chi_{\eta}
!    ! where 
!    !       \chi_{\eta} are the AOs
!    !       \vec{C}_{\eta} are fitted to mimic (j1e_type .eq. "Charge_Harmonizer")
!    !
!    ! The - sign is in the parameters \vec{C}_{\eta}
!
!    PROVIDE aos_grad_in_r_array
!
!    allocate(coef_fit3(ao_num,3))
!
!    if(mpi_master) then
!      call ezfio_has_jastrow_j1e_coef_ao3(exists)
!    endif 
!    IRP_IF MPI_DEBUG
!      print *,  irp_here, mpi_rank
!      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!    IRP_ENDIF
!    IRP_IF MPI
!      !include 'mpif.h'
!      call MPI_BCAST(coef_fit3, (ao_num*3), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
!      if (ierr /= MPI_SUCCESS) then
!        stop 'Unable to read j1e_coef_ao3 with MPI'
!      endif
!    IRP_ENDIF
!    if(exists) then
!      if(mpi_master) then
!        write(6,'(A)') '.. >>>>> [ IO READ: j1e_coef_ao3 ] <<<<< ..'
!        call ezfio_get_jastrow_j1e_coef_ao3(coef_fit3)
!        IRP_IF MPI
!          call MPI_BCAST(coef_fit3, (ao_num*3), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
!          if (ierr /= MPI_SUCCESS) then
!            stop 'Unable to read j1e_coef_ao3 with MPI'
!          endif
!        IRP_ENDIF
!      endif
!    else
!
!      call get_j1e_coef_fit_ao3(ao_num, coef_fit3)
!      call ezfio_set_jastrow_j1e_coef_ao3(coef_fit3)
!
!    endif
!
!    !$OMP PARALLEL                                &
!    !$OMP DEFAULT (NONE)                          &
!    !$OMP PRIVATE (i, ipoint, cx, cy, cz)         &
!    !$OMP SHARED (n_points_final_grid, ao_num,    &
!    !$OMP         aos_grad_in_r_array, coef_fit3, &
!    !$OMP         aos_in_r_array, j1e_gradx, j1e_grady, j1e_gradz)
!    !$OMP DO SCHEDULE (static)
!    do ipoint = 1, n_points_final_grid
!
!      j1e_gradx(ipoint) = 0.d0
!      j1e_grady(ipoint) = 0.d0
!      j1e_gradz(ipoint) = 0.d0
!      do i = 1, ao_num
!        cx = coef_fit3(i,1)
!        cy = coef_fit3(i,2)
!        cz = coef_fit3(i,3)
!
!        j1e_gradx(ipoint) += cx * aos_in_r_array(i,ipoint)
!        j1e_grady(ipoint) += cy * aos_in_r_array(i,ipoint)
!        j1e_gradz(ipoint) += cz * aos_in_r_array(i,ipoint)
!      enddo
!    enddo
!    !$OMP END DO
!    !$OMP END PARALLEL
!
!    deallocate(coef_fit3)

  else

    print *, ' Error in j1e_grad: Unknown j1e_type = ', j1e_type
    stop

  endif

  call wall_time(time1)
  print*, ' Wall time for j1e_grad (min) = ', (time1 - time0) / 60.d0
  call print_memory_usage()

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, j1e_lapl, (n_points_final_grid)]

  implicit none
  integer          :: ipoint, i, j, p
  double precision :: x, y, z, dx, dy, dz, d2
  double precision :: a, c, g, tmp

  if(j1e_type .eq. "None") then

    j1e_lapl = 0.d0

  elseif(j1e_type .eq. "Gauss") then

    ! - \sum_{A} (r - R_A) \sum_p c_{p_A} \exp(-\alpha_{p_A} (r - R_A)^2)

    PROVIDE j1e_size j1e_coef j1e_expo

    do ipoint = 1, n_points_final_grid

      x = final_grid_points(1,ipoint)
      y = final_grid_points(2,ipoint)
      z = final_grid_points(3,ipoint)

      tmp = 0.d0
      do j = 1, nucl_num

        dx = x - nucl_coord(j,1)
        dy = y - nucl_coord(j,2)
        dz = z - nucl_coord(j,3)
        d2 = dx*dx + dy*dy + dz*dz

        do p = 1, j1e_size

          c = j1e_coef(p,j)
          a = j1e_expo(p,j)
          g = c * a * dexp(-a*d2)

          tmp = tmp + (2.d0 * a * d2 - 3.d0) * g
        enddo
      enddo

      j1e_lapl(ipoint) = tmp
    enddo

  else

    print *, ' Error in j1e_lapl: Unknown j1e_type = ', j1e_type
    stop

  endif

END_PROVIDER

! ---

