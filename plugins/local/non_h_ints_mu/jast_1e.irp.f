
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
  integer                       :: ipoint, i, j, p
  double precision              :: x, y, z, dx, dy, dz, d2
  double precision              :: a, c, g, tmp_x, tmp_y, tmp_z
  double precision              :: time0, time1
  double precision, allocatable :: Pa(:,:), Pb(:,:), Pt(:,:)
  double precision, allocatable :: coef_fit(:)

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

          tmp_x = tmp_x - g * dx
          tmp_y = tmp_y - g * dy
          tmp_z = tmp_z - g * dz
        enddo
      enddo

      j1e_gradx(ipoint) = 2.d0 * tmp_x
      j1e_grady(ipoint) = 2.d0 * tmp_y
      j1e_gradz(ipoint) = 2.d0 * tmp_z
    enddo

  elseif(j1e_type .eq. "Charge_Harmonizer") then

    ! The - sign is in the integral over r2
    ! [(N-1)/2N] x \sum_{\mu,\nu} P_{\mu,\nu} \int dr2 [-1 * \grad_r1 J(r1,r2)] \phi_\mu(r2) \phi_nu(r2) 

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

    g = 0.5d0 * (dble(elec_num) - 1.d0) / dble(elec_num)

    call dgemv("T", ao_num*ao_num, n_points_final_grid, g, int2_grad1_u2e_ao(1,1,1,1), ao_num*ao_num, Pt, 1, 0.d0, j1e_gradx, 1)
    call dgemv("T", ao_num*ao_num, n_points_final_grid, g, int2_grad1_u2e_ao(1,1,1,2), ao_num*ao_num, Pt, 1, 0.d0, j1e_grady, 1)
    call dgemv("T", ao_num*ao_num, n_points_final_grid, g, int2_grad1_u2e_ao(1,1,1,3), ao_num*ao_num, Pt, 1, 0.d0, j1e_gradz, 1)

    FREE int2_grad1_u2e_ao

    deallocate(Pa, Pb, Pt)

  elseif(j1e_type .eq. "Charge_Harmonizer_AO") then

    ! \grad_1 \sum_{\eta} C_{\eta} \chi_{\eta}
    ! where 
    !       \chi_{\eta} are the AOs
    !       C_{\eta} are fitted to mimic (j1e_type .eq. "Charge_Harmonizer")
    !
    ! The - sign is in the parameters C_{\eta}

    PROVIDE aos_grad_in_r_array

    allocate(coef_fit(ao_num))

    call get_j1e_coef_fit_ao(ao_num, coef_fit)
    call ezfio_set_jastrow_j1e_coef_ao(coef_fit)

    !$OMP PARALLEL                                 &
    !$OMP DEFAULT (NONE)                           &
    !$OMP PRIVATE (i, ipoint, tmp_x, tmp_y, tmp_z, &
    !$OMP         c)                               &
    !$OMP SHARED (n_points_final_grid, ao_num,     &
    !$OMP         aos_grad_in_r_array, coef_fit,   &
    !$OMP         j1e_gradx, j1e_grady, j1e_gradz)
    !$OMP DO SCHEDULE (static)
    do ipoint = 1, n_points_final_grid

      tmp_x = 0.d0
      tmp_y = 0.d0
      tmp_z = 0.d0
      do i = 1, ao_num
        c = coef_fit(i)
        tmp_x = tmp_x + c * aos_grad_in_r_array(i,ipoint,1)
        tmp_y = tmp_y + c * aos_grad_in_r_array(i,ipoint,2)
        tmp_z = tmp_z + c * aos_grad_in_r_array(i,ipoint,3)
      enddo

      j1e_gradx(ipoint) = tmp_x
      j1e_grady(ipoint) = tmp_y
      j1e_gradz(ipoint) = tmp_z
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    deallocate(coef_fit)

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

