
! ---

 BEGIN_PROVIDER [ double precision, ao_deriv2_x, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, ao_deriv2_y, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, ao_deriv2_z, (ao_num, ao_num) ]

  BEGIN_DOC
  ! Second derivative matrix elements in the |AO| basis.
  !
  ! .. math::
  !
  !   {\tt ao\_deriv2\_x} =
  !   \langle \chi_i(x,y,z) | \frac{\partial^2}{\partial x^2} |\chi_j (x,y,z) \rangle
  !
  END_DOC

  implicit none
  integer          :: i, j, n, l, dim1, power_A(3), power_B(3)
  double precision :: overlap, overlap_y, overlap_z
  double precision :: overlap_x0, overlap_y0, overlap_z0
  double precision :: alpha, beta, c
  double precision :: A_center(3), B_center(3)
  double precision :: d_a_2,d_2

  if(use_cosgtos) then
    !print*, 'use_cosgtos for ao_kinetic_integrals ?', use_cosgtos

    do j = 1, ao_num
      do i = 1, ao_num
        ao_deriv2_x(i,j) = ao_deriv2_cosgtos_x(i,j)
        ao_deriv2_y(i,j) = ao_deriv2_cosgtos_y(i,j)
        ao_deriv2_z(i,j) = ao_deriv2_cosgtos_z(i,j)
      enddo
    enddo

  else

    dim1=100

    ! -- Dummy call to provide everything
    A_center(:) = 0.d0
    B_center(:) = 1.d0
    alpha = 1.d0
    beta  = .1d0
    power_A = 1
    power_B = 0
    call overlap_gaussian_xyz(A_center,B_center,alpha,beta,power_A,power_B,overlap_y,d_a_2,overlap_z,overlap,dim1)
    ! --

    !$OMP PARALLEL DO SCHEDULE(GUIDED) &
    !$OMP DEFAULT(NONE) &
    !$OMP PRIVATE(A_center,B_center,power_A,power_B,&
    !$OMP  overlap_y, overlap_z, overlap, &
    !$OMP  alpha, beta,i,j,c,d_a_2,d_2,deriv_tmp, &
    !$OMP  overlap_x0,overlap_y0,overlap_z0) &
    !$OMP SHARED(nucl_coord,ao_power,ao_prim_num, &
    !$OMP  ao_deriv2_x,ao_deriv2_y,ao_deriv2_z,ao_num,ao_coef_normalized_ordered_transp,ao_nucl, &
    !$OMP  ao_expo_ordered_transp,dim1)
    do j=1,ao_num
     A_center(1) = nucl_coord( ao_nucl(j), 1 )
     A_center(2) = nucl_coord( ao_nucl(j), 2 )
     A_center(3) = nucl_coord( ao_nucl(j), 3 )
     power_A(1)  = ao_power( j, 1 )
     power_A(2)  = ao_power( j, 2 )
     power_A(3)  = ao_power( j, 3 )
     do i= 1,ao_num
      ao_deriv2_x(i,j)= 0.d0
      ao_deriv2_y(i,j)= 0.d0
      ao_deriv2_z(i,j)= 0.d0
      B_center(1) = nucl_coord( ao_nucl(i), 1 )
      B_center(2) = nucl_coord( ao_nucl(i), 2 )
      B_center(3) = nucl_coord( ao_nucl(i), 3 )
      power_B(1)  = ao_power( i, 1 )
      power_B(2)  = ao_power( i, 2 )
      power_B(3)  = ao_power( i, 3 )
      do n = 1,ao_prim_num(j)
       alpha = ao_expo_ordered_transp(n,j)
       do l = 1, ao_prim_num(i)
        beta = ao_expo_ordered_transp(l,i)
        call overlap_gaussian_xyz(A_center,B_center,alpha,beta,power_A,power_B,overlap_x0,overlap_y0,overlap_z0,overlap,dim1)
        c = ao_coef_normalized_ordered_transp(n,j) * ao_coef_normalized_ordered_transp(l,i)

        power_A(1) = power_A(1)-2
        if (power_A(1)>-1) then
          call overlap_gaussian_xyz(A_center,B_center,alpha,beta,power_A,power_B,d_a_2,overlap_y,overlap_z,overlap,dim1)
        else
          d_a_2 = 0.d0
        endif
        power_A(1) = power_A(1)+4
        call overlap_gaussian_xyz(A_center,B_center,alpha,beta,power_A,power_B,d_2,overlap_y,overlap_z,overlap,dim1)
        power_A(1) = power_A(1)-2

        double precision :: deriv_tmp
        deriv_tmp = (-2.d0 * alpha * (2.d0 * power_A(1) +1.d0) * overlap_x0 &
        +power_A(1) * (power_A(1)-1.d0) * d_a_2 &
        +4.d0 * alpha * alpha * d_2   )*overlap_y0*overlap_z0

        ao_deriv2_x(i,j) += c*deriv_tmp
        power_A(2) = power_A(2)-2
        if (power_A(2)>-1) then
          call overlap_gaussian_xyz(A_center,B_center,alpha,beta,power_A,power_B,overlap_y,d_a_2,overlap_z,overlap,dim1)
        else
          d_a_2 = 0.d0
        endif
        power_A(2) = power_A(2)+4
        call overlap_gaussian_xyz(A_center,B_center,alpha,beta,power_A,power_B,overlap_y,d_2,overlap_z,overlap,dim1)
        power_A(2) = power_A(2)-2

        deriv_tmp = (-2.d0 * alpha * (2.d0 * power_A(2) +1.d0 ) * overlap_y0 &
        +power_A(2) * (power_A(2)-1.d0) * d_a_2 &
        +4.d0 * alpha * alpha * d_2   )*overlap_x0*overlap_z0
        ao_deriv2_y(i,j) += c*deriv_tmp

        power_A(3) = power_A(3)-2
        if (power_A(3)>-1) then
          call overlap_gaussian_xyz(A_center,B_center,alpha,beta,power_A,power_B,overlap_y,overlap_z,d_a_2,overlap,dim1)
        else
          d_a_2 = 0.d0
        endif
        power_A(3) = power_A(3)+4
        call overlap_gaussian_xyz(A_center,B_center,alpha,beta,power_A,power_B,overlap_y,overlap_z,d_2,overlap,dim1)
        power_A(3) = power_A(3)-2

        deriv_tmp = (-2.d0 * alpha * (2.d0 * power_A(3) +1.d0 ) * overlap_z0 &
        +power_A(3) * (power_A(3)-1.d0) * d_a_2 &
        +4.d0 * alpha * alpha * d_2   )*overlap_x0*overlap_y0
        ao_deriv2_z(i,j) += c*deriv_tmp

       enddo
      enddo
     enddo
    enddo
    !$OMP END PARALLEL DO

  endif

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, ao_kinetic_integrals, (ao_num,ao_num)]
  implicit none
  BEGIN_DOC
  ! Kinetic energy integrals in the |AO| basis.
  !
  ! $\langle \chi_i |\hat{T}| \chi_j \rangle$
  !
  END_DOC
  integer                        :: i,j,k,l

  if (read_ao_integrals_kinetic) then
    call ezfio_get_ao_one_e_ints_ao_integrals_kinetic(ao_kinetic_integrals)
    print *,  'AO kinetic integrals read from disk'
  else
    !$OMP PARALLEL DO DEFAULT(NONE) &
    !$OMP  PRIVATE(i,j) &
    !$OMP  SHARED(ao_num, ao_kinetic_integrals,ao_deriv2_x,ao_deriv2_y,ao_deriv2_z)
    do j = 1, ao_num
      do i = 1, ao_num
      ao_kinetic_integrals(i,j) = -0.5d0 * (ao_deriv2_x(i,j) + ao_deriv2_y(i,j) + ao_deriv2_z(i,j) )
      enddo
    enddo
    !$OMP END PARALLEL DO
  endif
  if (write_ao_integrals_kinetic) then
    call ezfio_set_ao_one_e_ints_ao_integrals_kinetic(ao_kinetic_integrals)
    print *,  'AO kinetic integrals written to disk'
  endif
END_PROVIDER

BEGIN_PROVIDER [double precision, ao_kinetic_integrals_imag, (ao_num,ao_num)]
  implicit none
  BEGIN_DOC
  ! Kinetic energy integrals in the |AO| basis.
  !
  ! $\langle \chi_i |\hat{T}| \chi_j \rangle$
  !
  END_DOC
  integer                        :: i,j,k,l

  if (read_ao_integrals_kinetic) then
    call ezfio_get_ao_one_e_ints_ao_integrals_kinetic(ao_kinetic_integrals_imag)
    print *,  'AO kinetic integrals read from disk'
  else
    print *,  irp_here, ': Not yet implemented'
  endif
  if (write_ao_integrals_kinetic) then
    call ezfio_set_ao_one_e_ints_ao_integrals_kinetic(ao_kinetic_integrals_imag)
    print *,  'AO kinetic integrals written to disk'
  endif
END_PROVIDER

