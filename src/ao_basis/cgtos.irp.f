
subroutine cgtos_not_implemented(feature)

  implicit none
  BEGIN_DOC
  ! Stops if use_cgtos is set: feature is computed with the real-Gaussian code
  ! (real parts of the exponents only, real normalization), which would silently
  ! give wrong results for cGTO AOs.
  END_DOC
  character*(*), intent(in) :: feature

  if (use_cgtos) then
    print *, ''
    print *, 'Error: '//trim(feature)//' is not implemented for cGTO AOs (ao_basis/use_cgtos = True).'
    print *, 'The real-Gaussian code would use the real parts of the exponents only.'
    stop 1
  endif

end

! ---

BEGIN_PROVIDER [logical, use_cgtos]

  implicit none

  BEGIN_DOC
  ! If true, use cgtos for AO integrals
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  use_cgtos = .False.
  if (mpi_master) then
    call ezfio_has_ao_basis_use_cgtos(has)
    if (has) then
!      write(6,'(A)') '.. >>>>> [ IO READ: use_cgtos ] <<<<< ..'
      call ezfio_get_ao_basis_use_cgtos(use_cgtos)
    else
      call ezfio_set_ao_basis_use_cgtos(use_cgtos)
    endif
  endif
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( use_cgtos, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read use_cgtos with MPI'
    endif
  IRP_ENDIF

!  call write_time(6)

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, ao_expo_im, (ao_num, ao_prim_num_max)]

  implicit none

  BEGIN_DOC
  ! Imaginary part $b$ of the exponents $a + i b$ of the cGTO primitives.
  ! Zero if absent from the EZFIO file.
  END_DOC

  logical                        :: has
  integer                        :: i, j
  PROVIDE ezfio_filename
  do j = 1, ao_prim_num_max
    do i = 1, ao_num
      ao_expo_im(i,j) = 0.d0
    enddo
  enddo
  if (mpi_master) then
    call ezfio_has_ao_basis_ao_expo_im(has)
    if (has) then
      call ezfio_get_ao_basis_ao_expo_im(ao_expo_im)
    endif
  endif
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( ao_expo_im, ao_num*ao_prim_num_max,               &
        MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read ao_expo_im with MPI'
    endif
  IRP_ENDIF

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, ao_expo_pw, (3, ao_num, ao_prim_num_max)]

  implicit none

  BEGIN_DOC
  ! Plane-wave vector $\mathbf{k}$ of the cGTO primitives.
  ! Zero (no plane wave) if absent from the EZFIO file.
  END_DOC

  logical                        :: has
  integer                        :: i, j, m
  PROVIDE ezfio_filename
  do j = 1, ao_prim_num_max
    do i = 1, ao_num
      do m = 1, 3
        ao_expo_pw(m,i,j) = 0.d0
      enddo
    enddo
  enddo
  if (mpi_master) then
    call ezfio_has_ao_basis_ao_expo_pw(has)
    if (has) then
      call ezfio_get_ao_basis_ao_expo_pw(ao_expo_pw)
    endif
  endif
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( ao_expo_pw, 3*ao_num*ao_prim_num_max,             &
        MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read ao_expo_pw with MPI'
    endif
  IRP_ENDIF

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, ao_expo_phase, (3, ao_num, ao_prim_num_max)]

  implicit none

  BEGIN_DOC
  ! Phase $\varphi$ (per direction) of the cGTO primitives.
  ! Zero if absent from the EZFIO file.
  END_DOC

  logical                        :: has
  integer                        :: i, j, m
  PROVIDE ezfio_filename
  do j = 1, ao_prim_num_max
    do i = 1, ao_num
      do m = 1, 3
        ao_expo_phase(m,i,j) = 0.d0
      enddo
    enddo
  enddo
  if (mpi_master) then
    call ezfio_has_ao_basis_ao_expo_phase(has)
    if (has) then
      call ezfio_get_ao_basis_ao_expo_phase(ao_expo_phase)
    endif
  endif
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( ao_expo_phase, 3*ao_num*ao_prim_num_max,          &
        MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read ao_expo_phase with MPI'
    endif
  IRP_ENDIF

END_PROVIDER

! ---

 BEGIN_PROVIDER [complex*16, ao_expo_cgtos_ord_transp, (ao_prim_num_max, ao_num)]
&BEGIN_PROVIDER [double precision, ao_expo_pw_ord_transp, (4, ao_prim_num_max, ao_num)]
&BEGIN_PROVIDER [double precision, ao_expo_phase_ord_transp, (4, ao_prim_num_max, ao_num)]

  implicit none
  BEGIN_DOC
  ! Transposed, sorted cGTO primitive data.
  ! ao_expo_cgtos_ord_transp(i,j) is the complex exponent of the
  ! i-th primitive of the j-th |AO|. Components 1 to 3 of
  ! ao_expo_pw_ord_transp are the plane wave part, component 4 its
  ! squared norm; ao_expo_phase_ord_transp uses the same layout
  ! for the phase and the total phase.
  END_DOC

  integer :: i, j, m

  do j = 1, ao_num
    do i = 1, ao_prim_num_max

      ao_expo_cgtos_ord_transp(i,j) = ao_expo_cgtos_ord(j,i)

      do m = 1, 4
        ao_expo_pw_ord_transp(m,i,j) = ao_expo_pw_ord(m,j,i)
        ao_expo_phase_ord_transp(m,i,j) = ao_expo_phase_ord(m,j,i)
      enddo
    enddo
  enddo

END_PROVIDER

! ---

 BEGIN_PROVIDER [double precision, ao_coef_norm_cgtos_ord, (ao_num, ao_prim_num_max)]
&BEGIN_PROVIDER [complex*16      , ao_expo_cgtos_ord, (ao_num, ao_prim_num_max)]
&BEGIN_PROVIDER [double precision, ao_expo_pw_ord, (4, ao_num, ao_prim_num_max)]
&BEGIN_PROVIDER [double precision, ao_expo_phase_ord, (4, ao_num, ao_prim_num_max)]

  implicit none
  BEGIN_DOC
  ! cGTO primitives of each |AO|, sorted by increasing real part
  ! of the exponent. Gives the complex exponent (ao_expo +
  ! i*ao_expo_im), the coefficients ao_coef_norm_cgtos, the plane
  ! wave part and the phase; component 4 of the last two is the
  ! squared norm of the plane wave part and the total phase.
  END_DOC

  integer          :: i, j, m
  integer          :: iorder(ao_prim_num_max)
  double precision :: d(ao_prim_num_max,11)

  ! entries j > ao_prim_num(i) are padding: set them to zero. They were left
  ! uninitialized, and use_pw takes the maxval of ao_expo_pw_ord_transp(4,:,:)
  ! over the whole array (a garbage value switches use_pw on)
  do j = 1, ao_prim_num_max
    do i = 1, ao_num
      ao_coef_norm_cgtos_ord(i,j) = 0.d0
      ao_expo_cgtos_ord(i,j)      = (0.d0, 0.d0)
      do m = 1, 4
        ao_expo_pw_ord(m,i,j)    = 0.d0
        ao_expo_phase_ord(m,i,j) = 0.d0
      enddo
    enddo
  enddo

  d = 0.d0

  do i = 1, ao_num

    do j = 1, ao_prim_num(i)
      iorder(j) = j
      d(j,1) = ao_expo(i,j)
      d(j,2) = ao_coef_norm_cgtos(i,j)
      d(j,3) = ao_expo_im(i,j)

      do m = 1, 3
        d(j,3+m) = ao_expo_pw(m,i,j)
      enddo
      d(j,7) = d(j,4) * d(j,4) + d(j,5) * d(j,5) + d(j,6) * d(j,6)

      do m = 1, 3
        d(j,7+m) = ao_expo_phase(m,i,j)
      enddo
      d(j,11) = d(j,8) + d(j,9) + d(j,10)
    enddo

    call dsort(d(1,1), iorder, ao_prim_num(i))
    do j = 2, 11
      call dset_order(d(1,j), iorder, ao_prim_num(i))
    enddo

    do j = 1, ao_prim_num(i)
      ao_expo_cgtos_ord     (i,j) = d(j,1) + (0.d0, 1.d0) * d(j,3)
      ao_coef_norm_cgtos_ord(i,j) = d(j,2)

      do m = 1, 4
        ao_expo_pw_ord(m,i,j) = d(j,3+m)
        ao_expo_phase_ord(m,i,j) = d(j,7+m)
      enddo
    enddo
  enddo

END_PROVIDER



! ---

BEGIN_PROVIDER [double precision, ao_coef_cgtos_norm_ord_transp, (ao_prim_num_max, ao_num)]

  implicit none
  BEGIN_DOC
  ! Transposed :c:data:`ao_coef_norm_cgtos_ord`
  END_DOC

  integer :: i, j

  do j = 1, ao_num
    do i = 1, ao_prim_num_max
      ao_coef_cgtos_norm_ord_transp(i,j) = ao_coef_norm_cgtos_ord(j,i)
    enddo
  enddo

END_PROVIDER


! ---

BEGIN_PROVIDER [double precision, ao_coef_norm_cgtos, (ao_num, ao_prim_num_max)]

  implicit none
  BEGIN_DOC
  ! Coefficients of the cGTO |AO|s, including the |AO|
  ! normalization, as :c:data:`ao_coef_normalized` does for
  ! real |AO|s. The self overlap of a primitive combines two
  ! overlap_cgaussian_xyz integrals, one with the complex Gaussian
  ! and one with its conjugate. The contracted |AO| is normalized
  ! when ao_normalized is set.
  END_DOC

  integer          :: i, j, k, ii, m, powA(3), kprim, nz
  double precision :: norm, norm_real, ov_rx, ov_ry, ov_rz
  double precision :: C_R(3)
  double precision :: kA2, phiA, KB2, phiB
  complex*16       :: expo, expo_inv, C_Ae(3), C_Ap(3)
  complex*16       :: expo_k, expo_k_inv, C_Be(3)
  complex*16       :: overlap_x, overlap_y, overlap_z
  complex*16       :: integ1, integ2, lC1(3), lC2(3)

  nz = 100
  C_R(1:3) = 0.d0

  ao_coef_norm_cgtos = 0.d0

  do i = 1, ao_num

    ii = ao_nucl(i)
    powA(1) = ao_power(i,1)
    powA(2) = ao_power(i,2)
    powA(3) = ao_power(i,3)
 
    if(primitives_normalized) then

      ! Normalization of the primitives. Each Cartesian component is normalized with its
      ! own powers (this defines the shape of a contracted cGTO). For ezfio_convention >=
      ! 20250211 the factor is scaled by prim_normalization_factor * N_real^-1 (real
      ! Gaussian, same powers), so that at b = 0 the AO is the one of ao_coef_normalized
      ! also when ao_normalized is false or the factors are read from the file (TREXIO);
      ! with the default factors this is a constant per component.
      if(ezfio_convention >= 20250211) then
        kprim = 1
        do while(kprim <= prim_num .and. shell_index(kprim) /= ao_shell(i))
          kprim = kprim + 1
        enddo
      endif

      do j = 1, ao_prim_num(i)

        expo = ao_expo(i,j) + (0.d0, 1.d0) * ao_expo_im(i,j)
        expo_inv = (1.d0, 0.d0) / expo
        do m = 1, 3
          C_Ap(m) = nucl_coord(ii,m)
          C_Ae(m) = nucl_coord(ii,m) - (0.d0, 0.5d0) * expo_inv * ao_expo_pw(m,i,j)
        enddo
        phiA = ao_expo_phase(1,i,j) + ao_expo_phase(2,i,j) + ao_expo_phase(3,i,j)
        KA2 = ao_expo_pw(1,i,j) * ao_expo_pw(1,i,j) &
            + ao_expo_pw(2,i,j) * ao_expo_pw(2,i,j) &
            + ao_expo_pw(3,i,j) * ao_expo_pw(3,i,j)

        ! logarithms of the plane-wave / phase prefactors of g g and g^* g, per direction,
        ! combined with the Gaussian-product factor before exponentiation (with plane
        ! waves both overflow / underflow separately when k^2 Re(1/zeta)/2 > 709)
        do m = 1, 3
          lC1(m) = -(0.d0, 2.d0) * ao_expo_phase(m,i,j) - 0.5d0 * expo_inv * ao_expo_pw(m,i,j)**2
          lC2(m) = -(0.5d0, 0.d0) * real(expo_inv) * ao_expo_pw(m,i,j)**2
        enddo

        call overlap_cgaussian_xyz_lpref(C_Ae, C_Ae, expo, expo, powA, powA, &
                                         C_Ap, C_Ap, lC1, overlap_x, overlap_y, overlap_z, integ1, nz)

        call overlap_cgaussian_xyz_lpref(conjg(C_Ae), C_Ae, conjg(expo), expo, powA, powA, &
                                         conjg(C_Ap), C_Ap, lC2, overlap_x, overlap_y, overlap_z, integ2, nz)

        norm = 2.d0 * real(integ1 + integ2)

        ! <g+g*|g+g*> = 2 Re(g|g) + 2 (g*|g) vanishes for a sin-type primitive (phase -pi/2)
        ! with b = 0, and loses about 2 log10(a/|b|) digits by cancellation when |b| << a
        if((norm <= 0.d0) .or. (norm /= norm)) then
          print *, 'Error: cGTO AO ', i, ' primitive ', j, ' has <g+g*|g+g*> = ', norm
          print *, '(a sin-type primitive, ao_expo_phase = -pi/2, needs ao_expo_im /= 0)'
          stop 1
        endif
        if(2.d0 * (abs(real(integ1)) + abs(real(integ2))) > 1.d6 * norm) then
          print *, 'Warning: cGTO AO ', i, ' primitive ', j, ': cancellation in <g+g*|g+g*>, ', &
                   'about ', int(log10(2.d0 * (abs(real(integ1)) + abs(real(integ2))) / norm)), &
                   ' digits lost in all its integrals (sin-type with |b| << a ?)'
        endif

        if(ezfio_convention >= 20250211) then
          call overlap_gaussian_xyz(C_R, C_R, ao_expo(i,j), ao_expo(i,j), powA, powA, &
                                    ov_rx, ov_ry, ov_rz, norm_real, nz)
          ao_coef_norm_cgtos(i,j) = ao_coef(i,j) * prim_normalization_factor(kprim+j-1) &
                                  * dsqrt(norm_real / norm)
        else
          ao_coef_norm_cgtos(i,j) = ao_coef(i,j) / dsqrt(norm)
        endif
      enddo

    else

      do j = 1, ao_prim_num(i)
        ao_coef_norm_cgtos(i,j) = ao_coef(i,j)
      enddo

    endif ! primitives_normalized

    if(ao_normalized) then

      ! Normalization of the contracted AO (as ao_coef_normalized does for real GTOs):
      ! <chi|chi> = sum_jk c_j c_k 2 Re[ C1 S(z_j, z_k) + C2 S(z_j^*, z_k) ]
      norm = 0.d0
      do j = 1, ao_prim_num(i)
        expo = ao_expo(i,j) + (0.d0, 1.d0) * ao_expo_im(i,j)
        expo_inv = (1.d0, 0.d0) / expo
        do m = 1, 3
          C_Ap(m) = nucl_coord(ii,m)
          C_Ae(m) = nucl_coord(ii,m) - (0.d0, 0.5d0) * expo_inv * ao_expo_pw(m,i,j)
        enddo
        phiA = ao_expo_phase(1,i,j) + ao_expo_phase(2,i,j) + ao_expo_phase(3,i,j)
        KA2 = ao_expo_pw(1,i,j) * ao_expo_pw(1,i,j) &
            + ao_expo_pw(2,i,j) * ao_expo_pw(2,i,j) &
            + ao_expo_pw(3,i,j) * ao_expo_pw(3,i,j)
        do k = 1, ao_prim_num(i)
          expo_k = ao_expo(i,k) + (0.d0, 1.d0) * ao_expo_im(i,k)
          expo_k_inv = (1.d0, 0.d0) / expo_k
          do m = 1, 3
            C_Be(m) = nucl_coord(ii,m) - (0.d0, 0.5d0) * expo_k_inv * ao_expo_pw(m,i,k)
          enddo
          phiB = ao_expo_phase(1,i,k) + ao_expo_phase(2,i,k) + ao_expo_phase(3,i,k)
          KB2 = ao_expo_pw(1,i,k) * ao_expo_pw(1,i,k) &
              + ao_expo_pw(2,i,k) * ao_expo_pw(2,i,k) &
              + ao_expo_pw(3,i,k) * ao_expo_pw(3,i,k)
          do m = 1, 3
            lC1(m) = -(0.d0, 1.d0) * (ao_expo_phase(m,i,j) + ao_expo_phase(m,i,k))   &
                   - 0.25d0 * (expo_inv * ao_expo_pw(m,i,j)**2 + expo_k_inv * ao_expo_pw(m,i,k)**2)
            lC2(m) =  (0.d0, 1.d0) * (ao_expo_phase(m,i,j) - ao_expo_phase(m,i,k))   &
                   - 0.25d0 * (conjg(expo_inv) * ao_expo_pw(m,i,j)**2 + expo_k_inv * ao_expo_pw(m,i,k)**2)
          enddo
          call overlap_cgaussian_xyz_lpref(C_Ae, C_Be, expo, expo_k, powA, powA, &
                                           C_Ap, C_Ap, lC1, overlap_x, overlap_y, overlap_z, integ1, nz)
          call overlap_cgaussian_xyz_lpref(conjg(C_Ae), C_Be, conjg(expo), expo_k, powA, powA, &
                                           conjg(C_Ap), C_Ap, lC2, overlap_x, overlap_y, overlap_z, integ2, nz)
          norm = norm + ao_coef_norm_cgtos(i,j) * ao_coef_norm_cgtos(i,k) &
                      * 2.d0 * real(integ1 + integ2)
        enddo
      enddo
      if((norm <= 0.d0) .or. (norm /= norm)) then
        print *, 'Error: contracted cGTO AO ', i, ' has <chi|chi> = ', norm
        stop 1
      endif
      do j = 1, ao_prim_num(i)
        ao_coef_norm_cgtos(i,j) = ao_coef_norm_cgtos(i,j) / dsqrt(norm)
      enddo

    endif ! ao_normalized

  enddo

END_PROVIDER


