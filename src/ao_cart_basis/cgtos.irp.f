
BEGIN_PROVIDER [logical, use_cgtos]

  implicit none

  BEGIN_DOC
  ! If true, use cgtos for AO integrals
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  use_cgtos = .False.
  if (mpi_master) then
    call ezfio_has_ao_cart_basis_use_cgtos(has)
    if (has) then
!      write(6,'(A)') '.. >>>>> [ IO READ: use_cgtos ] <<<<< ..'
      call ezfio_get_ao_cart_basis_use_cgtos(use_cgtos)
    else
      call ezfio_set_ao_cart_basis_use_cgtos(use_cgtos)
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

 BEGIN_PROVIDER [complex*16, ao_cart_expo_cgtos_ord_transp, (ao_cart_prim_num_max, ao_cart_num)]
&BEGIN_PROVIDER [double precision, ao_cart_expo_pw_ord_transp, (4, ao_cart_prim_num_max, ao_cart_num)]
&BEGIN_PROVIDER [double precision, ao_cart_expo_phase_ord_transp, (4, ao_cart_prim_num_max, ao_cart_num)]

  implicit none

  integer :: i, j, m

  do j = 1, ao_cart_num
    do i = 1, ao_cart_prim_num_max

      ao_cart_expo_cgtos_ord_transp(i,j) = ao_cart_expo_cgtos_ord(j,i)

      do m = 1, 4
        ao_cart_expo_pw_ord_transp(m,i,j) = ao_cart_expo_pw_ord(m,j,i)
        ao_cart_expo_phase_ord_transp(m,i,j) = ao_cart_expo_phase_ord(m,j,i)
      enddo
    enddo
  enddo

END_PROVIDER

! ---

 BEGIN_PROVIDER [double precision, ao_cart_coef_norm_cgtos_ord, (ao_cart_num, ao_cart_prim_num_max)]
&BEGIN_PROVIDER [complex*16      , ao_cart_expo_cgtos_ord, (ao_cart_num, ao_cart_prim_num_max)]
&BEGIN_PROVIDER [double precision, ao_cart_expo_pw_ord, (4, ao_cart_num, ao_cart_prim_num_max)]
&BEGIN_PROVIDER [double precision, ao_cart_expo_phase_ord, (4, ao_cart_num, ao_cart_prim_num_max)]

  implicit none

  integer          :: i, j, m
  integer          :: iorder(ao_cart_prim_num_max)
  double precision :: d(ao_cart_prim_num_max,11)

  d = 0.d0

  do i = 1, ao_cart_num

    do j = 1, ao_cart_prim_num(i)
      iorder(j) = j
      d(j,1) = ao_cart_expo(i,j)
      d(j,2) = ao_cart_coef_norm_cgtos(i,j)
      d(j,3) = ao_cart_expo_im(i,j)

      do m = 1, 3
        d(j,3+m) = ao_cart_expo_pw(m,i,j)
      enddo
      d(j,7) = d(j,4) * d(j,4) + d(j,5) * d(j,5) + d(j,6) * d(j,6)

      do m = 1, 3
        d(j,7+m) = ao_cart_expo_phase(m,i,j)
      enddo
      d(j,11) = d(j,8) + d(j,9) + d(j,10)
    enddo

    call dsort(d(1,1), iorder, ao_cart_prim_num(i))
    do j = 2, 11
      call dset_order(d(1,j), iorder, ao_cart_prim_num(i))
    enddo

    do j = 1, ao_cart_prim_num(i)
      ao_cart_expo_cgtos_ord     (i,j) = d(j,1) + (0.d0, 1.d0) * d(j,3)
      ao_cart_coef_norm_cgtos_ord(i,j) = d(j,2)

      do m = 1, 4
        ao_cart_expo_pw_ord(m,i,j) = d(j,3+m)
        ao_cart_expo_phase_ord(m,i,j) = d(j,7+m)
      enddo
    enddo
  enddo

END_PROVIDER



! ---

BEGIN_PROVIDER [double precision, ao_cart_coef_cgtos_norm_ord_transp, (ao_cart_prim_num_max, ao_cart_num)]

  implicit none

  integer :: i, j

  do j = 1, ao_cart_num
    do i = 1, ao_cart_prim_num_max
      ao_cart_coef_cgtos_norm_ord_transp(i,j) = ao_cart_coef_norm_cgtos_ord(j,i)
    enddo
  enddo

END_PROVIDER


! ---

BEGIN_PROVIDER [double precision, ao_cart_coef_norm_cgtos, (ao_cart_num, ao_cart_prim_num_max)]

  implicit none

  integer          :: i, j, ii, m, powA(3), nz
  double precision :: norm
  double precision :: kA2, phiA
  complex*16       :: expo, expo_inv, C_Ae(3), C_Ap(3)
  complex*16       :: overlap_x, overlap_y, overlap_z
  complex*16       :: integ1, integ2, C1, C2

  nz = 100

  ao_cart_coef_norm_cgtos = 0.d0

  do i = 1, ao_cart_num

    ii = ao_cart_nucl(i)
    powA(1) = ao_cart_power(i,1)
    powA(2) = ao_cart_power(i,2)
    powA(3) = ao_cart_power(i,3)
 
    if(primitives_normalized) then

      ! Normalization of the primitives
      do j = 1, ao_cart_prim_num(i)

        expo = ao_cart_expo(i,j) + (0.d0, 1.d0) * ao_cart_expo_im(i,j)
        expo_inv = (1.d0, 0.d0) / expo
        do m = 1, 3
          C_Ap(m) = nucl_coord(ii,m)
          C_Ae(m) = nucl_coord(ii,m) - (0.d0, 0.5d0) * expo_inv * ao_cart_expo_pw(m,i,j)
        enddo
        phiA = ao_cart_expo_phase(1,i,j) + ao_cart_expo_phase(2,i,j) + ao_cart_expo_phase(3,i,j)
        KA2 = ao_cart_expo_pw(1,i,j) * ao_cart_expo_pw(1,i,j) &
            + ao_cart_expo_pw(2,i,j) * ao_cart_expo_pw(2,i,j) &
            + ao_cart_expo_pw(3,i,j) * ao_cart_expo_pw(3,i,j)

        C1 = zexp(-(0.d0, 2.d0) * phiA - 0.5d0 * expo_inv * KA2)
        C2 = zexp(-(0.5d0, 0.d0) * real(expo_inv) * KA2)

        call overlap_cgaussian_xyz(C_Ae, C_Ae, expo, expo, powA, powA, &
                                   C_Ap, C_Ap, overlap_x, overlap_y, overlap_z, integ1, nz)

        call overlap_cgaussian_xyz(conjg(C_Ae), C_Ae, conjg(expo), expo, powA, powA, &
                                   conjg(C_Ap), C_Ap, overlap_x, overlap_y, overlap_z, integ2, nz)

        norm = 2.d0 * real(C1 * integ1 + C2 * integ2)

        !ao_cart_coef_norm_cgtos(i,j) = 1.d0 / dsqrt(norm)
        ao_cart_coef_norm_cgtos(i,j) = ao_cart_coef(i,j) / dsqrt(norm)
      enddo

    else

      do j = 1, ao_cart_prim_num(i)
        ao_cart_coef_norm_cgtos(i,j) = ao_cart_coef(i,j)
      enddo

    endif ! primitives_normalized

  enddo

END_PROVIDER


