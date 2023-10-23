subroutine export_trexio(update)
  use trexio
  implicit none
  BEGIN_DOC
  !     Exports the wave function in TREXIO format
  END_DOC

  logical, intent(in)            :: update
  integer(trexio_t)              :: f(N_states) ! TREXIO file handle
  integer(trexio_exit_code)      :: rc
  integer                        :: k
  double precision, allocatable  :: factor(:)
  character*(256)  :: filenames(N_states)
  character :: rw

  filenames(1) = trexio_filename
  do k=2,N_states
    write(filenames(k),'(A,I3.3)') trim(trexio_filename)//'.', k-1
  enddo

  do k=1,N_states
    print *, 'TREXIO file : ', trim(filenames(k))
    if (update) then
      call system('test -f '//trim(filenames(k))//' && cp -r '//trim(filenames(k))//' '//trim(filenames(k))//'.bak')
    else
      call system('test -f '//trim(filenames(k))//' && mv '//trim(filenames(k))//' '//trim(filenames(k))//'.bak')
    endif
  enddo
  print *, ''

  if (update) then
     rw = 'u'
  else
     rw = 'w'
  endif


  do k=1,N_states
    if (backend == 0) then
      f(k) = trexio_open(filenames(k), rw, TREXIO_HDF5, rc)
    else if (backend == 1) then
      f(k) = trexio_open(filenames(k), rw, TREXIO_TEXT, rc)
    endif
    if (f(k) == 0_8) then
      print *, 'Unable to open TREXIO file for writing'
      print *, 'rc = ', rc
      stop -1
    endif
  enddo
  call ezfio_set_trexio_trexio_file(trexio_filename)

! ------------------------------------------------------------------------------

! Electrons
! ---------

  print *, 'Electrons'

  rc = trexio_write_electron_up_num(f(1), elec_alpha_num)
  call trexio_assert(rc, TREXIO_SUCCESS)

  rc = trexio_write_electron_dn_num(f(1), elec_beta_num)
  call trexio_assert(rc, TREXIO_SUCCESS)


! Nuclei
! ------

  print *, 'Nuclei'

  rc = trexio_write_nucleus_num(f(1), nucl_num)
  call trexio_assert(rc, TREXIO_SUCCESS)

  rc = trexio_write_nucleus_charge(f(1), nucl_charge)
  call trexio_assert(rc, TREXIO_SUCCESS)

  rc = trexio_write_nucleus_coord(f(1), nucl_coord_transp)
  call trexio_assert(rc, TREXIO_SUCCESS)

  rc = trexio_write_nucleus_label(f(1), nucl_label, 32)
  call trexio_assert(rc, TREXIO_SUCCESS)

  rc = trexio_write_nucleus_repulsion(f(1), nuclear_repulsion)
  call trexio_assert(rc, TREXIO_SUCCESS)


! Pseudo-potentials
! -----------------

  if (do_pseudo) then

    print *, 'ECP'
    integer                       :: num

    num = 0
    do k=1,pseudo_klocmax
      do i=1,nucl_num
        if (pseudo_dz_k(i,k) /= 0.d0) then
          num = num+1
        end if
      end do
    end do

    do l=0,pseudo_lmax
      do k=1,pseudo_kmax
        do i=1,nucl_num
          if (pseudo_dz_kl(i,k,l) /= 0.d0) then
            num = num+1
          end if
        end do
      end do
    end do

    integer, allocatable          :: ang_mom(:), nucleus_index(:), power(:), lmax(:)
    double precision, allocatable :: exponent(:), coefficient(:)

    allocate(ang_mom(num), nucleus_index(num), exponent(num), coefficient(num), power(num), &
            lmax(nucl_num) )

    do i=1,nucl_num
      lmax(i) = -1
      do l=0,pseudo_lmax
        do k=1,pseudo_kmax
          if (pseudo_dz_kl_transp(k,l,i) /= 0.d0) then
            lmax(i) = max(lmax(i), l)
          end if
        end do
      end do
    end do

    j = 0
    do i=1,nucl_num
      do k=1,pseudo_klocmax
        if (pseudo_dz_k_transp(k,i) /= 0.d0) then
          j = j+1
          ang_mom(j) = lmax(i)+1
          nucleus_index(j) = i
          exponent(j) = pseudo_dz_k_transp(k,i)
          coefficient(j) = pseudo_v_k_transp(k,i)
          power(j) = pseudo_n_k_transp(k,i)
        end if
      end do

      do l=0,lmax(i)
        do k=1,pseudo_kmax
          if (pseudo_dz_kl_transp(k,l,i) /= 0.d0) then
            j = j+1
            ang_mom(j) = l
            nucleus_index(j) = i
            exponent(j) = pseudo_dz_kl_transp(k,l,i)
            coefficient(j) = pseudo_v_kl_transp(k,l,i)
            power(j) = pseudo_n_kl_transp(k,l,i)
          end if
        end do
      end do
    end do


    lmax(:) = lmax(:)+1
    rc = trexio_write_ecp_max_ang_mom_plus_1(f(1), lmax)
    call trexio_assert(rc, TREXIO_SUCCESS)

    rc = trexio_write_ecp_z_core(f(1), int(nucl_charge_remove))
    call trexio_assert(rc, TREXIO_SUCCESS)

    rc = trexio_write_ecp_num(f(1), num)
    call trexio_assert(rc, TREXIO_SUCCESS)

    rc = trexio_write_ecp_ang_mom(f(1), ang_mom)
    call trexio_assert(rc, TREXIO_SUCCESS)

    rc = trexio_write_ecp_nucleus_index(f(1), nucleus_index)
    call trexio_assert(rc, TREXIO_SUCCESS)

    rc = trexio_write_ecp_exponent(f(1), exponent)
    call trexio_assert(rc, TREXIO_SUCCESS)

    rc = trexio_write_ecp_coefficient(f(1), coefficient)
    call trexio_assert(rc, TREXIO_SUCCESS)

    rc = trexio_write_ecp_power(f(1), power)
    call trexio_assert(rc, TREXIO_SUCCESS)

  endif


  if (export_basis) then

! Basis
! -----

    print *, 'Basis'

    rc = trexio_write_basis_type(f(1), 'Gaussian', len('Gaussian'))
    call trexio_assert(rc, TREXIO_SUCCESS)

    rc = trexio_write_basis_prim_num(f(1), prim_num)
    call trexio_assert(rc, TREXIO_SUCCESS)

     rc = trexio_write_basis_shell_num(f(1), shell_num)
     call trexio_assert(rc, TREXIO_SUCCESS)

     rc = trexio_write_basis_nucleus_index(f(1), basis_nucleus_index)
     call trexio_assert(rc, TREXIO_SUCCESS)

     rc = trexio_write_basis_shell_ang_mom(f(1), shell_ang_mom)
     call trexio_assert(rc, TREXIO_SUCCESS)

     allocate(factor(shell_num))
!     if (ao_normalized) then
!       factor(1:shell_num) = shell_normalization_factor(1:shell_num)
!     else
       factor(1:shell_num) = 1.d0
!     endif
     rc = trexio_write_basis_shell_factor(f(1), factor)
     call trexio_assert(rc, TREXIO_SUCCESS)

     deallocate(factor)

    rc = trexio_write_basis_shell_index(f(1), shell_index)
    call trexio_assert(rc, TREXIO_SUCCESS)

    rc = trexio_write_basis_exponent(f(1), prim_expo)
    call trexio_assert(rc, TREXIO_SUCCESS)

    rc = trexio_write_basis_coefficient(f(1), prim_coef)
    call trexio_assert(rc, TREXIO_SUCCESS)

    allocate(factor(prim_num))
    if (primitives_normalized) then
      factor(1:prim_num) = prim_normalization_factor(1:prim_num)
    else
      factor(1:prim_num) = 1.d0
    endif
    rc = trexio_write_basis_prim_factor(f(1), factor)
    call trexio_assert(rc, TREXIO_SUCCESS)
    deallocate(factor)


! Atomic orbitals
! ---------------

    print *, 'AOs'

    rc = trexio_write_ao_num(f(1), ao_num)
    call trexio_assert(rc, TREXIO_SUCCESS)

    rc = trexio_write_ao_cartesian(f(1), 1)
    call trexio_assert(rc, TREXIO_SUCCESS)

    rc = trexio_write_ao_shell(f(1), ao_shell)
    call trexio_assert(rc, TREXIO_SUCCESS)

    integer :: i, pow0(3), powA(3), j, l, nz
    double precision :: normA, norm0, C_A(3), overlap_x, overlap_z, overlap_y, c
    nz=100

    C_A(1) = 0.d0
    C_A(2) = 0.d0
    C_A(3) = 0.d0

    allocate(factor(ao_num))
    if (ao_normalized) then
      do i=1,ao_num
        l = ao_first_of_shell(ao_shell(i))
        factor(i) = (ao_coef_normalized(i,1)+tiny(1.d0))/(ao_coef_normalized(l,1)+tiny(1.d0))
      enddo
    else
      factor(:) = 1.d0
    endif
    rc = trexio_write_ao_normalization(f(1), factor)
    call trexio_assert(rc, TREXIO_SUCCESS)
    deallocate(factor)

  endif

! One-e AO integrals
! ------------------

  if (export_ao_one_e_ints) then
    print *, 'AO one-e integrals'

    rc = trexio_write_ao_1e_int_overlap(f(1),ao_overlap)
    call trexio_assert(rc, TREXIO_SUCCESS)

    rc = trexio_write_ao_1e_int_kinetic(f(1),ao_kinetic_integrals)
    call trexio_assert(rc, TREXIO_SUCCESS)

    rc = trexio_write_ao_1e_int_potential_n_e(f(1),ao_integrals_n_e)
    call trexio_assert(rc, TREXIO_SUCCESS)

    if (do_pseudo) then
      rc = trexio_write_ao_1e_int_ecp(f(1), ao_pseudo_integrals_local + ao_pseudo_integrals_non_local)
      call trexio_assert(rc, TREXIO_SUCCESS)
    endif

    rc = trexio_write_ao_1e_int_core_hamiltonian(f(1),ao_one_e_integrals)
    call trexio_assert(rc, TREXIO_SUCCESS)
  end if

! Two-e AO integrals
! ------------------

  if (export_ao_two_e_ints) then
    print *, 'AO two-e integrals'
    PROVIDE ao_two_e_integrals_in_map

    integer(8), parameter :: BUFSIZE=100000_8
    double precision :: eri_buffer(BUFSIZE), integral
    integer(4) :: eri_index(4,BUFSIZE)
    integer(8) :: icount, offset

    double precision, external :: get_ao_two_e_integral


    icount = 0_8
    offset = 0_8
    do l=1,ao_num
      do k=1,ao_num
        do j=l,ao_num
          do i=k,ao_num
            if (i==j .and. k<l) cycle
            if (i<j) cycle
            integral = get_ao_two_e_integral(i,j,k,l,ao_integrals_map)
            if (integral == 0.d0) cycle
            icount += 1_8
            eri_buffer(icount) = integral
            eri_index(1,icount) = i
            eri_index(2,icount) = j
            eri_index(3,icount) = k
            eri_index(4,icount) = l
            if (icount == BUFSIZE) then
              rc = trexio_write_ao_2e_int_eri(f(1), offset, icount, eri_index, eri_buffer)
              call trexio_assert(rc, TREXIO_SUCCESS)
              offset += icount
              icount = 0_8
            end if
          end do
        end do
      end do
    end do

    if (icount >= 0_8) then
      rc = trexio_write_ao_2e_int_eri(f(1), offset, icount, eri_index, eri_buffer)
      call trexio_assert(rc, TREXIO_SUCCESS)
    end if
  end if

! Two-e AO integrals - Cholesky
! -----------------------------

  integer(4) :: chol_index(3,BUFSIZE)
  double precision :: chol_buffer(BUFSIZE)

  if (export_ao_two_e_ints_cholesky) then
    print *, 'AO two-e integrals Cholesky'

    rc = trexio_write_ao_2e_int_eri_cholesky_num(f(1), cholesky_ao_num)
    call trexio_assert(rc, TREXIO_SUCCESS)

    icount = 0_8
    offset = 0_8
    do k=1,cholesky_ao_num
     do j=1,ao_num
      do i=1,ao_num
         integral = cholesky_ao(i,j,k)
         if (integral == 0.d0) cycle
         icount += 1_8
         chol_buffer(icount) = integral
         chol_index(1,icount) = i
         chol_index(2,icount) = j
         chol_index(3,icount) = k
         if (icount == BUFSIZE) then
           rc = trexio_write_ao_2e_int_eri_cholesky(f(1), offset, icount, chol_index, chol_buffer)
           call trexio_assert(rc, TREXIO_SUCCESS)
           offset += icount
           icount = 0_8
         end if
      end do
     end do
    end do

    if (icount > 0_8) then
      rc = trexio_write_ao_2e_int_eri_cholesky(f(1), offset, icount, chol_index, chol_buffer)
      call trexio_assert(rc, TREXIO_SUCCESS)
    end if
  end if



! Molecular orbitals
! ------------------

  if (export_mos) then
    print *, 'MOs'

    rc = trexio_write_mo_type(f(1), mo_label, len(trim(mo_label)))
    call trexio_assert(rc, TREXIO_SUCCESS)

    do k=1,N_states
      rc = trexio_write_mo_num(f(k), mo_num)
      call trexio_assert(rc, TREXIO_SUCCESS)
    enddo

    rc = trexio_write_mo_coefficient(f(1), mo_coef)
    call trexio_assert(rc, TREXIO_SUCCESS)

    if ( (trim(mo_label) == 'Canonical').and. &
         (export_mo_two_e_ints_cholesky.or.export_mo_two_e_ints) ) then
      rc = trexio_write_mo_energy(f(1), fock_matrix_diag_mo)
      call trexio_assert(rc, TREXIO_SUCCESS)
    endif

    rc = trexio_write_mo_class(f(1), mo_class, len(mo_class(1)))
    call trexio_assert(rc, TREXIO_SUCCESS)
  endif

! One-e MO integrals
! ------------------

  if (export_mo_one_e_ints) then
    print *, 'MO one-e integrals'

    rc = trexio_write_mo_1e_int_kinetic(f(1),mo_kinetic_integrals)
    call trexio_assert(rc, TREXIO_SUCCESS)

    rc = trexio_write_mo_1e_int_potential_n_e(f(1),mo_integrals_n_e)
    call trexio_assert(rc, TREXIO_SUCCESS)

    if (do_pseudo) then
      rc = trexio_write_mo_1e_int_ecp(f(1),mo_pseudo_integrals_local)
      call trexio_assert(rc, TREXIO_SUCCESS)
    endif

    rc = trexio_write_mo_1e_int_core_hamiltonian(f(1),mo_one_e_integrals)
    call trexio_assert(rc, TREXIO_SUCCESS)
  end if

! Two-e MO integrals
! ------------------

  if (export_mo_two_e_ints) then
    print *, 'MO two-e integrals'
    PROVIDE mo_two_e_integrals_in_map

    double precision, external :: mo_two_e_integral


    icount = 0_8
    offset = 0_8
    do l=1,mo_num
      do k=1,mo_num
        do j=l,mo_num
          do i=k,mo_num
            if (i==j .and. k<l) cycle
            if (i<j) cycle
            integral = mo_two_e_integral(i,j,k,l)
            if (integral == 0.d0) cycle
            icount += 1_8
            eri_buffer(icount) = integral
            eri_index(1,icount) = i
            eri_index(2,icount) = j
            eri_index(3,icount) = k
            eri_index(4,icount) = l
            if (icount == BUFSIZE) then
              rc = trexio_write_mo_2e_int_eri(f(1), offset, icount, eri_index, eri_buffer)
              call trexio_assert(rc, TREXIO_SUCCESS)
              offset += icount
              icount = 0_8
            end if
          end do
        end do
      end do
    end do

    if (icount > 0_8) then
      rc = trexio_write_mo_2e_int_eri(f(1), offset, icount, eri_index, eri_buffer)
      call trexio_assert(rc, TREXIO_SUCCESS)
    end if
  end if

! Two-e MO integrals - Cholesky
! -----------------------------

  if (export_mo_two_e_ints_cholesky) then
    print *, 'MO two-e integrals Cholesky'

    rc = trexio_write_mo_2e_int_eri_cholesky_num(f(1), cholesky_ao_num)
    call trexio_assert(rc, TREXIO_SUCCESS)

    icount = 0_8
    offset = 0_8
    do k=1,cholesky_ao_num
     do j=1,mo_num
      do i=1,mo_num
         integral = cholesky_mo(i,j,k)
         if (integral == 0.d0) cycle
         icount += 1_8
         chol_buffer(icount) = integral
         chol_index(1,icount) = i
         chol_index(2,icount) = j
         chol_index(3,icount) = k
         if (icount == BUFSIZE) then
           rc = trexio_write_mo_2e_int_eri_cholesky(f(1), offset, icount, chol_index, chol_buffer)
           call trexio_assert(rc, TREXIO_SUCCESS)
           offset += icount
           icount = 0_8
         end if
      end do
     end do
    end do

    if (icount > 0_8) then
      rc = trexio_write_mo_2e_int_eri_cholesky(f(1), offset, icount, chol_index, chol_buffer)
      call trexio_assert(rc, TREXIO_SUCCESS)
    end if
  end if


! One-e RDM
! ---------

  do k=1,N_states
    rc = trexio_write_rdm_1e(f(k),one_e_dm_mo_alpha(:,:,k) + one_e_dm_mo_beta(:,:,k))
    call trexio_assert(rc, TREXIO_SUCCESS)

    rc = trexio_write_rdm_1e_up(f(k),one_e_dm_mo_alpha(:,:,k))
    call trexio_assert(rc, TREXIO_SUCCESS)

    rc = trexio_write_rdm_1e_dn(f(k),one_e_dm_mo_beta(:,:,k))
    call trexio_assert(rc, TREXIO_SUCCESS)
  enddo


! Two-e RDM
! ---------

  if (export_rdm) then
    PROVIDE two_e_dm_mo
    print *, 'Two-e RDM'

    icount = 0_8
    offset = 0_8
    do l=1,mo_num
      do k=1,mo_num
        do j=1,mo_num
          do i=1,mo_num
            integral = two_e_dm_mo(i,j,k,l)
            if (integral == 0.d0) cycle
            icount += 1_8
            eri_buffer(icount) = integral
            eri_index(1,icount) = i
            eri_index(2,icount) = j
            eri_index(3,icount) = k
            eri_index(4,icount) = l
            if (icount == BUFSIZE) then
              rc = trexio_write_rdm_2e(f(1), offset, icount, eri_index, eri_buffer)
              call trexio_assert(rc, TREXIO_SUCCESS)
              offset += icount
              icount = 0_8
            end if
          end do
        end do
      end do
    end do

    if (icount >= 0_8) then
      rc = trexio_write_rdm_2e(f(1), offset, icount, eri_index, eri_buffer)
      call trexio_assert(rc, TREXIO_SUCCESS)
    end if
  end if


! ------------------------------------------------------------------------------

  ! Determinants
  ! ------------

    integer*8, allocatable :: det_buffer(:,:,:)
    double precision, allocatable :: coef_buffer(:,:)
    integer :: nint

    rc = trexio_get_int64_num(f(1), nint)
    call trexio_assert(rc, TREXIO_SUCCESS)
!    nint = N_int
    if (nint /= N_int) then
       stop 'Problem with N_int'
    endif
    allocate ( det_buffer(nint, 2, BUFSIZE), coef_buffer(BUFSIZE, n_states) )

    do k=1, N_states
      icount = 0_8
      offset = 0_8
      rc = trexio_write_state_num(f(k), n_states)
      call trexio_assert(rc, TREXIO_SUCCESS)

! Will need to be updated with TREXIO 2.4
!      rc = trexio_write_state_id(f(k), k-1)
      rc = trexio_write_state_id(f(k), k)
      call trexio_assert(rc, TREXIO_SUCCESS)

      rc = trexio_write_state_file_name(f(k), filenames, len(filenames(1)))
      call trexio_assert(rc, TREXIO_SUCCESS)
    enddo

    do k=1,n_det
       icount += 1_8
       det_buffer(1:nint, 1:2, icount) = psi_det(1:N_int, 1:2, k)
       coef_buffer(icount,1:N_states) = psi_coef(k,1:N_states)
       if (icount == BUFSIZE) then
         do i=1,N_states
           rc = trexio_write_determinant_list(f(i), offset, icount, det_buffer)
           call trexio_assert(rc, TREXIO_SUCCESS)
           rc = trexio_write_determinant_coefficient(f(i), offset, icount, coef_buffer(1,i))
           call trexio_assert(rc, TREXIO_SUCCESS)
         end do
         offset += icount
         icount = 0_8
       end if
    end do

  if (icount >= 0_8) then
    do i=1,N_states
      rc = trexio_write_determinant_list(f(i), offset, icount, det_buffer)
      call trexio_assert(rc, TREXIO_SUCCESS)
      rc = trexio_write_determinant_coefficient(f(i), offset, icount, coef_buffer(1,i))
      call trexio_assert(rc, TREXIO_SUCCESS)
    end do
  end if

  deallocate ( det_buffer, coef_buffer )

  do k=1,N_states
    rc = trexio_close(f(k))
    call trexio_assert(rc, TREXIO_SUCCESS)
  enddo

end


! -*- mode: f90 -*-
