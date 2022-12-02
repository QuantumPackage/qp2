use bitmasks

BEGIN_PROVIDER [ double precision, psi_l_coef_bi_ortho, (psi_det_size,N_states) ]
  implicit none
  BEGIN_DOC
  ! The wave function coefficients. Initialized with Hartree-Fock if the |EZFIO| file
  ! is empty.
  END_DOC

  integer                        :: i,k, N_int2
  logical                        :: exists
  character*(64)                 :: label

  PROVIDE read_wf N_det mo_label ezfio_filename nproc
  psi_l_coef_bi_ortho = 0.d0
  do i=1,min(N_states,N_det)
    psi_l_coef_bi_ortho(i,i) = 1.d0
  enddo

  if (mpi_master) then
    if (read_wf) then
      call ezfio_has_tc_bi_ortho_psi_l_coef_bi_ortho(exists)
!      if (exists) then
!        call ezfio_has_tc_bi_ortho_mo_label(exists)
!        if (exists) then
!          call ezfio_get_tc_bi_ortho_mo_label(label)
!          exists = (label == mo_label)
!        endif
!      endif

      if (exists) then

        double precision, allocatable  :: psi_l_coef_bi_ortho_read(:,:)
        allocate (psi_l_coef_bi_ortho_read(N_det,N_states))
        print *,  'Read psi_l_coef_bi_ortho', N_det, N_states
        call ezfio_get_tc_bi_ortho_psi_l_coef_bi_ortho(psi_l_coef_bi_ortho_read)
        do k=1,N_states
          do i=1,N_det
            psi_l_coef_bi_ortho(i,k) = psi_l_coef_bi_ortho_read(i,k)
          enddo
        enddo
        deallocate(psi_l_coef_bi_ortho_read)

      else

        print*, 'psi_l_coef_bi_ortho are psi_coef'
        do k=1,N_states
          do i=1,N_det
            psi_l_coef_bi_ortho(i,k) = psi_coef(i,k)
          enddo
        enddo

      endif
    endif
  endif
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer                        :: ierr
    call     MPI_BCAST( psi_l_coef_bi_ortho, size(psi_l_coef_bi_ortho), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read psi_l_coef_bi_ortho with MPI'
    endif
  IRP_ENDIF
END_PROVIDER


BEGIN_PROVIDER [ double precision, psi_r_coef_bi_ortho, (psi_det_size,N_states) ]
  implicit none
  BEGIN_DOC
  ! The wave function coefficients. Initialized with Hartree-Fock if the |EZFIO| file
  ! is empty.
  END_DOC

  integer                        :: i,k, N_int2
  logical                        :: exists
  character*(64)                 :: label

  PROVIDE read_wf N_det mo_label ezfio_filename nproc
  psi_r_coef_bi_ortho = 0.d0
  do i=1,min(N_states,N_det)
    psi_r_coef_bi_ortho(i,i) = 1.d0
  enddo

  if (mpi_master) then
    if (read_wf) then
      call ezfio_has_tc_bi_ortho_psi_r_coef_bi_ortho(exists)
!      if (exists) then
!        call ezfio_has_tc_bi_ortho_mo_label(exists)
!        if (exists) then
!          call ezfio_get_tc_bi_ortho_mo_label(label)
!          exists = (label == mo_label)
!        endif
!      endif

      if (exists) then

        double precision, allocatable  :: psi_r_coef_bi_ortho_read(:,:)
        allocate (psi_r_coef_bi_ortho_read(N_det,N_states))
        print *,  'Read psi_r_coef_bi_ortho', N_det, N_states
        call ezfio_get_tc_bi_ortho_psi_r_coef_bi_ortho(psi_r_coef_bi_ortho_read)
        do k=1,N_states
          do i=1,N_det
            psi_r_coef_bi_ortho(i,k) = psi_r_coef_bi_ortho_read(i,k)
          enddo
        enddo
        deallocate(psi_r_coef_bi_ortho_read)

      else

        print*, 'psi_r_coef_bi_ortho are psi_coef'
        do k=1,N_states
          do i=1,N_det
            psi_r_coef_bi_ortho(i,k) = psi_coef(i,k)
          enddo
        enddo

      endif
    endif
  endif
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer                        :: ierr
    call     MPI_BCAST( psi_r_coef_bi_ortho, size(psi_r_coef_bi_ortho), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read psi_r_coef_bi_ortho with MPI'
    endif
  IRP_ENDIF
END_PROVIDER


subroutine save_tc_wavefunction_general(ndet,nstates,psidet,dim_psicoef,psilcoef,psircoef)
  implicit none
  BEGIN_DOC
  !  Save the wave function into the |EZFIO| file
  END_DOC
  use bitmasks
  include 'constants.include.F'
  integer, intent(in)            :: ndet,nstates,dim_psicoef
  integer(bit_kind), intent(in)  :: psidet(N_int,2,ndet)
  double precision, intent(in)   :: psilcoef(dim_psicoef,nstates)
  double precision, intent(in)   :: psircoef(dim_psicoef,nstates)
  integer*8, allocatable         :: psi_det_save(:,:,:)
  double precision, allocatable  :: psil_coef_save(:,:)
  double precision, allocatable  :: psir_coef_save(:,:)

  double precision               :: accu_norm
  integer                        :: i,j,k, ndet_qp_edit

  if (mpi_master) then
    ndet_qp_edit = min(ndet,N_det_qp_edit)

    call ezfio_set_determinants_N_int(N_int)
    call ezfio_set_determinants_bit_kind(bit_kind)
    call ezfio_set_determinants_N_det(ndet)
    call ezfio_set_determinants_N_det_qp_edit(ndet_qp_edit)
    call ezfio_set_determinants_n_states(nstates)
    call ezfio_set_determinants_mo_label(mo_label)

    allocate (psi_det_save(N_int,2,ndet))
    do i=1,ndet
      do j=1,2
        do k=1,N_int
          psi_det_save(k,j,i) = transfer(psidet(k,j,i),1_8)
        enddo
      enddo
    enddo
    call ezfio_set_determinants_psi_det(psi_det_save)
    call ezfio_set_determinants_psi_det_qp_edit(psi_det_save)
    deallocate (psi_det_save)

    allocate (psil_coef_save(ndet,nstates),psir_coef_save(ndet,nstates))
    do k=1,nstates
      do i=1,ndet
        psil_coef_save(i,k) = psilcoef(i,k)
        psir_coef_save(i,k) = psircoef(i,k)
      enddo
    enddo

    call ezfio_set_tc_bi_ortho_psi_l_coef_bi_ortho(psil_coef_save)
    call ezfio_set_tc_bi_ortho_psi_r_coef_bi_ortho(psir_coef_save)
    deallocate (psil_coef_save,psir_coef_save)

!    allocate (psi_coef_save(ndet_qp_edit,nstates))
!    do k=1,nstates
!      do i=1,ndet_qp_edit
!        psi_coef_save(i,k) = psicoef(i,k)
!      enddo
!    enddo
!
!    call ezfio_set_determinants_psi_coef_qp_edit(psi_coef_save)
!    deallocate (psi_coef_save)

    call write_int(6,ndet,'Saved determinantsi and psi_r/psi_l coef')
  endif
end

subroutine save_tc_bi_ortho_wavefunction
 implicit none
 call save_tc_wavefunction_general(N_det,N_states,psi_det,size(psi_l_coef_bi_ortho, 1),psi_l_coef_bi_ortho,psi_r_coef_bi_ortho)
 call routine_save_right_bi_ortho
end

subroutine routine_save_right_bi_ortho
 implicit none
 double precision, allocatable :: coef_tmp(:,:)
 integer :: i
 N_states = 1
 allocate(coef_tmp(N_det, N_states))
 do i = 1, N_det
  coef_tmp(i,1) = psi_r_coef_bi_ortho(i,1)
 enddo
 call save_wavefunction_general_unormalized(N_det,N_states,psi_det,size(coef_tmp,1),coef_tmp(1,1))
end                     

subroutine routine_save_left_right_bi_ortho
 implicit none
 double precision, allocatable :: coef_tmp(:,:)
 integer :: i,n_states_tmp
 n_states_tmp = 2
 allocate(coef_tmp(N_det, n_states_tmp))
 do i = 1, N_det
  coef_tmp(i,1) = psi_r_coef_bi_ortho(i,1)
  coef_tmp(i,2) = psi_l_coef_bi_ortho(i,1)
 enddo
 call save_wavefunction_general_unormalized(N_det,n_states_tmp,psi_det,size(coef_tmp,1),coef_tmp(1,1))
end

