use bitmasks




BEGIN_PROVIDER [ complex*16, psi_coef_complex, (psi_det_size,N_states) ]
  implicit none
  BEGIN_DOC
  ! The wave function coefficients. Initialized with Hartree-Fock if the |EZFIO| file
  ! is empty.
  END_DOC

  integer                        :: i,k, N_int2
  logical                        :: exists
  character*(64)                 :: label

  PROVIDE read_wf N_det mo_label ezfio_filename
  psi_coef = (0.d0,0.d0)
  do i=1,min(N_states,psi_det_size)
    psi_coef(i,i) = (1.d0,0.d0)
  enddo

  if (mpi_master) then
    if (read_wf) then
      call ezfio_has_determinants_psi_coef_complex(exists)
      if (exists) then
        call ezfio_has_determinants_mo_label(exists)
        if (exists) then
          call ezfio_get_determinants_mo_label(label)
          exists = (label == mo_label)
        endif
      endif

      if (exists) then

        complex*16, allocatable  :: psi_coef_read(:,:)
        allocate (psi_coef_read(N_det,N_states))
        print *,  'Read psi_coef_complex', N_det, N_states
        call ezfio_get_determinants_psi_coef_complex(psi_coef_read)
        do k=1,N_states
          do i=1,N_det
            psi_coef_complex(i,k) = psi_coef_read(i,k)
          enddo
        enddo
        deallocate(psi_coef_read)

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
    call     MPI_BCAST( psi_coef_complex, size(psi_coef_complex), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read psi_coef_complex with MPI'
    endif
  IRP_ENDIF



END_PROVIDER

!==============================================================================!
!                                                                              !
!                               Sorting providers                              !
!                                                                              !
!==============================================================================!

!TODO: implement for complex (new psi_det_sorted? reuse? combine complex provider with real?)
 BEGIN_PROVIDER [ integer(bit_kind), psi_det_sorted_complex, (N_int,2,psi_det_size) ]
&BEGIN_PROVIDER [ complex*16, psi_coef_sorted_complex, (psi_det_size,N_states) ]
&BEGIN_PROVIDER [ double precision, psi_average_norm_contrib_sorted_complex, (psi_det_size) ]
&BEGIN_PROVIDER [ integer, psi_det_sorted_order_complex, (psi_det_size) ]
   implicit none
   BEGIN_DOC
   ! Wave function sorted by determinants contribution to the norm (state-averaged)
   !
   ! psi_det_sorted_order(i) -> k : index in psi_det
   END_DOC
   integer                        :: i,j,k
   integer, allocatable           :: iorder(:)
   allocate ( iorder(N_det) )
   do i=1,N_det
     psi_average_norm_contrib_sorted_complex(i) = -psi_average_norm_contrib(i)
     iorder(i) = i
   enddo
   call dsort(psi_average_norm_contrib_sorted_complex,iorder,N_det)
   do i=1,N_det
     do j=1,N_int
       psi_det_sorted_complex(j,1,i) = psi_det(j,1,iorder(i))
       psi_det_sorted_complex(j,2,i) = psi_det(j,2,iorder(i))
     enddo
     do k=1,N_states
       psi_coef_sorted_complex(i,k) = psi_coef_complex(iorder(i),k)
     enddo
     psi_average_norm_contrib_sorted_complex(i) = -psi_average_norm_contrib_sorted_complex(i)
   enddo
   do i=1,N_det
     psi_det_sorted_order_complex(iorder(i)) = i
   enddo

   psi_det_sorted_complex(:,:,N_det+1:psi_det_size) = 0_bit_kind
   psi_coef_sorted_complex(N_det+1:psi_det_size,:) = (0.d0,0.d0)
   psi_average_norm_contrib_sorted_complex(N_det+1:psi_det_size) = 0.d0
   psi_det_sorted_order_complex(N_det+1:psi_det_size) = 0

   deallocate(iorder)

END_PROVIDER

 BEGIN_PROVIDER [ integer(bit_kind), psi_det_sorted_bit_complex, (N_int,2,psi_det_size) ]
&BEGIN_PROVIDER [ complex*16, psi_coef_sorted_bit_complex, (psi_det_size,N_states) ]
   implicit none
   BEGIN_DOC
   ! Determinants on which we apply $\langle i|H|psi \rangle$ for perturbation.
   ! They are sorted by determinants interpreted as integers. Useful
   ! to accelerate the search of a random determinant in the wave
   ! function.
   END_DOC

   call sort_dets_by_det_search_key_complex(N_det, psi_det, psi_coef_complex, &
                                 size(psi_coef_complex,1), psi_det_sorted_bit_complex, &
                                 psi_coef_sorted_bit_complex, N_states)

END_PROVIDER

subroutine sort_dets_by_det_search_key_complex(Ndet, det_in, coef_in, sze, det_out, coef_out, N_st)
   use bitmasks
   implicit none
   integer, intent(in)            :: Ndet, N_st, sze
   integer(bit_kind), intent(in)  :: det_in  (N_int,2,sze)
   complex*16 , intent(in)  :: coef_in(sze,N_st)
   integer(bit_kind), intent(out) :: det_out (N_int,2,sze)
   complex*16 , intent(out) :: coef_out(sze,N_st)
   BEGIN_DOC
   ! Determinants are sorted according to their :c:func:`det_search_key`.
   ! Useful to accelerate the search of a random determinant in the wave
   ! function.
   !
   ! /!\ The first dimension of coef_out and coef_in need to be psi_det_size
   !
   END_DOC
   integer                        :: i,j,k
   integer, allocatable           :: iorder(:)
   integer*8, allocatable         :: bit_tmp(:)
   integer*8, external            :: det_search_key

   allocate ( iorder(Ndet), bit_tmp(Ndet) )

   do i=1,Ndet
     iorder(i) = i
     !$DIR FORCEINLINE
     bit_tmp(i) = det_search_key(det_in(1,1,i),N_int)
   enddo
   call i8sort(bit_tmp,iorder,Ndet)
   !DIR$ IVDEP
   do i=1,Ndet
     do j=1,N_int
       det_out(j,1,i) = det_in(j,1,iorder(i))
       det_out(j,2,i) = det_in(j,2,iorder(i))
     enddo
     do k=1,N_st
       coef_out(i,k) = coef_in(iorder(i),k)
     enddo
   enddo

   deallocate(iorder, bit_tmp)

end


! TODO:complex? only keep abs max/min? real max/min?
! BEGIN_PROVIDER [ double precision, psi_coef_max, (N_states) ]
!&BEGIN_PROVIDER [ double precision, psi_coef_min, (N_states) ]
!&BEGIN_PROVIDER [ double precision, abs_psi_coef_max, (N_states) ]
!&BEGIN_PROVIDER [ double precision, abs_psi_coef_min, (N_states) ]
!   implicit none
!   BEGIN_DOC
!   ! Max and min values of the coefficients
!   END_DOC
!   integer                        :: i
!   do i=1,N_states
!     psi_coef_min(i) = minval(psi_coef(:,i))
!     psi_coef_max(i) = maxval(psi_coef(:,i))
!     abs_psi_coef_min(i) = minval( dabs(psi_coef(:,i)) )
!     abs_psi_coef_max(i) = maxval( dabs(psi_coef(:,i)) )
!     call write_double(6,psi_coef_max(i), 'Max coef')
!     call write_double(6,psi_coef_min(i), 'Min coef')
!     call write_double(6,abs_psi_coef_max(i), 'Max abs coef')
!     call write_double(6,abs_psi_coef_min(i), 'Min abs coef')
!   enddo
!
!END_PROVIDER


!==============================================================================!
!                                                                              !
!                             Read/write routines                              !
!                                                                              !
!==============================================================================!



subroutine save_wavefunction_general_complex(ndet,nstates,psidet,dim_psicoef,psicoef)
  implicit none
  BEGIN_DOC
  !  Save the wave function into the |EZFIO| file
  END_DOC
  use bitmasks
  include 'constants.include.F'
  integer, intent(in)            :: ndet,nstates,dim_psicoef
  integer(bit_kind), intent(in)  :: psidet(N_int,2,ndet)
  complex*16, intent(in)   :: psicoef(dim_psicoef,nstates)
  integer*8, allocatable         :: psi_det_save(:,:,:)
  complex*16, allocatable  :: psi_coef_save(:,:)

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

    allocate (psi_coef_save(ndet,nstates))
    do k=1,nstates
      do i=1,ndet
        psi_coef_save(i,k) = psicoef(i,k)
      enddo
      call normalize_complex(psi_coef_save(1,k),ndet)
    enddo

    call ezfio_set_determinants_psi_coef_complex(psi_coef_save)
    deallocate (psi_coef_save)

    allocate (psi_coef_save(ndet_qp_edit,nstates))
    do k=1,nstates
      do i=1,ndet_qp_edit
        psi_coef_save(i,k) = psicoef(i,k)
      enddo
      call normalize_complex(psi_coef_save(1,k),ndet_qp_edit)
    enddo

    call ezfio_set_determinants_psi_coef_complex_qp_edit(psi_coef_save)
    deallocate (psi_coef_save)

    call write_int(6,ndet,'Saved determinants')
  endif
end



subroutine save_wavefunction_specified_complex(ndet,nstates,psidet,psicoef,ndetsave,index_det_save)
  implicit none
  BEGIN_DOC
  !  Save the wave function into the |EZFIO| file
  END_DOC
  use bitmasks
  integer, intent(in)            :: ndet,nstates
  integer(bit_kind), intent(in)  :: psidet(N_int,2,ndet)
  complex*16, intent(in)   :: psicoef(ndet,nstates)
  integer, intent(in)            :: index_det_save(ndet)
  integer, intent(in)            :: ndetsave
  integer*8, allocatable         :: psi_det_save(:,:,:)
  complex*16, allocatable  :: psi_coef_save(:,:)
  integer*8                      :: det_8(100)
  integer(bit_kind)              :: det_bk((100*8)/bit_kind)
  integer                        :: N_int2
  equivalence (det_8, det_bk)

  integer                        :: i,j,k, ndet_qp_edit

  if (mpi_master) then
    ndet_qp_edit = min(ndetsave,N_det_qp_edit)
    call ezfio_set_determinants_N_int(N_int)
    call ezfio_set_determinants_bit_kind(bit_kind)
    call ezfio_set_determinants_N_det(ndetsave)
    call ezfio_set_determinants_N_det_qp_edit(ndet_qp_edit)
    call ezfio_set_determinants_n_states(nstates)
    call ezfio_set_determinants_mo_label(mo_label)

    N_int2 = (N_int*bit_kind)/8
    allocate (psi_det_save(N_int2,2,ndetsave))
    do i=1,ndetsave
      do k=1,N_int
        det_bk(k) = psidet(k,1,index_det_save(i))
      enddo
      do k=1,N_int2
        psi_det_save(k,1,i) = det_8(k)
      enddo
      do k=1,N_int
        det_bk(k) = psidet(k,2,index_det_save(i))
      enddo
      do k=1,N_int2
        psi_det_save(k,2,i) = det_8(k)
      enddo
    enddo
    call ezfio_set_determinants_psi_det(psi_det_save)
    call ezfio_set_determinants_psi_det_qp_edit(psi_det_save)
    deallocate (psi_det_save)

    allocate (psi_coef_save(ndetsave,nstates))
    double precision               :: accu_norm(nstates)
    accu_norm = 0.d0
    do k=1,nstates
      do i=1,ndetsave
        accu_norm(k) = accu_norm(k) + cdabs(psicoef(index_det_save(i),k) * psicoef(index_det_save(i),k))
        psi_coef_save(i,k) = psicoef(index_det_save(i),k)
      enddo
    enddo
    do k = 1, nstates
      accu_norm(k) = 1.d0/dsqrt(accu_norm(k))
    enddo
    do k=1,nstates
      do i=1,ndetsave
        psi_coef_save(i,k) = psi_coef_save(i,k) * accu_norm(k)
      enddo
    enddo

    call ezfio_set_determinants_psi_coef_complex(psi_coef_save)
    deallocate (psi_coef_save)

    allocate (psi_coef_save(ndet_qp_edit,nstates))
    accu_norm = 0.d0
    do k=1,nstates
      do i=1,ndet_qp_edit
        accu_norm(k) = accu_norm(k) + cdabs(psicoef(index_det_save(i),k) * psicoef(index_det_save(i),k))
        psi_coef_save(i,k) = psicoef(index_det_save(i),k)
      enddo
    enddo
    do k = 1, nstates
      accu_norm(k) = 1.d0/dsqrt(accu_norm(k))
    enddo
    do k=1,nstates
      do i=1,ndet_qp_edit
        psi_coef_save(i,k) = psi_coef_save(i,k) * accu_norm(k)
      enddo
    enddo
    !TODO: should this be psi_coef_complex_qp_edit?
    call ezfio_set_determinants_psi_coef_complex(psi_coef_save)
    deallocate (psi_coef_save)

    call write_int(6,ndet,'Saved determinants')
  endif
end

