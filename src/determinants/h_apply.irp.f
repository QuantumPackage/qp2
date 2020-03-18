use bitmasks
use omp_lib

type H_apply_buffer_type
  integer                        :: N_det
  integer                        :: sze
  integer(bit_kind), pointer     :: det(:,:,:)
  double precision , pointer     :: coef(:,:)
  double precision , pointer     :: e2(:,:)
end type H_apply_buffer_type

type(H_apply_buffer_type), pointer :: H_apply_buffer(:)


 BEGIN_PROVIDER [ logical, H_apply_buffer_allocated ]
&BEGIN_PROVIDER [ integer(omp_lock_kind), H_apply_buffer_lock, (64,0:nproc-1) ]
  use omp_lib
  implicit none
  BEGIN_DOC
  ! Buffer of determinants/coefficients/perturbative energy for H_apply.
  ! Uninitialized. Filled by H_apply subroutines.
  END_DOC
  integer                        :: iproc, sze
  sze = 10000
  if (.not.associated(H_apply_buffer)) then
    allocate(H_apply_buffer(0:nproc-1))
    iproc = 0
    !$OMP PARALLEL PRIVATE(iproc) DEFAULT(NONE)  &
    !$OMP SHARED(H_apply_buffer,N_int,sze,N_states,H_apply_buffer_lock)
    !$   iproc = omp_get_thread_num()
    H_apply_buffer(iproc)%N_det = 0
    H_apply_buffer(iproc)%sze = sze
    allocate (                                                       &
        H_apply_buffer(iproc)%det(N_int,2,sze),                      &
        H_apply_buffer(iproc)%coef(sze,N_states),                    &
        H_apply_buffer(iproc)%e2(sze,N_states)                       &
        )
    H_apply_buffer(iproc)%det  = 0_bit_kind
    H_apply_buffer(iproc)%coef = 0.d0
    H_apply_buffer(iproc)%e2   = 0.d0
    call omp_init_lock(H_apply_buffer_lock(1,iproc))
    !$OMP END PARALLEL
  endif
  do iproc=2,nproc-1
    if (.not.associated(H_apply_buffer(iproc)%det)) then
      print *,  ' ===================== Error =================== '
      print *,  'H_apply_buffer_allocated  should be provided outside'
      print *,  'of an OpenMP section'
      print *,  ' =============================================== '
      stop
    endif
  enddo

END_PROVIDER


subroutine resize_H_apply_buffer(new_size,iproc)
  implicit none
  integer, intent(in)            :: new_size, iproc
  integer(bit_kind), pointer     :: buffer_det(:,:,:)
  double precision,  pointer     :: buffer_coef(:,:)
  double precision,  pointer     :: buffer_e2(:,:)
  integer                        :: i,j,k
  integer                        :: Ndet

  BEGIN_DOC
! Resizes the H_apply buffer of proc iproc. The buffer lock should
! be set before calling this function.
  END_DOC
  PROVIDE H_apply_buffer_allocated

  ASSERT (new_size > 0)
  ASSERT (iproc >= 0)
  ASSERT (iproc < nproc)

  allocate ( buffer_det(N_int,2,new_size),                           &
      buffer_coef(new_size,N_states),                                &
      buffer_e2(new_size,N_states) )
  buffer_coef = 0.d0
  buffer_e2 = 0.d0
  do i=1,min(new_size,H_apply_buffer(iproc)%N_det)
    do k=1,N_int
      buffer_det(k,1,i) = H_apply_buffer(iproc)%det(k,1,i)
      buffer_det(k,2,i) = H_apply_buffer(iproc)%det(k,2,i)
    enddo
    ASSERT (sum(popcnt(H_apply_buffer(iproc)%det(:,1,i))) == elec_alpha_num)
    ASSERT (sum(popcnt(H_apply_buffer(iproc)%det(:,2,i))) == elec_beta_num )
  enddo
  deallocate(H_apply_buffer(iproc)%det)
  H_apply_buffer(iproc)%det => buffer_det

  do k=1,N_states
    do i=1,min(new_size,H_apply_buffer(iproc)%N_det)
      buffer_coef(i,k) = H_apply_buffer(iproc)%coef(i,k)
    enddo
  enddo
  deallocate(H_apply_buffer(iproc)%coef)
  H_apply_buffer(iproc)%coef => buffer_coef

  do k=1,N_states
    do i=1,min(new_size,H_apply_buffer(iproc)%N_det)
      buffer_e2(i,k)   = H_apply_buffer(iproc)%e2(i,k)
    enddo
  enddo
  deallocate(H_apply_buffer(iproc)%e2)
  H_apply_buffer(iproc)%e2 => buffer_e2

  H_apply_buffer(iproc)%sze = new_size
  H_apply_buffer(iproc)%N_det = min(new_size,H_apply_buffer(iproc)%N_det)

end

subroutine copy_H_apply_buffer_to_wf
  use omp_lib
  implicit none
  BEGIN_DOC
! Copies the H_apply buffer to psi_coef.
! After calling this subroutine, N_det, psi_det and psi_coef need to be touched
  END_DOC
  integer(bit_kind), allocatable :: buffer_det(:,:,:)
  double precision, allocatable  :: buffer_coef(:,:)
  integer                        :: i,j,k
  integer                        :: N_det_old

  PROVIDE H_apply_buffer_allocated


  ASSERT (N_int > 0)
  ASSERT (N_det > 0)

  allocate ( buffer_det(N_int,2,N_det), buffer_coef(N_det,N_states) )

  ! Backup determinants
  j=0
  do i=1,N_det
    if (pruned(i)) cycle  ! Pruned determinants
    j+=1
    ASSERT (sum(popcnt(psi_det(:,1,i))) == elec_alpha_num)
    ASSERT (sum(popcnt(psi_det(:,2,i))) == elec_beta_num)
    buffer_det(:,:,j) = psi_det(:,:,i)
  enddo
  N_det_old = j

  ! Backup coefficients
  do k=1,N_states
    j=0
    do i=1,N_det
      if (pruned(i)) cycle  ! Pruned determinants
      j += 1
      buffer_coef(j,k) = psi_coef(i,k)
    enddo
    ASSERT ( j == N_det_old )
  enddo

  ! Update N_det
  N_det = N_det_old
  do j=0,nproc-1
    N_det = N_det + H_apply_buffer(j)%N_det
  enddo

  ! Update array sizes
  if (psi_det_size < N_det) then
    psi_det_size = N_det
    TOUCH psi_det_size
  endif

  ! Restore backup in resized array
  do i=1,N_det_old
    psi_det(:,:,i) = buffer_det(:,:,i)
    ASSERT (sum(popcnt(psi_det(:,1,i))) == elec_alpha_num)
    ASSERT (sum(popcnt(psi_det(:,2,i))) == elec_beta_num )
  enddo
  do k=1,N_states
    do i=1,N_det_old
      psi_coef(i,k) = buffer_coef(i,k)
    enddo
  enddo

  ! Copy new buffers

  !$OMP PARALLEL DEFAULT(SHARED)                                     &
      !$OMP PRIVATE(j,k,i) FIRSTPRIVATE(N_det_old)                   &
      !$OMP SHARED(N_int,H_apply_buffer,psi_det,psi_coef,N_states,psi_det_size)
  j=0
  !$ j=omp_get_thread_num()
  do k=0,j-1
    N_det_old += H_apply_buffer(k)%N_det
  enddo
  do i=1,H_apply_buffer(j)%N_det
    do k=1,N_int
      psi_det(k,1,i+N_det_old) = H_apply_buffer(j)%det(k,1,i)
      psi_det(k,2,i+N_det_old) = H_apply_buffer(j)%det(k,2,i)
    enddo
    ASSERT (sum(popcnt(psi_det(:,1,i+N_det_old))) == elec_alpha_num)
    ASSERT (sum(popcnt(psi_det(:,2,i+N_det_old))) == elec_beta_num )
  enddo
  do k=1,N_states
    do i=1,H_apply_buffer(j)%N_det
      psi_coef(i+N_det_old,k) = H_apply_buffer(j)%coef(i,k)
    enddo
  enddo
  !$OMP BARRIER
  H_apply_buffer(j)%N_det = 0
  !$OMP END PARALLEL
  SOFT_TOUCH N_det psi_det psi_coef

  logical :: found_duplicates
  call remove_duplicates_in_psi_det(found_duplicates)
  do k=1,N_states
    call normalize(psi_coef(1,k),N_det)
  enddo
  SOFT_TOUCH N_det psi_det psi_coef

end

subroutine remove_duplicates_in_psi_det(found_duplicates)
  implicit none
  logical, intent(out) :: found_duplicates
  BEGIN_DOC
! Removes duplicate determinants in the wave function.
  END_DOC
  integer                        :: i,j,k
  integer(bit_kind), allocatable :: bit_tmp(:)
  logical,allocatable            :: duplicate(:)
  logical                        :: dup

  allocate (duplicate(N_det), bit_tmp(N_det))

  found_duplicates = .False.

  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,dup)

  !$OMP DO
  do i=1,N_det
    integer, external            :: det_search_key
    !$DIR FORCEINLINE
    bit_tmp(i) = det_search_key(psi_det_sorted_bit(1,1,i),N_int)
    duplicate(i) = .False.
  enddo
  !$OMP END DO

  !$OMP DO schedule(dynamic,1024)
  do i=1,N_det-1
    if (duplicate(i)) then
      cycle
    endif
    j = i+1
    do while (bit_tmp(j)==bit_tmp(i))
      if (duplicate(j)) then
        j = j+1
        if (j > N_det) then
          exit
        else
          cycle
        endif
      endif
      dup = .True.
      do k=1,N_int
        if ( (psi_det_sorted_bit(k,1,i) /= psi_det_sorted_bit(k,1,j) ) &
        .or. (psi_det_sorted_bit(k,2,i) /= psi_det_sorted_bit(k,2,j) ) ) then
          dup = .False.
          exit
        endif
      enddo
      if (dup) then
        duplicate(j) = .True.
        found_duplicates = .True.
      endif
      j += 1
      if (j > N_det) then
        exit
      endif
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  if (found_duplicates) then
    k=0
    do i=1,N_det
      if (.not.duplicate(i)) then
        k += 1
        psi_det(:,:,k) = psi_det_sorted_bit (:,:,i)
        psi_coef(k,:)  = psi_coef_sorted_bit(i,:)
      else
        if (sum(abs(psi_coef_sorted_bit(i,:))) /= 0.d0 ) then
          psi_coef(k,:)  = psi_coef_sorted_bit(i,:)
        endif
      endif
    enddo
    N_det = k
    psi_det_sorted_bit(:,:,1:N_det) = psi_det(:,:,1:N_det)
    psi_coef_sorted_bit(1:N_det,:) = psi_coef(1:N_det,:)
    TOUCH N_det psi_det psi_coef psi_det_sorted_bit psi_coef_sorted_bit c0_weight
  endif
  psi_det = psi_det_sorted
  psi_coef = psi_coef_sorted
  SOFT_TOUCH psi_det psi_coef psi_det_sorted_bit psi_coef_sorted_bit
  deallocate (duplicate,bit_tmp)
end


subroutine fill_H_apply_buffer_no_selection(n_selected,det_buffer,Nint,iproc)
  use bitmasks
  implicit none
  BEGIN_DOC
  !  Fill the H_apply buffer with determiants for |CISD|
  END_DOC

  integer, intent(in)            :: n_selected, Nint, iproc
  integer(bit_kind), intent(in)  :: det_buffer(Nint,2,n_selected)
  integer                        :: i,j,k
  integer                        :: new_size
  PROVIDE H_apply_buffer_allocated
  call omp_set_lock(H_apply_buffer_lock(1,iproc))
  new_size = H_apply_buffer(iproc)%N_det + n_selected
  if (new_size > H_apply_buffer(iproc)%sze) then
    call resize_h_apply_buffer(max(2*H_apply_buffer(iproc)%sze,new_size),iproc)
  endif
  do i=1,H_apply_buffer(iproc)%N_det
    ASSERT (sum(popcnt(H_apply_buffer(iproc)%det(:,1,i)) )== elec_alpha_num)
    ASSERT (sum(popcnt(H_apply_buffer(iproc)%det(:,2,i))) == elec_beta_num)
  enddo
  do i=1,n_selected
    do j=1,N_int
      H_apply_buffer(iproc)%det(j,1,i+H_apply_buffer(iproc)%N_det) = det_buffer(j,1,i)
      H_apply_buffer(iproc)%det(j,2,i+H_apply_buffer(iproc)%N_det) = det_buffer(j,2,i)
    enddo
    ASSERT (sum(popcnt(H_apply_buffer(iproc)%det(:,1,i+H_apply_buffer(iproc)%N_det)) )== elec_alpha_num)
    ASSERT (sum(popcnt(H_apply_buffer(iproc)%det(:,2,i+H_apply_buffer(iproc)%N_det))) == elec_beta_num)
  enddo
  do j=1,N_states
    do i=1,N_selected
      H_apply_buffer(iproc)%coef(i+H_apply_buffer(iproc)%N_det,j) = 0.d0
    enddo
  enddo
  H_apply_buffer(iproc)%N_det = new_size
  do i=1,H_apply_buffer(iproc)%N_det
    ASSERT (sum(popcnt(H_apply_buffer(iproc)%det(:,1,i)) )== elec_alpha_num)
    ASSERT (sum(popcnt(H_apply_buffer(iproc)%det(:,2,i))) == elec_beta_num)
  enddo
  call omp_unset_lock(H_apply_buffer_lock(1,iproc))
end

