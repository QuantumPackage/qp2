
use bitmasks

subroutine copy_H_apply_buffer_to_wf_tc
  use omp_lib
  implicit none
  BEGIN_DOC
! Copies the H_apply buffer to psi_coef.
! After calling this subroutine, N_det, psi_det and psi_coef need to be touched
  END_DOC
  integer(bit_kind), allocatable :: buffer_det(:,:,:)
  double precision, allocatable  :: buffer_r_coef(:,:), buffer_l_coef(:,:)
  integer                        :: i,j,k
  integer                        :: N_det_old

  PROVIDE H_apply_buffer_allocated


  ASSERT (N_int > 0)
  ASSERT (N_det > 0)

  allocate ( buffer_det(N_int,2,N_det), buffer_r_coef(N_det,N_states), buffer_l_coef(N_det,N_states)  )

  ! Backup determinants
  j=0
  do i=1,N_det
!    if (pruned(i)) cycle  ! Pruned determinants
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
!      if (pruned(i)) cycle  ! Pruned determinants
      j += 1
      buffer_r_coef(j,k) = psi_r_coef_bi_ortho(i,k)
      buffer_l_coef(j,k) = psi_l_coef_bi_ortho(i,k)
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
      psi_r_coef_bi_ortho(i,k) = buffer_r_coef(i,k)
      psi_l_coef_bi_ortho(i,k) = buffer_l_coef(i,k)
    enddo
  enddo

  ! Copy new buffers

  !$OMP PARALLEL DEFAULT(SHARED)                                     &
      !$OMP PRIVATE(j,k,i) FIRSTPRIVATE(N_det_old)                   &
      !$OMP SHARED(N_int,H_apply_buffer,psi_det,psi_r_coef_bi_ortho,psi_l_coef_bi_ortho,N_states,psi_det_size)
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
      psi_r_coef_bi_ortho(i+N_det_old,k) = H_apply_buffer(j)%coef(i,k)
      psi_l_coef_bi_ortho(i+N_det_old,k) = 0.d0
    enddo
  enddo
  !$OMP BARRIER
  H_apply_buffer(j)%N_det = 0
  !$OMP END PARALLEL
  SOFT_TOUCH N_det psi_det psi_r_coef_bi_ortho psi_l_coef_bi_ortho 

  logical :: found_duplicates
  call remove_duplicates_in_psi_det_tc(found_duplicates)
  call bi_normalize(psi_l_coef_bi_ortho,psi_r_coef_bi_ortho,size(psi_l_coef_bi_ortho,1),N_det,N_states)
  SOFT_TOUCH N_det psi_det psi_r_coef_bi_ortho psi_l_coef_bi_ortho

end

subroutine remove_duplicates_in_psi_det_tc(found_duplicates)
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
    bit_tmp(i) = det_search_key(psi_det_sorted_bit_tc(1,1,i),N_int)
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
        if ( (psi_det_sorted_bit_tc(k,1,i) /= psi_det_sorted_bit_tc(k,1,j) ) &
        .or. (psi_det_sorted_bit_tc(k,2,i) /= psi_det_sorted_bit_tc(k,2,j) ) ) then
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
        psi_det(:,:,k) = psi_det_sorted_bit_tc (:,:,i)
        psi_r_coef_bi_ortho(k,:)  = psi_r_coef_sorted_bit(i,:)
        psi_l_coef_bi_ortho(k,:)  = psi_l_coef_sorted_bit(i,:)
      else
        if (sum(abs(psi_r_coef_sorted_bit(i,:))) /= 0.d0 ) then
          psi_r_coef_bi_ortho(k,:)  = psi_r_coef_sorted_bit(i,:)
          psi_l_coef_bi_ortho(k,:)  = psi_l_coef_sorted_bit(i,:)
        endif
      endif
    enddo
    N_det = k
    psi_det_sorted_bit_tc(:,:,1:N_det) = psi_det(:,:,1:N_det)
    psi_r_coef_sorted_bit(1:N_det,:) = psi_r_coef_bi_ortho(1:N_det,:)
    psi_l_coef_sorted_bit(1:N_det,:) = psi_l_coef_bi_ortho(1:N_det,:)
    TOUCH N_det psi_det psi_det_sorted_bit_tc c0_weight psi_r_coef_sorted_bit psi_l_coef_sorted_bit
  endif
  psi_det = psi_det_sorted_tc
  psi_r_coef_bi_ortho = psi_r_coef_sorted_bi_ortho
  psi_l_coef_bi_ortho = psi_l_coef_sorted_bi_ortho
  SOFT_TOUCH psi_det psi_r_coef_bi_ortho psi_l_coef_bi_ortho psi_det_sorted_bit_tc psi_r_coef_sorted_bit psi_l_coef_sorted_bit
  deallocate (duplicate,bit_tmp)
end


 BEGIN_PROVIDER [ integer(bit_kind), psi_det_sorted_bit_tc, (N_int,2,psi_det_size) ]
&BEGIN_PROVIDER [ double precision, psi_r_coef_sorted_bit, (N_det,N_states) ]
&BEGIN_PROVIDER [ double precision, psi_l_coef_sorted_bit, (N_det,N_states) ]
   implicit none
   BEGIN_DOC
   ! Determinants on which we apply $\langle i|H|psi \rangle$ for perturbation.
   ! They are sorted by determinants interpreted as integers. Useful
   ! to accelerate the search of a random determinant in the wave
   ! function.
   END_DOC

   call sort_dets_by_det_search_key(N_det, psi_det, psi_r_coef_bi_ortho, size(psi_r_coef_bi_ortho,1),       &
       psi_det_sorted_bit_tc, psi_r_coef_sorted_bit, N_states)
   call sort_dets_by_det_search_key(N_det, psi_det, psi_l_coef_bi_ortho, size(psi_l_coef_bi_ortho,1),       &
       psi_det_sorted_bit_tc, psi_l_coef_sorted_bit, N_states)

END_PROVIDER
