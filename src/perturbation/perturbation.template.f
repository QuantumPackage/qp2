BEGIN_SHELL [ /usr/bin/env python2 ]
import perturbation
END_SHELL


subroutine perturb_buffer_$PERT(i_generator,buffer,buffer_size,e_2_pert_buffer,coef_pert_buffer,sum_e_2_pert,sum_norm_pert,sum_H_pert_diag,N_st,Nint,key_mask,fock_diag_tmp,electronic_energy)
  implicit none
  BEGIN_DOC
  ! Apply pertubration ``$PERT`` to the buffer of determinants generated in the H_apply
  ! routine.
  END_DOC

  integer, intent(in)            :: Nint, N_st, buffer_size, i_generator
  integer(bit_kind), intent(in)  :: buffer(Nint,2,buffer_size)
  integer(bit_kind),intent(in)    :: key_mask(Nint,2)
  double precision, intent(in)    :: fock_diag_tmp(2,0:mo_num)
  double precision, intent(in)    :: electronic_energy(N_st)
  double precision, intent(inout) :: sum_norm_pert(N_st),sum_e_2_pert(N_st)
  double precision, intent(inout) :: coef_pert_buffer(N_st,buffer_size),e_2_pert_buffer(N_st,buffer_size),sum_H_pert_diag(N_st)
  double precision               :: c_pert(N_st), e_2_pert(N_st),  H_pert_diag(N_st)
  integer                        :: i,k,l, c_ref, ni, ex
  integer, external              :: connected_to_ref
  logical, external              :: is_in_wavefunction

  integer(bit_kind), allocatable :: minilist(:,:,:)
  integer, allocatable           :: idx_minilist(:)
  integer                        :: N_minilist

  integer(bit_kind), allocatable :: minilist_gen(:,:,:)
  integer :: N_minilist_gen
  logical :: fullMatch
  logical, external :: is_connected_to

  integer(bit_kind), allocatable :: microlist(:,:,:), microlist_zero(:,:,:)
  integer, allocatable           :: idx_microlist(:), N_microlist(:), ptr_microlist(:), idx_microlist_zero(:)
  integer :: mobiles(2), smallerlist


  integer(bit_kind), allocatable :: microlist_gen(:,:,:)
  integer, allocatable           :: idx_microlist_gen(:), N_microlist_gen(:), ptr_microlist_gen(:)

  allocate( minilist(Nint,2,N_det_selectors),                        &
      minilist_gen(Nint,2,N_det_generators),                         &
      idx_minilist(N_det_selectors))



  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)
  ASSERT (buffer_size >= 0)
  ASSERT (minval(sum_norm_pert) >= 0.d0)
  ASSERT (N_st > 0)


  call create_minilist_find_previous(key_mask, psi_det_generators, miniList_gen, i_generator-1, N_minilist_gen, fullMatch, Nint)


  if(fullMatch) then
    deallocate( minilist, minilist_gen, idx_minilist )
    return
  end if

  call create_minilist(key_mask, psi_selectors, minilist, idx_miniList, N_det_selectors, N_minilist, Nint)
  allocate(   microlist(Nint,2,N_minilist*4),               &
       idx_microlist(N_minilist*4),                  &
       ptr_microlist(0:mo_num*2+1),  &
       N_microlist(0:mo_num*2) )

  allocate(   microlist_gen(Nint,2,N_minilist_gen*4),               &
      idx_microlist_gen(N_minilist_gen*4 ),                  &
      ptr_microlist_gen(0:mo_num*2+1),  &
       N_microlist_gen(0:mo_num*2) )

  if(key_mask(1,1) /= 0) then

    call create_microlist(minilist, N_minilist, key_mask, microlist, idx_microlist, N_microlist, ptr_microlist, Nint)
    call create_microlist(minilist_gen, N_minilist_gen, key_mask, microlist_gen, idx_microlist_gen, N_microlist_gen,ptr_microlist_gen,Nint)

    allocate(microlist_zero(Nint,2,N_minilist))
    allocate(idx_microlist_zero(N_minilist))


    do i=0,mo_num*2
      do k=ptr_microlist(i),ptr_microlist(i+1)-1
        idx_microlist(k) = idx_minilist(idx_microlist(k))
      end do
    end do


    if(N_microlist(0) > 0) then
! TODO OLD
!      microlist_zero(:,:,1:N_microlist(0)) = microlist(:,:,1:N_microlist(0))
!      idx_microlist_zero(1:N_microlist(0)) = idx_microlist(1:N_microlist(0))
! TODO OLD
      ASSERT (N_microlist(0) <= N_minilist)
      do l=1,N_microlist(0)
        do k=1,Nint
          microlist_zero(k,1,l) = microlist(k,1,l)
          microlist_zero(k,2,l) = microlist(k,2,l)
        enddo
        idx_microlist_zero(l) = idx_microlist(l)
      enddo
    end if

  end if

  do i=1,buffer_size

    if (is_in_wavefunction(buffer(1,1,i),Nint)) then
      cycle
    endif

    if(key_mask(1,1) /= 0) then
      call getMobiles(buffer(1,1,i), key_mask, mobiles, Nint)
      if(N_microlist(mobiles(1)) < N_microlist(mobiles(2))) then
        smallerlist = mobiles(1)
      else
        smallerlist = mobiles(2)
      end if

      if(N_microlist_gen(smallerlist) > 0) then
! TODO OLD
!        if(is_connected_to(buffer(1,1,i), microlist_gen(:,:,ptr_microlist_gen(smallerlist):ptr_microlist_gen(smallerlist+1)-1), Nint, N_microlist_gen(smallerlist))) then
! TODO OLD
      ASSERT (ptr_microlist_gen(smallerlist) <= N_minilist_gen*4)
        if(is_connected_to(buffer(1,1,i), microlist_gen(1,1,ptr_microlist_gen(smallerlist)), Nint, N_microlist_gen(smallerlist))) then
          cycle
        end if
      end if
      if(N_microlist_gen(0) > 0) then
! TODO OLD
!        if(is_connected_to(buffer(1,1,i), microlist_gen(:,:,1:ptr_microlist_gen(1)-1), Nint, N_microlist_gen(0))) then
! TODO OLD
        if(is_connected_to(buffer(1,1,i), microlist_gen(1,1,1), Nint, N_microlist_gen(0))) then
          cycle
        end if
      end if

      if(N_microlist(smallerlist) > 0) then
! TODO OLD
!         microlist_zero(:,:,ptr_microlist(1):ptr_microlist(1)+N_microlist(smallerlist)-1) = microlist(:,:,ptr_microlist(smallerlist):ptr_microlist(smallerlist+1)-1)
!         idx_microlist_zero(ptr_microlist(1):ptr_microlist(1)+N_microlist(smallerlist)-1) = idx_microlist(ptr_microlist(smallerlist):ptr_microlist(smallerlist+1)-1)
! TODO OLD
        ASSERT ( ptr_microlist(1)+N_microlist(smallerlist)-1 <= N_minilist )
        ASSERT ( ptr_microlist(smallerlist)+N_microlist(smallerlist)-1 <= N_minilist*4 )
        do l=0, N_microlist(smallerlist)-1
          do k=1,Nint
            microlist_zero(k,1,ptr_microlist(1)+l) = microlist(k,1,ptr_microlist(smallerlist)+l)
            microlist_zero(k,2,ptr_microlist(1)+l) = microlist(k,2,ptr_microlist(smallerlist)+l)
          enddo
          idx_microlist_zero(ptr_microlist(1)+l) = idx_microlist(ptr_microlist(smallerlist)+l)
        enddo
      end if
      call pt2_$PERT(electronic_energy,psi_det_generators(1,1,i_generator),buffer(1,1,i), fock_diag_tmp,     &
          c_pert,e_2_pert,H_pert_diag,Nint,N_microlist(smallerlist)+N_microlist(0),        &
          n_st,microlist_zero,idx_microlist_zero,N_microlist(smallerlist)+N_microlist(0))
    else
      ASSERT (N_minilist_gen <= N_det_generators)
      if(is_connected_to(buffer(1,1,i), miniList_gen, Nint, N_minilist_gen)) then
        cycle
      end if

      call pt2_$PERT(electronic_energy,psi_det_generators(1,1,i_generator),buffer(1,1,i), fock_diag_tmp,        &
        c_pert,e_2_pert,H_pert_diag,Nint,N_minilist,n_st,minilist,idx_minilist,N_minilist)
    end if

!     call pt2_$PERT(electronic_energy,psi_det_generators(1,1,i_generator),buffer(1,1,i), fock_diag_tmp,        &
!          c_pert,e_2_pert,H_pert_diag,Nint,N_minilist,n_st,minilist,idx_minilist,N_minilist)

    do k = 1,N_st
      e_2_pert_buffer(k,i)   = e_2_pert(k)
      coef_pert_buffer(k,i)  = c_pert(k)
      sum_norm_pert(k)       = sum_norm_pert(k)   + c_pert(k) * c_pert(k)
      sum_e_2_pert(k)        = sum_e_2_pert(k)    + e_2_pert(k)
      sum_H_pert_diag(k)     = sum_H_pert_diag(k) + H_pert_diag(k)
    enddo

  enddo
  deallocate( minilist, minilist_gen, idx_minilist,                  &
      microlist, idx_microlist, N_microlist,ptr_microlist,           &
      microlist_gen, idx_microlist_gen,N_microlist_gen,ptr_microlist_gen )
end


subroutine perturb_buffer_by_mono_$PERT(i_generator,buffer,buffer_size,e_2_pert_buffer,coef_pert_buffer,sum_e_2_pert,sum_norm_pert,sum_H_pert_diag,N_st,Nint,key_mask,fock_diag_tmp,electronic_energy)
  implicit none
  BEGIN_DOC
  ! Apply pertubration ``$PERT`` to the buffer of determinants generated in the H_apply
  ! routine.
  END_DOC

  integer, intent(in)            :: Nint, N_st, buffer_size, i_generator
  integer(bit_kind), intent(in)  :: buffer(Nint,2,buffer_size)
  integer(bit_kind),intent(in)    :: key_mask(Nint,2)
  double precision, intent(in)    :: fock_diag_tmp(2,0:mo_num)
  double precision, intent(in)    :: electronic_energy(N_st)
  double precision, intent(inout) :: sum_norm_pert(N_st),sum_e_2_pert(N_st)
  double precision, intent(inout) :: coef_pert_buffer(N_st,buffer_size),e_2_pert_buffer(N_st,buffer_size),sum_H_pert_diag(N_st)
  double precision               :: c_pert(N_st), e_2_pert(N_st),  H_pert_diag(N_st)
  integer                        :: i,k, c_ref, ni, ex
  integer, external              :: connected_to_ref_by_single
  logical, external              :: is_in_wavefunction

  integer(bit_kind), allocatable :: minilist(:,:,:)
  integer, allocatable           :: idx_minilist(:)
  integer                        :: N_minilist

  integer(bit_kind), allocatable :: minilist_gen(:,:,:)
  integer :: N_minilist_gen
  logical :: fullMatch
  logical, external :: is_connected_to

  allocate( minilist(Nint,2,N_det_selectors),                        &
      minilist_gen(Nint,2,N_det_generators),                         &
      idx_minilist(N_det_selectors) )


  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)
  ASSERT (buffer_size >= 0)
  ASSERT (minval(sum_norm_pert) >= 0.d0)
  ASSERT (N_st > 0)

  call create_minilist(key_mask, psi_selectors, miniList, idx_miniList, N_det_selectors, N_minilist, Nint)
  call create_minilist_find_previous(key_mask, psi_det_generators, miniList_gen, i_generator-1, N_minilist_gen, fullMatch, Nint)

  if(fullMatch) then
    deallocate( minilist, minilist_gen, idx_minilist )
    return
  end if


  do i=1,buffer_size

    c_ref = connected_to_ref_by_single(buffer(1,1,i),psi_det_generators,Nint,i_generator,N_det)

    if (c_ref /= 0) then
      cycle
    endif

    if (is_in_wavefunction(buffer(1,1,i),Nint)) then
      cycle
    endif

    call pt2_$PERT(electronic_energy,psi_det_generators(1,1,i_generator),buffer(1,1,i), fock_diag_tmp,        &
         c_pert,e_2_pert,H_pert_diag,Nint,N_minilist,n_st,minilist,idx_minilist,N_minilist)

    do k = 1,N_st
      e_2_pert_buffer(k,i)   = e_2_pert(k)
      coef_pert_buffer(k,i)  = c_pert(k)
      sum_norm_pert(k)       = sum_norm_pert(k)   + c_pert(k) * c_pert(k)
      sum_e_2_pert(k)        = sum_e_2_pert(k)    + e_2_pert(k)
      sum_H_pert_diag(k)     = sum_H_pert_diag(k) + H_pert_diag(k)
    enddo

  enddo
  deallocate( minilist, minilist_gen, idx_minilist )

end

