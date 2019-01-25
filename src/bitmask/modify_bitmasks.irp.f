
use bitmasks
subroutine initialize_bitmask_to_restart_ones
 implicit none
 integer :: i,j,k,l,m
 integer :: ispin
  BEGIN_DOC
 ! Initialization of the generators_bitmask to the restart bitmask
  END_DOC
 do i = 1, N_int
   do k=1,N_generators_bitmask
     do ispin=1,2
       generators_bitmask(i,ispin,s_hole ,k) =   generators_bitmask_restart(i,ispin,s_hole ,k)
       generators_bitmask(i,ispin,s_part ,k) =   generators_bitmask_restart(i,ispin,s_part ,k)
       generators_bitmask(i,ispin,d_hole1,k) =   generators_bitmask_restart(i,ispin,d_hole1,k)
       generators_bitmask(i,ispin,d_part1,k) =   generators_bitmask_restart(i,ispin,d_part1,k)
       generators_bitmask(i,ispin,d_hole2,k) =   generators_bitmask_restart(i,ispin,d_hole2,k)
       generators_bitmask(i,ispin,d_part2,k) =   generators_bitmask_restart(i,ispin,d_part2,k)
     enddo
   enddo
 enddo
end


subroutine modify_bitmasks_for_hole(i_hole)
 implicit none
 integer, intent(in) :: i_hole
 integer :: i,j,k,l,m
 integer :: ispin
 BEGIN_DOC
! modify the generators_bitmask in order that one can only excite
! the electrons occupying i_hole
 END_DOC

 ! Set to Zero the holes
 do k=1,N_generators_bitmask
  do l = 1, 3
   i = index_holes_bitmask(l)
    do ispin=1,2
     do j = 1, N_int
      generators_bitmask(j,ispin,i,k) =   0_bit_kind
     enddo
    enddo
   enddo
 enddo

 k = shiftr(i_hole-1,bit_kind_shift)+1
 j = i_hole-shiftl(k-1,bit_kind_shift)-1
 do m = 1, N_generators_bitmask
  do l = 1, 3
   i = index_holes_bitmask(l)
   do ispin=1,2
    generators_bitmask(k,ispin,i,m) = ibset(generators_bitmask(k,ispin,i,m),j)
   enddo
  enddo
 enddo

end

subroutine modify_bitmasks_for_hole_in_out(i_hole)
 implicit none
 integer, intent(in) :: i_hole
 integer :: i,j,k,l,m
 integer :: ispin
 BEGIN_DOC
! modify the generators_bitmask in order that one can only excite
! the electrons occupying i_hole
 END_DOC

 k = shiftr(i_hole-1,bit_kind_shift)+1
 j = i_hole-shiftl(k-1,bit_kind_shift)-1
 do m = 1, N_generators_bitmask
  do l = 1, 3
   i = index_holes_bitmask(l)
   do ispin=1,2
    generators_bitmask(k,ispin,i,m) = ibset(generators_bitmask(k,ispin,i,m),j)
   enddo
  enddo
 enddo

end

subroutine modify_bitmasks_for_particl(i_part)
 implicit none
 integer, intent(in) :: i_part
 integer :: i,j,k,l,m
 integer :: ispin
 BEGIN_DOC
! modify the generators_bitmask in order that one can only excite
! the electrons to the orbital i_part
 END_DOC

 ! Set to Zero the particles
 do k=1,N_generators_bitmask
  do l = 1, 3
   i = index_particl_bitmask(l)
    do ispin=1,2
     do j = 1, N_int
      generators_bitmask(j,ispin,i,k) =   0_bit_kind
     enddo
    enddo
   enddo
 enddo

 k = shiftr(i_part-1,bit_kind_shift)+1
 j = i_part-shiftl(k-1,bit_kind_shift)-1
 do m = 1, N_generators_bitmask
  do l = 1, 3
   i = index_particl_bitmask(l)
   do ispin=1,2
    generators_bitmask(k,ispin,i,m) = ibset(generators_bitmask(k,ispin,i,m),j)
   enddo
  enddo
 enddo

end


subroutine set_bitmask_particl_as_input(input_bimask)
 implicit none
 integer(bit_kind), intent(in) :: input_bimask(N_int,2)
 integer :: i,j,k,l,m
 integer :: ispin
 BEGIN_DOC
! set the generators_bitmask for the particles
! as the input_bimask
 END_DOC

 do k=1,N_generators_bitmask
  do l = 1, 3
   i = index_particl_bitmask(l)
    do ispin=1,2
     do j = 1, N_int
      generators_bitmask(j,ispin,i,k) =  input_bimask(j,ispin)
     enddo
    enddo
   enddo
 enddo
 touch generators_bitmask

end


subroutine set_bitmask_hole_as_input(input_bimask)
 implicit none
 integer(bit_kind), intent(in) :: input_bimask(N_int,2)
 integer :: i,j,k,l,m
 integer :: ispin
 BEGIN_DOC
! set the generators_bitmask for the holes
! as the input_bimask
 END_DOC

 do k=1,N_generators_bitmask
  do l = 1, 3
   i = index_holes_bitmask(l)
    do ispin=1,2
     do j = 1, N_int
      generators_bitmask(j,ispin,i,k) =  input_bimask(j,ispin)
     enddo
    enddo
   enddo
 enddo
 touch generators_bitmask

end


subroutine print_generators_bitmasks_holes
 implicit none
 integer :: i,j,k,l
 integer(bit_kind),allocatable :: key_tmp(:,:)

 allocate(key_tmp(N_int,2))
 do l = 1, 3
  k = 1
  i = index_holes_bitmask(l)
   do j = 1, N_int
     key_tmp(j,1) = generators_bitmask(j,1,i,k)
     key_tmp(j,2) = generators_bitmask(j,2,i,k)
   enddo
   print*,''
   print*,'index hole  = ',i
   call print_det(key_tmp,N_int)
   print*,''
 enddo
 deallocate(key_tmp)

end

subroutine print_generators_bitmasks_particles
 implicit none
 integer :: i,j,k,l
 integer(bit_kind),allocatable :: key_tmp(:,:)

 allocate(key_tmp(N_int,2))
 do l = 1, 3
  k = 1
  i = index_particl_bitmask(l)
   do j = 1, N_int
     key_tmp(j,1) = generators_bitmask(j,1,i,k)
     key_tmp(j,2) = generators_bitmask(j,2,i,k)
   enddo
   print*,''
   print*,'index particl ',i
   call print_det(key_tmp,N_int)
   print*,''
 enddo
 deallocate(key_tmp)

end

subroutine print_generators_bitmasks_holes_for_one_generator(i_gen)
 implicit none
 integer, intent(in) :: i_gen
 integer :: i,j,k,l
 integer(bit_kind),allocatable :: key_tmp(:,:)

 allocate(key_tmp(N_int,2))
 do l = 1, 3
  k = i_gen
  i = index_holes_bitmask(l)
   do j = 1, N_int
     key_tmp(j,1) = generators_bitmask(j,1,i,k)
     key_tmp(j,2) = generators_bitmask(j,2,i,k)
   enddo
   print*,''
   print*,'index hole  = ',i
   call print_det(key_tmp,N_int)
   print*,''
 enddo
 deallocate(key_tmp)

end

subroutine print_generators_bitmasks_particles_for_one_generator(i_gen)
 implicit none
 integer, intent(in) :: i_gen
 integer :: i,j,k,l
 integer(bit_kind),allocatable :: key_tmp(:,:)

 allocate(key_tmp(N_int,2))
 do l = 1, 3
  k = i_gen
  i = index_particl_bitmask(l)
   do j = 1, N_int
     key_tmp(j,1) = generators_bitmask(j,1,i,k)
     key_tmp(j,2) = generators_bitmask(j,2,i,k)
   enddo
   print*,''
   print*,'index particl ',i
   call print_det(key_tmp,N_int)
   print*,''
 enddo
 deallocate(key_tmp)

end


 BEGIN_PROVIDER [integer,  index_holes_bitmask, (3)]
 implicit none
 BEGIN_DOC
! Index of the holes in the generators_bitmasks
 END_DOC
 index_holes_bitmask(1) = d_hole1
 index_holes_bitmask(2) = d_hole2
 index_holes_bitmask(3) = s_hole

 END_PROVIDER

 BEGIN_PROVIDER [integer,  index_particl_bitmask, (3)]
 implicit none
 BEGIN_DOC
! Index of the holes in the generators_bitmasks
 END_DOC
 index_particl_bitmask(1) = d_part1
 index_particl_bitmask(2) = d_part2
 index_particl_bitmask(3) = s_part

 END_PROVIDER
