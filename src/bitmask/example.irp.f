subroutine example_bitmask
 use bitmasks ! you need to include the bitmasks_module.f90 features
 implicit none
 BEGIN_DOC
! subroutine that illustrates the main features available in bitmask
 END_DOC
 integer :: i,j
 print*,''
 print*,'**************'
 print*,'**************'
 print*,'MO class: to set the various type of MO class, see the following exectuable'
 print*,'qp_set_mo_class'
 print*,'**************'
 print*,'number of core orbitals = ',n_core_orb
 print*,'list of the core orbitals '
 do i = 1, n_core_orb
  write(*,'(2(I3,X))')i,list_core(i)
 enddo

 print*,'number of inact orbitals = ',n_inact_orb
 print*,'list of the inact orbitals '
 do i = 1, n_inact_orb
  write(*,'(2(I3,X))')i,list_inact(i)
 enddo

 print*,'number of act orbitals = ',n_act_orb
 print*,'list of the act orbitals '
 do i = 1, n_act_orb
  write(*,'(2(I3,X))')i,list_act(i)
 enddo

 print*,'number of virt orbitals = ',n_virt_orb
 print*,'list of the virt orbitals '
 do i = 1, n_virt_orb
  write(*,'(2(I3,X))')i,list_virt(i)
 enddo
 print*,''
 print*,'**************'
 print*,'**************'
 print*,'manipulating bitstrings (usefull for determinant representation)'
 print*,'**************'
 integer(bit_kind), allocatable :: key(:)
 print*,'Size of the integers used to represent all the orbitals '
 print*,'bit_kind = ',bit_kind
 print*,'Number of bits in the integers ',bit_kind_size
 print*,'Number of integers to represent all the orbitals on integer'
 print*,'N_int = ',N_int
 allocate(key(N_int))
 print*,'**** '
 print*,' initialize a bistring to zero '
 do i = 1, N_int
  key(i) = 0_bit_kind
 enddo
 print*,'print a human readable representation of the bitstring'
 call bitstring_to_str( output, key, N_int )
 print *,  trim(output)
 integer :: i_orb
 character*(2048) :: output

 do i_orb = 1, min(4,mo_num) ! you set the first four bits to 1 in key
  call set_bit_to_integer(i_orb,key,N_int)
 enddo
 print*,'print a human readable representation of the bitstring'
 call bitstring_to_str( output, key, N_int )
 print *,  trim(output)
 print*,''
 integer :: n_elements
 integer, allocatable :: list_occ(:)
 allocate(list_occ(N_int*bit_kind_size))
 call bitstring_to_list( key, list_occ, n_elements, N_int)
 print*,'number of bits set to 1 = ',n_elements
 print*,'list of bits set to 1 '
 do i = 1, n_elements
  write(*,'(2(I3,X))')i,list_occ(i)
 enddo
 call clear_bit_to_integer(2,key,N_int) ! you set to 0 the second bit
 print*,'print a human readable representation of the bitstring'
 call bitstring_to_str( output, key, N_int )
 print *,  trim(output)

end
