use bitmasks

 BEGIN_PROVIDER [integer(bit_kind), holes_operators, (N_int,2)]
&BEGIN_PROVIDER [integer(bit_kind), particles_operators, (N_int,2)]

 BEGIN_DOC
 ! holes_operators represents an array of integers where all the holes have
 ! been done going from psi_ref to psi_non_ref
 ! particles_operators represents an array of integers where all the particles have
 ! been done going from psi_ref to psi_non_ref
 END_DOC
 holes_operators = 0_bit_kind
 particles_operators = 0_bit_kind
 implicit none
 integer(bit_kind), allocatable :: key_test(:,:)
 integer(bit_kind), allocatable :: holes(:,:),particles(:,:)
 allocate(key_test(N_int,2))
 allocate(holes(N_int,2),particles(N_int,2))
 integer :: i,j,k
 print*,'providing holes_operators and particles_operators'
 do i = 1, N_det_ref
  do j = 1, N_det_non_ref
   do k = 1, N_int
    key_test(k,1) = xor(psi_ref(k,1,i),psi_non_ref(k,1,j))
    key_test(k,2) = xor(psi_ref(k,2,i),psi_non_ref(k,2,j))
   enddo
   do k = 1,N_int
    holes(k,1) = iand(psi_ref(k,1,i),key_test(k,1))
    holes(k,2) = iand(psi_ref(k,2,i),key_test(k,2))
    particles(k,1) = iand(psi_non_ref(k,1,j),key_test(k,1))
    particles(k,2) = iand(psi_non_ref(k,2,j),key_test(k,2))
   enddo
   do k = 1, N_int
    holes_operators(k,1) = ior(holes_operators(k,1),holes(k,1))
    holes_operators(k,2) = ior(holes_operators(k,2),holes(k,2))
    particles_operators(k,1) = ior(particles_operators(k,1),particles(k,1))
    particles_operators(k,2) = ior(particles_operators(k,2),particles(k,2))
   enddo
  enddo
 enddo

 deallocate(key_test)
 deallocate(holes,particles)

END_PROVIDER
