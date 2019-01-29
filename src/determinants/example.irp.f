subroutine example_determinants
  use bitmasks ! you need to include the bitmasks_module.f90 features
  implicit none
  BEGIN_DOC
  ! subroutine that illustrates the main features available in determinants
  END_DOC
  print*,'a determinant is stored as a binary representation of the occupancy of the spatial orbitals'
  print*,'see the bitmask module for more information about that '
  print*,'a spin determinant is an array of (N_int) integers of type bit_kind (see bitmask for more information)'
  print*,'A determinant containing alpha and beta electrons is an array of dimension (2,N_int)'
  integer(bit_kind), allocatable :: det_i(:,:)
  allocate(det_i(N_int,2))
  print*,'det_i(1,:) alpha spins '
  print*,'det_i(2,:) beta  spins '
  integer                        :: i,j
  print*,'initialize det_i to an electron occupation corresponding RHF or ROHF: ref_bitmask '
  do i = 1, N_int
    det_i(i,1) = ref_bitmask(i,1)
    det_i(i,2) = ref_bitmask(i,2)
  enddo
  print*,''
  print*,'print a human readable representation of the determinant '
  call print_det(det_i,N_int)
  print*,'doing a single excitation on top of det_i'
  integer                        :: h1,p1,s1,i_ok
  h1 = 1
  p1 = elec_alpha_num + 1
  s1 = 1
  print*,'h1 --> p1 of spin s1'
  print*,'i_ok == +1 : excitation is possible '
  print*,'i_ok == -1 : excitation is NOT possible '
  call do_mono_excitation(det_i,h1,p1,s1,i_ok)
  print*,'h1,p1,s1,i_ok'
  print*, h1,p1,s1,i_ok
  if(i_ok == -1)then
    print*,'excitation was not possible '
    stop
  endif
  call debug_det(det_i,N_int)
  print*,'computing the interaction between ref_determinant and det_i '
  double precision               :: h0i,hii,h00
  call i_H_j(det_i,det_i,N_int,h0i)
  print*,'  < ref | H | det_i > = ',h0i
  print*,'computing the diagonal Hamiltonian matrix element of det_i '
  call i_H_j(ref_bitmask,det_i,N_int,hii)
  print*,'< det_i | H | det_i > = ',hii
  print*,'computing the first-order coefficient of det_i with H0=EN '
  double precision               :: c_i
  call i_H_j(ref_bitmask,ref_bitmask,N_int,h00)
  c_i = h0i/(h00 - hii)
  print*,'c_i^{(1)} = ',c_i
  print*,''
  print*,'doing another single excitation on top of det_i'
  h1 = elec_alpha_num
  p1 = elec_alpha_num + 1
  s1 = 2
  call do_mono_excitation(det_i,h1,p1,s1,i_ok)
  print*,'h1,p1,s1,i_ok'
  print*, h1,p1,s1,i_ok
  call i_H_j(det_i,det_i,N_int,h0i)
  print*,'  < ref | H | det_i > = ',h0i
  print*,'computing the diagonal Hamiltonian matrix element of det_i '
  call i_H_j(ref_bitmask,ref_bitmask,N_int,h00)
  c_i = h0i/(h00 - hii)
  print*,'c_i^{(1)} = ',c_i
  print*,''
  print*,'Finding the excitation degree between two arbitrary determinants '
  integer                        :: exc(0:2,2,2)
  double precision               :: phase
  integer                        :: h2,p2,s2,degree
  call get_excitation_degree(ref_bitmask,det_i,degree,N_int)
  print*,'degree = ',degree
  print*,'Finding the differences in terms of holes and particles, together with the fermionic phase '
  call get_excitation(ref_bitmask,det_i,exc,degree,phase,N_int)
  print*,'Fermionic phase for the excitation from ref_bitmask to det_i'
  print*,phase
  print*,'put the excitation information in a human readable format'
  call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
  print*,'s1',s1
  print*,'h1,p1 = ',h1,p1
  print*,'s2',s2
  print*,'h2,p2 = ',h2,p2
  print*,''
  print*,'Finding the occupancy of det_i'
  integer, allocatable           :: occ(:,:)
  integer                        :: n_occ_ab(2)
  allocate(occ(N_int*bit_kind_size,2))
  call bitstring_to_list_ab(det_i, occ, n_occ_ab, N_int)
  print*,'alpha electrons orbital occupancy'
  do i = 1, n_occ_ab(1) ! browsing the alpha electrons
    print*,occ(i,1)
  enddo
  print*,'beta  electrons orbital occupancy'
  do i = 1, n_occ_ab(2) ! browsing the beta  electrons
    print*,occ(i,2)
  enddo
end


subroutine example_determinants_psi_det
  use bitmasks ! you need to include the bitmasks_module.f90 features
  implicit none
  BEGIN_DOC
  ! subroutine that illustrates the main features available in determinants using the psi_det/psi_coef
  END_DOC
  read_wf = .True.
  touch read_wf
  ! you force the wave function to be set to the one in the EZFIO directory
  call routine_example_psi_det
end

subroutine routine_example_psi_det
  use bitmasks ! you need to include the bitmasks_module.f90 features
  implicit none
  BEGIN_DOC
  ! subroutine that illustrates the main features available in determinants using many determinants
  END_DOC
  integer                        :: i,j
  integer, allocatable           :: degree_list(:)
  integer, allocatable           :: idx(:)
  allocate(degree_list(N_det),idx(0:N_det))

  print*,'Number of determinants in the wave function'
  print*,'N_det = ',N_det
  print*,''
  print*,'Printing in a human readable format all Slater determinants '
  do i = 1, N_det
    call debug_det(psi_det(1,1,i),N_int)
  enddo
  print*,''
  print*,'Number of states computed '
  print*,'N_states = ',N_states
  print*,'Printing the coefficients for all states for all Slater determinants '
  do j = 1, N_states
    print*,'State = ',j
    do i = 1, N_det
      write(*,'(I9,X,F16.10)')i,psi_coef(i,j)
    enddo
  enddo
  print*,''
  print*,'Finding the connection through a two-electron operator in the wave function'
  print*,'You want to know the connections of the first determinant '
  !                                 wave function  determinant    exc degree              list
  call get_excitation_degree_vector(   psi_det  ,  psi_det(1,1,1),degree_list,N_int,N_det,idx)
  double precision               :: hij
  double precision, allocatable  :: i_H_psi(:)
  allocate(i_H_psi(N_states))
  i_H_psi = 0.d0
  print*,'Computing <psi_det(1) | H | psi_det > = \sum_I c_I <psi_det(1)| H | psi_det(I)>'
  do i = 1, idx(0) ! number of Slater determinants connected to the first one
    print*,'Determinant connected'
    call debug_det(psi_det(1,1,idx(i)),N_int)
    print*,'excitation degree = ',degree_list(i)
    call i_H_j(psi_det(1,1,1) , psi_det(1,1,idx(i)),hij,N_int)
    do j = 1, N_states
      i_H_psi(j) += hij * psi_coef(idx(i),j)
    enddo
  enddo
  print*,'i_H_psi = ',i_H_psi
end





