use bitmasks

BEGIN_PROVIDER [ integer, N_det_generators ]
 implicit none
 BEGIN_DOC
 ! For Single reference wave functions, the number of generators is 1 : the
 ! Hartree-Fock determinant
 END_DOC
 integer :: i
 double precision :: norm
 call write_time(6)
 norm = 1.d0
 N_det_generators = N_det
 do i=1,N_det
   norm = norm - psi_average_norm_contrib_sorted_tc(i)
   if (norm - 1.d-10 < 1.d0 - threshold_generators) then
     N_det_generators = i
     exit
   endif
 enddo
 N_det_generators = max(N_det_generators,1)
 call write_int(6,N_det_generators,'Number of generators')
END_PROVIDER

 BEGIN_PROVIDER [ integer(bit_kind), psi_det_generators, (N_int,2,psi_det_size) ]
&BEGIN_PROVIDER [ double precision, psi_coef_generators, (psi_det_size,N_states) ]
 implicit none
 BEGIN_DOC
 ! For Single reference wave functions, the generator is the
 ! Hartree-Fock determinant
 END_DOC
 psi_det_generators(1:N_int,1:2,1:N_det) = psi_det_sorted_tc(1:N_int,1:2,1:N_det)
 psi_coef_generators(1:N_det,1:N_states) = psi_coef_sorted_tc(1:N_det,1:N_states)

END_PROVIDER

 BEGIN_PROVIDER [ integer(bit_kind), psi_det_sorted_gen, (N_int,2,psi_det_size) ]
&BEGIN_PROVIDER [ double precision, psi_coef_sorted_gen, (psi_det_size,N_states) ]
&BEGIN_PROVIDER [ integer, psi_det_sorted_gen_order,     (psi_det_size)  ]

 implicit none
 BEGIN_DOC
 ! For Single reference wave functions, the generator is the
 ! Hartree-Fock determinant
 END_DOC
 psi_det_sorted_gen =  psi_det_sorted_tc
 psi_coef_sorted_gen = psi_coef_sorted_tc
 psi_det_sorted_gen_order = psi_det_sorted_tc_order
END_PROVIDER


BEGIN_PROVIDER [integer, degree_max_generators]
 implicit none
 BEGIN_DOC
! Max degree of excitation (respect to HF) of the generators
 END_DOC
 integer :: i,degree
  degree_max_generators = 0
  do i = 1, N_det_generators
   call get_excitation_degree(HF_bitmask,psi_det_generators(1,1,i),degree,N_int)
   if(degree .gt. degree_max_generators)then
    degree_max_generators = degree
   endif
  enddo
END_PROVIDER

BEGIN_PROVIDER [ integer, size_select_max]
 implicit none
 BEGIN_DOC
 ! Size of the select_max array
 END_DOC
 size_select_max = 10000
END_PROVIDER

BEGIN_PROVIDER [ double precision, select_max, (size_select_max) ]
 implicit none
 BEGIN_DOC
 ! Memo to skip useless selectors
 END_DOC
 select_max = huge(1.d0)
END_PROVIDER

