
use bitmasks

BEGIN_PROVIDER [ integer, N_det_generators_HF ]
 implicit none
 BEGIN_DOC
 ! For Single reference wave functions, the number of generators is 1 : the
 ! Hartree-Fock determinant
 END_DOC
 N_det_generators_HF = 1
END_PROVIDER

 BEGIN_PROVIDER [ integer(bit_kind), psi_det_generators_HF, (N_int,2,psi_det_size) ]
&BEGIN_PROVIDER [ double precision, psi_coef_generators_HF, (psi_det_size,N_states) ]
 implicit none
 BEGIN_DOC
 ! For Single reference wave functions, the generator is the
 ! Hartree-Fock determinant
 END_DOC
 psi_det_generators_HF = 0_bit_kind
 integer :: i,j
 integer :: degree

 do i=1,N_int
   psi_det_generators_HF(i,1,1) = HF_bitmask(i,1)
   psi_det_generators_HF(i,2,1) = HF_bitmask(i,2)
 enddo

 do j=1,N_det
   call get_excitation_degree(HF_bitmask,psi_det(1,1,j),degree,N_int)
   if (degree == 0) then
     exit
   endif
 end do

 psi_det_generators_HF(:,:,1) = psi_det(:,:,j)
 psi_coef_generators_HF(1,:) = 1.d0

END_PROVIDER

 BEGIN_PROVIDER [ integer          , HF_index ]
 implicit none
 integer :: j,degree
 do j=1,N_det
   call get_excitation_degree(HF_bitmask,psi_det_sorted(1,1,j),degree,N_int)
   if (degree == 0) then
     HF_index = j 
     exit
   endif
 end do
END_PROVIDER 
