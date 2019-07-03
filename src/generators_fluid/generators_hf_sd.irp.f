
use bitmasks

BEGIN_PROVIDER [ integer, N_det_generators_HF_SD ]
 implicit none
 BEGIN_DOC
 ! For Single reference wave functions, the number of generators is 1 : the
 ! Hartree-Fock determinant
 END_DOC
 N_det_generators_HF_SD = 0
 integer :: i,degree
 double precision :: thr
 double precision :: accu
 accu = 0.d0
 thr = threshold_generators
 do i = 1, N_det
  call get_excitation_degree(HF_bitmask,psi_det_sorted(1,1,i),degree,N_int)
  if(degree.le.2.and. accu .le. thr )then
   accu += psi_coef_sorted(i,1)**2
   N_det_generators_HF_SD += 1
  endif
 enddo
!print*,''
!print*,'N_det_generators_HF_SD = ',N_det_generators_HF_SD
END_PROVIDER

 BEGIN_PROVIDER [ integer(bit_kind), psi_det_generators_HF_SD, (N_int,2,psi_det_size) ]
&BEGIN_PROVIDER [ double precision, psi_coef_generators_HF_SD, (psi_det_size,N_states) ]
&BEGIN_PROVIDER [ integer(bit_kind), psi_det_sorted_gen_HF_SD, (N_int,2,psi_det_size) ]
&BEGIN_PROVIDER [ double precision, psi_coef_sorted_gen_HF_SD, (psi_det_size,N_states) ]
&BEGIN_PROVIDER [ integer, psi_det_sorted_gen_HF_SD_order, (psi_det_size) ]
 implicit none
 BEGIN_DOC
 ! For Single reference wave functions, the generator is the
 ! Hartree-Fock determinant
 END_DOC
 psi_det_generators_HF_SD = 0_bit_kind
 integer :: i,j,k
 integer :: degree
 double precision :: thr
 double precision :: accu
 integer, allocatable :: nongen(:)
 integer :: inongen

 allocate(nongen(N_det))

 thr = threshold_generators

 accu = 0.d0
 k = 0
 inongen = 0
 do j=1,N_det
   call get_excitation_degree(HF_bitmask,psi_det_sorted(1,1,j),degree,N_int)
   if(degree.le.2.and. accu.le.thr )then
    accu += psi_coef_sorted(j,1)**2
    k += 1
    psi_det_sorted_gen_HF_SD_order(j) = k
    do i = 1, N_int
     psi_det_generators_HF_SD(i,1,k) = psi_det_sorted(i,1,j)
     psi_det_generators_HF_SD(i,2,k) = psi_det_sorted(i,2,j)
    enddo
    do i = 1, N_states
     psi_coef_generators_HF_SD(k,i) = psi_coef_sorted(j,i)
    enddo
   else
    inongen += 1
    nongen(inongen) = j
   endif
 end do

 psi_det_sorted_gen_HF_SD(:,:,:N_det_generators_HF_SD) = psi_det_generators_HF_SD(:,:,:N_det_generators_HF_SD)
 psi_coef_sorted_gen_HF_SD(:N_det_generators_HF_SD, :) = psi_coef_generators_HF_SD(:N_det_generators_HF_SD, :)
 do i=1,inongen
   psi_det_sorted_gen_HF_SD_order(nongen(i)) = N_det_generators_HF_SD+i
   psi_det_sorted_gen_HF_SD(:,:,N_det_generators_HF_SD+i) = psi_det_sorted(:,:,nongen(i))
   psi_coef_sorted_gen_HF_SD(N_det_generators_HF_SD+i, :) = psi_coef_sorted(nongen(i),:)
 end do

END_PROVIDER

