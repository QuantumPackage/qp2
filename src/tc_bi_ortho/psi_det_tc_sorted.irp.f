use bitmasks

BEGIN_PROVIDER [ double precision, psi_average_norm_contrib_tc, (psi_det_size) ]                                                                
  implicit none
  BEGIN_DOC
  ! Contribution of determinants to the state-averaged density.
  END_DOC
  integer                        :: i,j,k
  double precision               :: f
    
  psi_average_norm_contrib_tc(:) = 0.d0
  do k=1,N_states
    do i=1,N_det
      psi_average_norm_contrib_tc(i) = psi_average_norm_contrib_tc(i) +    &
          dabs(psi_l_coef_bi_ortho(i,k)*psi_r_coef_bi_ortho(i,k))*state_average_weight(k)
    enddo
  enddo
  f = 1.d0/sum(psi_average_norm_contrib_tc(1:N_det))
  do i=1,N_det
    psi_average_norm_contrib_tc(i) = psi_average_norm_contrib_tc(i)*f
  enddo
END_PROVIDER


 BEGIN_PROVIDER [ integer(bit_kind), psi_det_sorted_tc, (N_int,2,psi_det_size) ]
&BEGIN_PROVIDER [ double precision, psi_coef_sorted_tc, (psi_det_size,N_states) ]
&BEGIN_PROVIDER [ double precision, psi_average_norm_contrib_sorted_tc, (psi_det_size) ]
&BEGIN_PROVIDER [ integer, psi_det_sorted_tc_order, (psi_det_size) ]
   implicit none
   BEGIN_DOC
   ! Wave function sorted by determinants contribution to the norm (state-averaged)
   !
   ! psi_det_sorted_tc_order(i) -> k : index in psi_det
   END_DOC
   integer                        :: i,j,k
   integer, allocatable           :: iorder(:)
   print *, 'providing psi_det_sorted_tc'
   allocate ( iorder(N_det) )
   print*,'before '
   do i=1,N_det
     psi_average_norm_contrib_sorted_tc(i) = -psi_average_norm_contrib_tc(i)
     iorder(i) = i
     print*,i,iorder(i),psi_average_norm_contrib_sorted_tc(i)
   enddo
   call dsort(psi_average_norm_contrib_sorted_tc,iorder,N_det)
   print*,'after '
   do i=1,N_det
!     iorder(i) = i
     print*,i,iorder(i),psi_average_norm_contrib_sorted_tc(i)
     do j=1,N_int
       psi_det_sorted_tc(j,1,i) = psi_det(j,1,iorder(i))
       psi_det_sorted_tc(j,2,i) = psi_det(j,2,iorder(i))
     enddo
    psi_average_norm_contrib_sorted_tc(i) = -psi_average_norm_contrib_sorted_tc(i)
    psi_det_sorted_tc_order(iorder(i)) = i
   enddo
   double precision :: accu
   do k=1,N_states
    accu = 0.d0
    do i=1,N_det
     psi_coef_sorted_tc(i,k) = dsqrt(dabs(psi_l_coef_bi_ortho(iorder(i),k)*psi_r_coef_bi_ortho(iorder(i),k)))
     accu += psi_coef_sorted_tc(i,k)**2
    enddo
    accu = 1.d0/dsqrt(accu)
    do i=1,N_det
     psi_coef_sorted_tc(i,k) *= accu
    enddo
   enddo

   psi_det_sorted_tc(:,:,N_det+1:psi_det_size) = 0_bit_kind
   psi_coef_sorted_tc(N_det+1:psi_det_size,:) = 0.d0
   psi_average_norm_contrib_sorted_tc(N_det+1:psi_det_size) = 0.d0
   psi_det_sorted_tc_order(N_det+1:psi_det_size) = 0

   deallocate(iorder)
  logical :: pouet
 pouet  = .true.
 do i = 1, N_det
  if(psi_average_norm_contrib_sorted_tc(i) == 0.d0)then
   pouet = .False.
   exit
  endif
 enddo
 
 if(pouet.and.N_det.ne.1)then
 print*,'writing sorted'
 do i = 1, N_det
!   call debug_det(psi_det_sorted_tc(1,1,i),N_int)
   print*,i,psi_average_norm_contrib_sorted_tc(i)
   write(35,*)psi_det_sorted_tc(1,1,i),psi_det_sorted_tc(1,2,i)
 enddo
 endif

END_PROVIDER

 BEGIN_PROVIDER [double precision, psi_r_coef_sorted_bi_ortho, (psi_det_size, N_states)]
&BEGIN_PROVIDER [double precision, psi_l_coef_sorted_bi_ortho, (psi_det_size, N_states)]
 BEGIN_DOC
 ! psi_r_coef_sorted_bi_ortho : right coefficients corresponding to psi_det_sorted_tc
 ! psi_l_coef_sorted_bi_ortho : left  coefficients corresponding to psi_det_sorted_tc
 END_DOC
   implicit none
   integer                       :: i, j, k
   psi_r_coef_sorted_bi_ortho = 0.d0
   psi_l_coef_sorted_bi_ortho = 0.d0
   do i = 1, N_det
    psi_r_coef_sorted_bi_ortho(i,1) = psi_r_coef_bi_ortho(psi_det_sorted_tc_order(i),1)
    psi_l_coef_sorted_bi_ortho(i,1) = psi_l_coef_bi_ortho(psi_det_sorted_tc_order(i),1)
   enddo
   logical :: pouet
   pouet  = .true.
   do i = 1, N_det
    if(psi_l_coef_sorted_bi_ortho(i,1) == 0.d0)then
     pouet = .False.
     exit
    endif
   enddo
   if(pouet.and.N_det.ne.1)then
   print*,'psi_r_coef_sorted_bi_ortho'
    do i = 1, N_det
     print*,psi_r_coef_bi_ortho(psi_det_sorted_tc_order(i),1)
     write(42,'(F10.7)')dabs(psi_r_coef_sorted_bi_ortho(i,1)*psi_l_coef_sorted_bi_ortho(i,1))
    enddo
   endif

END_PROVIDER


 BEGIN_PROVIDER [ integer(bit_kind), psi_det_sorted_tc_bit, (N_int,2,psi_det_size) ]
&BEGIN_PROVIDER [ double precision, psi_coef_sorted_tc_bit, (psi_det_size,N_states) ]
   implicit none
   BEGIN_DOC
   ! Determinants on which we apply $\langle i|H|psi \rangle$ for perturbation.
   ! They are sorted by determinants interpreted as integers. Useful
   ! to accelerate the search of a random determinant in the wave
   ! function.
   END_DOC

   call sort_dets_by_det_search_key(N_det, psi_det, psi_coef, size(psi_coef,1),       &
       psi_det_sorted_tc_bit, psi_coef_sorted_tc_bit, N_states)

END_PROVIDER


 BEGIN_PROVIDER [ integer(bit_kind), psi_det_sorted_tc_right, (N_int,2,N_det) ]
&BEGIN_PROVIDER [double precision, psi_r_coef_sorted_bi_ortho_right, (N_det)]
 implicit none
 BEGIN_DOC
 ! psi_det_sorted_tc_right : Slater determinants sorted by decreasing value of |right- coefficients|
 ! 
 ! psi_r_coef_sorted_bi_ortho_right : right wave function according to psi_det_sorted_tc_right
 END_DOC
 integer, allocatable           :: iorder(:)
 double precision, allocatable :: coef(:)
 integer :: i,j
 allocate ( iorder(N_det) , coef(N_det))
 do i=1,N_det
   coef(i) = -dabs(psi_r_coef_bi_ortho(i,1)/psi_r_coef_bi_ortho(1,1))
   iorder(i) = i
 enddo
 call dsort(coef,iorder,N_det)
 do i=1,N_det
   do j=1,N_int
     psi_det_sorted_tc_right(j,1,i) = psi_det(j,1,iorder(i))
     psi_det_sorted_tc_right(j,2,i) = psi_det(j,2,iorder(i))
   enddo
  psi_r_coef_sorted_bi_ortho_right(i) = psi_r_coef_bi_ortho(iorder(i),1)/psi_r_coef_bi_ortho(iorder(1),1)
 enddo
END_PROVIDER 

 BEGIN_PROVIDER [ integer(bit_kind), psi_det_sorted_tc_left, (N_int,2,N_det) ]
&BEGIN_PROVIDER [double precision, psi_l_coef_sorted_bi_ortho_left, (N_det)]
 implicit none
 BEGIN_DOC
 ! psi_det_sorted_tc_left : Slater determinants sorted by decreasing value of |LEFTt- coefficients|
 ! 
 ! psi_r_coef_sorted_bi_ortho_left : LEFT wave function according to psi_det_sorted_tc_left
 END_DOC
 integer, allocatable           :: iorder(:)
 double precision, allocatable :: coef(:)
 integer :: i,j
 allocate ( iorder(N_det) , coef(N_det))
 do i=1,N_det
   coef(i) = -dabs(psi_l_coef_bi_ortho(i,1)/psi_r_coef_bi_ortho(1,1))
   iorder(i) = i
 enddo
 call dsort(coef,iorder,N_det)
 do i=1,N_det
   do j=1,N_int
     psi_det_sorted_tc_left(j,1,i) = psi_det(j,1,iorder(i))
     psi_det_sorted_tc_left(j,2,i) = psi_det(j,2,iorder(i))
   enddo
  psi_l_coef_sorted_bi_ortho_left(i) = psi_l_coef_bi_ortho(iorder(i),1)/psi_l_coef_bi_ortho(iorder(1),1)
 enddo
END_PROVIDER 
