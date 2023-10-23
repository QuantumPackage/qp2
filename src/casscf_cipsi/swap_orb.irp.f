 BEGIN_PROVIDER [double precision, SXvector_lowest, (nMonoEx)]
 implicit none
 integer :: i 
  do i=2,nMonoEx+1
   SXvector_lowest(i-1)=SXeigenvec(i,1)
  enddo
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, thresh_overlap_switch]
 implicit none
 thresh_overlap_switch = 0.5d0
 END_PROVIDER 

 BEGIN_PROVIDER [integer, max_overlap, (nMonoEx)]
&BEGIN_PROVIDER [integer, n_max_overlap]
&BEGIN_PROVIDER [integer, dim_n_max_overlap]
 implicit none
 double precision, allocatable :: vec_tmp(:)
 integer, allocatable :: iorder(:)
 allocate(vec_tmp(nMonoEx),iorder(nMonoEx))
 integer :: i
 do i = 1, nMonoEx
  iorder(i)  = i
  vec_tmp(i) = -dabs(SXvector_lowest(i))
 enddo
 call dsort(vec_tmp,iorder,nMonoEx)
 n_max_overlap = 0
 do i = 1, nMonoEx
  if(dabs(vec_tmp(i)).gt.thresh_overlap_switch)then
   n_max_overlap += 1
   max_overlap(n_max_overlap) = iorder(i)
  endif
 enddo
 dim_n_max_overlap = max(1,n_max_overlap)
 END_PROVIDER 

 BEGIN_PROVIDER [integer, orb_swap, (2,dim_n_max_overlap)]
&BEGIN_PROVIDER [integer, index_orb_swap, (dim_n_max_overlap)]
&BEGIN_PROVIDER [integer, n_orb_swap ]
 implicit none
 use bitmasks ! you need to include the bitmasks_module.f90 features
 integer :: i,imono,iorb,jorb,j
 n_orb_swap = 0
 do i = 1, n_max_overlap
  imono = max_overlap(i)
  iorb = excit(1,imono)
  jorb = excit(2,imono)
  if (excit_class(imono) == "c-a" .and.hessmat(imono,imono).gt.0.d0)then ! core --> active rotation
   n_orb_swap += 1
   orb_swap(1,n_orb_swap) = iorb ! core 
   orb_swap(2,n_orb_swap) = jorb ! active
   index_orb_swap(n_orb_swap) = imono
  else if (excit_class(imono) == "a-v" .and.hessmat(imono,imono).gt.0.d0)then ! active --> virtual rotation
   n_orb_swap += 1
   orb_swap(1,n_orb_swap) = jorb ! virtual 
   orb_swap(2,n_orb_swap) = iorb ! active
   index_orb_swap(n_orb_swap) = imono
  endif
 enddo

 integer,allocatable :: orb_swap_tmp(:,:)
 allocate(orb_swap_tmp(2,dim_n_max_overlap))
 do i = 1, n_orb_swap
  orb_swap_tmp(1,i) = orb_swap(1,i)
  orb_swap_tmp(2,i) = orb_swap(2,i)
 enddo
 
 integer(bit_kind), allocatable :: det_i(:),det_j(:)
 allocate(det_i(N_int),det_j(N_int))
 logical, allocatable :: good_orb_rot(:)
 allocate(good_orb_rot(n_orb_swap))
 integer, allocatable ::  index_orb_swap_tmp(:) 
 allocate(index_orb_swap_tmp(dim_n_max_overlap))
 index_orb_swap_tmp = index_orb_swap
 good_orb_rot = .True.
 integer :: icount,k
 do i = 1, n_orb_swap
  if(.not.good_orb_rot(i))cycle
  det_i = 0_bit_kind 
  call set_bit_to_integer(orb_swap(1,i),det_i,N_int)
  call set_bit_to_integer(orb_swap(2,i),det_i,N_int)
  do j = i+1, n_orb_swap
   det_j = 0_bit_kind 
   call set_bit_to_integer(orb_swap(1,j),det_j,N_int)
   call set_bit_to_integer(orb_swap(2,j),det_j,N_int)
   icount = 0
   do k = 1, N_int
    icount += popcnt(ior(det_i(k),det_j(k)))
   enddo
   if (icount.ne.4)then
    good_orb_rot(i) = .False.
    good_orb_rot(j) = .False.
    exit
   endif
  enddo
 enddo
 icount = n_orb_swap
 n_orb_swap = 0
 do i = 1, icount
  if(good_orb_rot(i))then
   n_orb_swap += 1
   index_orb_swap(n_orb_swap) = index_orb_swap_tmp(i)
   orb_swap(1,n_orb_swap) = orb_swap_tmp(1,i)
   orb_swap(2,n_orb_swap) = orb_swap_tmp(2,i)
  endif
 enddo

 if(n_orb_swap.gt.0)then
  print*,'n_orb_swap = ',n_orb_swap
 endif
 do i = 1, n_orb_swap
  print*,'imono = ',index_orb_swap(i)
  print*,orb_swap(1,i),'-->',orb_swap(2,i)
 enddo
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, switch_mo_coef, (ao_num,mo_num)]
  implicit none
  integer :: i,j,iorb,jorb
  switch_mo_coef = NatOrbsFCI
  do i = 1, n_orb_swap
   iorb = orb_swap(1,i)
   jorb = orb_swap(2,i)
   do j = 1, ao_num
    switch_mo_coef(j,jorb) = NatOrbsFCI(j,iorb)
   enddo
   do j = 1, ao_num
    switch_mo_coef(j,iorb) = NatOrbsFCI(j,jorb)
   enddo
  enddo

 END_PROVIDER 
