use bitmasks


  BEGIN_PROVIDER [ integer, n_core_orb]
 &BEGIN_PROVIDER [ integer, n_inact_orb ]
 &BEGIN_PROVIDER [ integer, n_act_orb]
 &BEGIN_PROVIDER [ integer, n_virt_orb ]
 &BEGIN_PROVIDER [ integer, n_del_orb ]
 implicit none
 BEGIN_DOC
 ! inact_bitmask : Bitmask of the inactive orbitals which are supposed to be doubly excited 
 ! in post CAS methods
 ! n_inact_orb   : Number of inactive orbitals
 ! virt_bitmask  : Bitmaks of vritual orbitals which are supposed to be recieve electrons 
 ! in post CAS methods
 ! n_virt_orb    : Number of virtual orbitals
 ! list_inact : List of the inactive orbitals which are supposed to be doubly excited 
 ! in post CAS methods
 ! list_virt  : List of vritual orbitals which are supposed to be recieve electrons 
 ! in post CAS methods
 ! list_inact_reverse : reverse list of inactive orbitals 
 ! list_inact_reverse(i) = 0 ::> not an inactive 
 ! list_inact_reverse(i) = k ::> IS the kth inactive 
 ! list_virt_reverse : reverse list of virtual orbitals 
 ! list_virt_reverse(i) = 0 ::> not an virtual 
 ! list_virt_reverse(i) = k ::> IS the kth virtual 
 ! list_act(i) = index of the ith active orbital
 ! 
 ! list_act_reverse : reverse list of active orbitals 
 ! list_act_reverse(i) = 0 ::> not an active 
 ! list_act_reverse(i) = k ::> IS the kth active orbital
 END_DOC
 logical                        :: exists
 integer                        :: j,i

 n_core_orb = 0
 n_inact_orb = 0
 n_act_orb = 0
 n_virt_orb = 0
 n_del_orb = 0
 do i = 1, mo_num
  if(mo_class(i) == 'Core')then
   n_core_orb += 1
  else if (mo_class(i) == 'Inactive')then
   n_inact_orb += 1
  else if (mo_class(i) == 'Active')then
   n_act_orb += 1
  else if (mo_class(i) == 'Virtual')then
   n_virt_orb += 1
  else if (mo_class(i) == 'Deleted')then
   n_del_orb += 1
  endif
 enddo


 call write_int(6,n_core_orb, 'Number of core     MOs')
 call write_int(6,n_inact_orb,'Number of inactive MOs')
 call write_int(6,n_act_orb, 'Number of active   MOs')
 call write_int(6,n_virt_orb, 'Number of virtual  MOs')
 call write_int(6,n_del_orb, 'Number of deleted  MOs')

 END_PROVIDER


 BEGIN_PROVIDER [integer, dim_list_core_orb]
&BEGIN_PROVIDER [integer, dim_list_inact_orb]
&BEGIN_PROVIDER [integer, dim_list_virt_orb]
&BEGIN_PROVIDER [integer, dim_list_act_orb]
&BEGIN_PROVIDER [integer, dim_list_del_orb]
 implicit none
 BEGIN_DOC
! dimensions for the allocation of list_inact, list_virt, list_core and list_act
! it is at least 1
 END_DOC
 dim_list_core_orb = max(n_core_orb,1)
 dim_list_inact_orb = max(n_inact_orb,1)
 dim_list_virt_orb = max(n_virt_orb,1)
 dim_list_act_orb = max(n_act_orb,1)
 dim_list_del_orb = max(n_del_orb,1)
END_PROVIDER 

 BEGIN_PROVIDER [ integer, list_inact, (dim_list_inact_orb)]
&BEGIN_PROVIDER [ integer, list_virt, (dim_list_virt_orb)]
&BEGIN_PROVIDER [ integer, list_inact_reverse, (mo_num)]
&BEGIN_PROVIDER [ integer, list_virt_reverse, (mo_num)]
&BEGIN_PROVIDER [ integer, list_del_reverse, (mo_num)]
&BEGIN_PROVIDER [ integer, list_del, (mo_num)]
&BEGIN_PROVIDER [integer, list_core, (dim_list_core_orb)]
&BEGIN_PROVIDER [integer, list_core_reverse, (mo_num)]
&BEGIN_PROVIDER [integer, list_act, (dim_list_act_orb)]
&BEGIN_PROVIDER [integer, list_act_reverse, (mo_num)]
&BEGIN_PROVIDER [ integer(bit_kind), core_bitmask,  (N_int,2)]
&BEGIN_PROVIDER [ integer(bit_kind), inact_bitmask, (N_int,2) ]
&BEGIN_PROVIDER [ integer(bit_kind), act_bitmask,   (N_int,2) ]
&BEGIN_PROVIDER [ integer(bit_kind), virt_bitmask,  (N_int,2) ]
&BEGIN_PROVIDER [ integer(bit_kind), del_bitmask,  (N_int,2) ]
 implicit none
 BEGIN_DOC
 ! inact_bitmask : Bitmask of the inactive orbitals which are supposed to be doubly excited 
 ! in post CAS methods
 ! n_inact_orb   : Number of inactive orbitals
 ! virt_bitmask  : Bitmaks of vritual orbitals which are supposed to be recieve electrons 
 ! in post CAS methods
 ! n_virt_orb    : Number of virtual orbitals
 ! list_inact : List of the inactive orbitals which are supposed to be doubly excited 
 ! in post CAS methods
 ! list_virt  : List of vritual orbitals which are supposed to be recieve electrons 
 ! in post CAS methods
 ! list_inact_reverse : reverse list of inactive orbitals 
 ! list_inact_reverse(i) = 0 ::> not an inactive 
 ! list_inact_reverse(i) = k ::> IS the kth inactive 
 ! list_virt_reverse : reverse list of virtual orbitals 
 ! list_virt_reverse(i) = 0 ::> not an virtual 
 ! list_virt_reverse(i) = k ::> IS the kth virtual 
 ! list_act(i) = index of the ith active orbital
 ! 
 ! list_act_reverse : reverse list of active orbitals 
 ! list_act_reverse(i) = 0 ::> not an active 
 ! list_act_reverse(i) = k ::> IS the kth active orbital
 END_DOC
 logical                        :: exists
 integer                        :: j,i
 integer                        :: n_core_orb_tmp, n_inact_orb_tmp, n_act_orb_tmp, n_virt_orb_tmp,n_del_orb_tmp
 integer                        :: list_core_tmp(N_int*bit_kind_size)
 integer                        :: list_inact_tmp(N_int*bit_kind_size)
 integer                        :: list_act_tmp(N_int*bit_kind_size)
 integer                        :: list_virt_tmp(N_int*bit_kind_size)
 integer                        :: list_del_tmp(N_int*bit_kind_size)
 list_core = 0
 list_inact = 0
 list_act = 0
 list_virt = 0
 list_del = 0
 list_core_reverse = 0
 list_inact_reverse = 0
 list_act_reverse = 0
 list_virt_reverse = 0
 list_del_reverse = 0
 n_core_orb_tmp = 0
 n_inact_orb_tmp = 0
 n_act_orb_tmp = 0
 n_virt_orb_tmp = 0
 n_del_orb_tmp = 0
 do i = 1, mo_num
  if(mo_class(i) == 'Core')then
   n_core_orb_tmp += 1
   list_core(n_core_orb_tmp) = i
   list_core_tmp(n_core_orb_tmp) = i
   list_core_reverse(i) = n_core_orb_tmp
  else if (mo_class(i) == 'Inactive')then
   n_inact_orb_tmp += 1
   list_inact(n_inact_orb_tmp) = i
   list_inact_tmp(n_inact_orb_tmp) = i
   list_inact_reverse(i) = n_inact_orb_tmp
  else if (mo_class(i) == 'Active')then
   n_act_orb_tmp += 1
   list_act(n_act_orb_tmp) = i
   list_act_tmp(n_act_orb_tmp) = i
   list_act_reverse(i) = n_act_orb_tmp
  else if (mo_class(i) == 'Virtual')then
   n_virt_orb_tmp += 1
   list_virt(n_virt_orb_tmp) = i
   list_virt_tmp(n_virt_orb_tmp) = i
   list_virt_reverse(i) = n_virt_orb_tmp
  else if (mo_class(i) == 'Deleted')then
   n_del_orb_tmp += 1
   list_del(n_del_orb_tmp) = i
   list_del_tmp(n_del_orb_tmp) = i
   list_del_reverse(i) = n_del_orb_tmp
  endif
 enddo

 if(n_core_orb.ne.0)then
  call list_to_bitstring( core_bitmask(1,1), list_core, n_core_orb, N_int)
  call list_to_bitstring( core_bitmask(1,2), list_core, n_core_orb, N_int)
 endif
 if(n_inact_orb.ne.0)then
  call list_to_bitstring( inact_bitmask(1,1), list_inact, n_inact_orb, N_int)
  call list_to_bitstring( inact_bitmask(1,2), list_inact, n_inact_orb, N_int)
 endif
 if(n_act_orb.ne.0)then
  call list_to_bitstring( act_bitmask(1,1), list_act, n_act_orb, N_int)
  call list_to_bitstring( act_bitmask(1,2), list_act, n_act_orb, N_int)
 endif
 if(n_virt_orb.ne.0)then
  call list_to_bitstring( virt_bitmask(1,1), list_virt, n_virt_orb, N_int)
  call list_to_bitstring( virt_bitmask(1,2), list_virt, n_virt_orb, N_int)
 endif
 if(n_del_orb.ne.0)then
  call list_to_bitstring( del_bitmask(1,1), list_del, n_del_orb, N_int)
  call list_to_bitstring( del_bitmask(1,2), list_del, n_del_orb, N_int)
 endif


END_PROVIDER 

BEGIN_PROVIDER [integer, n_inact_act ]
 implicit none
 n_inact_act = (n_inact_orb+n_act_orb)

END_PROVIDER 

BEGIN_PROVIDER [integer, list_inact_act, (n_inact_act)]
 integer :: i,itmp
 itmp = 0
 do i = 1, n_inact_orb
  itmp += 1
  list_inact_act(itmp) = list_inact(i) 
 enddo
 do i = 1, n_act_orb
  itmp += 1
  list_inact_act(itmp) = list_act(i) 
 enddo
END_PROVIDER 

BEGIN_PROVIDER [integer, n_core_inact_act_orb ]
 implicit none
 n_core_inact_act_orb = (n_core_orb + n_inact_orb + n_act_orb)

END_PROVIDER 

 BEGIN_PROVIDER [integer, list_core_inact_act, (n_core_inact_act_orb)]
&BEGIN_PROVIDER [ integer, list_core_inact_act_reverse, (n_core_inact_act_orb)]
 integer :: i,itmp
 itmp = 0
 do i = 1, n_core_orb
  itmp += 1
  list_core_inact_act(itmp) = list_core(i) 
 enddo
 do i = 1, n_inact_orb
  itmp += 1
  list_core_inact_act(itmp) = list_inact(i) 
 enddo
 do i = 1, n_act_orb
  itmp += 1
  list_core_inact_act(itmp) = list_act(i) 
 enddo

 integer :: occ_inact(N_int*bit_kind_size)
 occ_inact = 0
 call bitstring_to_list(reunion_of_core_inact_act_bitmask(1,1), occ_inact(1), itest, N_int)
 list_inact_reverse = 0
 do i = 1, n_core_inact_act_orb
  list_core_inact_act_reverse(occ_inact(i)) = i
 enddo
END_PROVIDER 
