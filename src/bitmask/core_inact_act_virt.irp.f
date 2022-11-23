use bitmasks

BEGIN_PROVIDER [ integer, n_core_orb]
  implicit none
  BEGIN_DOC
  ! Number of core MOs
  END_DOC
  integer                        :: i
  
  n_core_orb = 0
  do i = 1, mo_num
    if(mo_class(i) == 'Core')then
      n_core_orb += 1
    endif
  enddo
  
  call write_int(6,n_core_orb, 'Number of core     MOs')
   
END_PROVIDER

BEGIN_PROVIDER [ integer, n_inact_orb ]
  implicit none
  BEGIN_DOC
  ! Number of inactive MOs
  END_DOC
  integer                        :: i
  
  n_inact_orb = 0
  do i = 1, mo_num
    if (mo_class(i) == 'Inactive')then
      n_inact_orb += 1
    endif
  enddo
  
  call write_int(6,n_inact_orb,'Number of inactive MOs')
  
END_PROVIDER

BEGIN_PROVIDER [ integer, n_act_orb]
  implicit none
  BEGIN_DOC
  ! Number of active MOs
  END_DOC
  integer                        :: i
  
  n_act_orb = 0
  do i = 1, mo_num
    if (mo_class(i) == 'Active')then
      n_act_orb += 1
    endif
  enddo
  call write_int(6,n_act_orb, 'Number of active   MOs')
  if (mpi_master) then
    call ezfio_set_bitmask_n_act_orb(n_act_orb)
  endif
END_PROVIDER

BEGIN_PROVIDER [ integer, n_virt_orb ]
  implicit none
  BEGIN_DOC
  ! Number of virtual MOs
  END_DOC
  integer                        :: i
  
  n_virt_orb = 0
  do i = 1, mo_num
    if (mo_class(i) == 'Virtual')then
      n_virt_orb += 1
    endif
  enddo
  
  call write_int(6,n_virt_orb, 'Number of virtual  MOs')
   
END_PROVIDER

BEGIN_PROVIDER [ integer, n_del_orb ]
  implicit none
  BEGIN_DOC
  ! Number of deleted MOs
  END_DOC
  integer                        :: i
  
  n_del_orb = 0
  do i = 1, mo_num
    if (mo_class(i) == 'Deleted')then
      n_del_orb += 1
    endif
  enddo
  
  call write_int(6,n_del_orb, 'Number of deleted  MOs')
   
END_PROVIDER


BEGIN_PROVIDER [ integer, n_core_inact_orb ]
  implicit none
  BEGIN_DOC
  ! n_core + n_inact
  END_DOC
  integer                        :: i
  n_core_inact_orb = 0
  do i = 1, N_int
    n_core_inact_orb += popcnt(reunion_of_core_inact_bitmask(i,1))
  enddo
END_PROVIDER

BEGIN_PROVIDER [integer, n_inact_act_orb ]
   implicit none
  BEGIN_DOC
  ! n_inact + n_act
  END_DOC
   n_inact_act_orb = (n_inact_orb+n_act_orb)
END_PROVIDER
 
BEGIN_PROVIDER [integer, dim_list_core_orb]
  implicit none
  BEGIN_DOC
  ! dimensions for the allocation of list_core.
  ! it is at least 1
  END_DOC
   dim_list_core_orb = max(n_core_orb,1)
END_PROVIDER

BEGIN_PROVIDER [integer, dim_list_inact_orb]
   implicit none
   BEGIN_DOC
   ! dimensions for the allocation of list_inact.
   ! it is at least 1
   END_DOC
   dim_list_inact_orb = max(n_inact_orb,1)
END_PROVIDER

BEGIN_PROVIDER [integer, dim_list_core_inact_orb]
  implicit none
  BEGIN_DOC
  ! dimensions for the allocation of list_core.
  ! it is at least 1
  END_DOC
   dim_list_core_inact_orb = max(n_core_inact_orb,1)
END_PROVIDER

BEGIN_PROVIDER [integer, dim_list_act_orb]
   implicit none
   BEGIN_DOC
   ! dimensions for the allocation of list_act.
   ! it is at least 1
   END_DOC
   dim_list_act_orb = max(n_act_orb,1)
END_PROVIDER

BEGIN_PROVIDER [integer, dim_list_virt_orb]
   implicit none
   BEGIN_DOC
   ! dimensions for the allocation of list_virt.
   ! it is at least 1
   END_DOC
   dim_list_virt_orb = max(n_virt_orb,1)
END_PROVIDER

BEGIN_PROVIDER [integer, dim_list_del_orb]
   implicit none
   BEGIN_DOC
   ! dimensions for the allocation of list_del.
   ! it is at least 1
   END_DOC
   dim_list_del_orb = max(n_del_orb,1)
END_PROVIDER

BEGIN_PROVIDER [integer, n_core_inact_act_orb ]
  implicit none
  BEGIN_DOC
  !  Number of core inactive and active MOs
  END_DOC
  n_core_inact_act_orb = (n_core_orb + n_inact_orb + n_act_orb)
END_PROVIDER

 BEGIN_PROVIDER [ integer(bit_kind), core_bitmask , (N_int,2) ]
   implicit none
   BEGIN_DOC
   ! Bitmask identifying the core MOs 
   END_DOC
   core_bitmask  = 0_bit_kind
   if(n_core_orb > 0)then
     call list_to_bitstring( core_bitmask(1,1), list_core, n_core_orb, N_int)
     call list_to_bitstring( core_bitmask(1,2), list_core, n_core_orb, N_int)
   endif
 END_PROVIDER

 BEGIN_PROVIDER [ integer(bit_kind), inact_bitmask, (N_int,2) ]
   implicit none
   BEGIN_DOC
   ! Bitmask identifying the  inactive MOs 
   END_DOC
   inact_bitmask = 0_bit_kind
   if(n_inact_orb > 0)then
     call list_to_bitstring( inact_bitmask(1,1), list_inact, n_inact_orb, N_int)
     call list_to_bitstring( inact_bitmask(1,2), list_inact, n_inact_orb, N_int)
   endif
 END_PROVIDER

 BEGIN_PROVIDER [ integer(bit_kind), act_bitmask  , (N_int,2) ]
   implicit none
   BEGIN_DOC
   ! Bitmask identifying the active MOs 
   END_DOC
   act_bitmask   = 0_bit_kind
   if(n_act_orb > 0)then
     call list_to_bitstring( act_bitmask(1,1), list_act, n_act_orb, N_int)
     call list_to_bitstring( act_bitmask(1,2), list_act, n_act_orb, N_int)
   endif
  END_PROVIDER

 BEGIN_PROVIDER [ integer(bit_kind), virt_bitmask , (N_int,2) ]
   implicit none
   BEGIN_DOC
   ! Bitmask identifying the virtual MOs 
   END_DOC
   virt_bitmask  = 0_bit_kind
   if(n_virt_orb > 0)then
     call list_to_bitstring( virt_bitmask(1,1), list_virt, n_virt_orb, N_int)
     call list_to_bitstring( virt_bitmask(1,2), list_virt, n_virt_orb, N_int)
   endif
 END_PROVIDER

 BEGIN_PROVIDER [ integer(bit_kind), del_bitmask  , (N_int,2) ]
   implicit none
   BEGIN_DOC
   ! Bitmask identifying the deleted MOs 
   END_DOC

   del_bitmask   = 0_bit_kind
   
   if(n_del_orb > 0)then
     call list_to_bitstring( del_bitmask(1,1), list_del, n_del_orb, N_int)
     call list_to_bitstring( del_bitmask(1,2), list_del, n_del_orb, N_int)
   endif
   
 END_PROVIDER




 
 BEGIN_PROVIDER [ integer, list_core        , (dim_list_core_orb) ]
&BEGIN_PROVIDER [ integer, list_core_reverse, (mo_num) ]
   implicit none
   BEGIN_DOC
   ! List of MO indices which are in the core.
   END_DOC
   integer                        :: i, n
   list_core = 0
   list_core_reverse = 0

   n=0
   do i = 1, mo_num
     if(mo_class(i) == 'Core')then
       n += 1
       list_core(n) = i
       list_core_reverse(i) = n
     endif
   enddo
   print *,  'Core MOs:'
   print *,  list_core(1:n_core_orb)
   
END_PROVIDER
 
 BEGIN_PROVIDER [ integer, list_inact        , (dim_list_inact_orb) ]
&BEGIN_PROVIDER [ integer, list_inact_reverse, (mo_num) ]
   implicit none
   BEGIN_DOC
   ! List of MO indices which are inactive.
   END_DOC
   integer                        :: i, n
   list_inact = 0
   list_inact_reverse = 0

   n=0
   do i = 1, mo_num
     if (mo_class(i) == 'Inactive')then
       n += 1
       list_inact(n) = i
       list_inact_reverse(i) = n
     endif
   enddo
   print *,  'Inactive MOs:'
   print *,  list_inact(1:n_inact_orb)
   
END_PROVIDER
 
 BEGIN_PROVIDER [ integer, list_virt        , (dim_list_virt_orb) ]
&BEGIN_PROVIDER [ integer, list_virt_reverse, (mo_num) ]
   implicit none
   BEGIN_DOC
   ! List of MO indices which are virtual
   END_DOC
   integer                        :: i, n
   list_virt = 0
   list_virt_reverse = 0

   n=0
   do i = 1, mo_num
     if (mo_class(i) == 'Virtual')then
       n += 1
       list_virt(n) = i
       list_virt_reverse(i) = n
     endif
   enddo
   print *,  'Virtual MOs:'
   print *,  list_virt(1:n_virt_orb)
   
END_PROVIDER
 
 BEGIN_PROVIDER [ integer, list_del        , (dim_list_del_orb) ]
&BEGIN_PROVIDER [ integer, list_del_reverse, (mo_num) ]
   implicit none
   BEGIN_DOC
   ! List of MO indices which are deleted.
   END_DOC
   integer                        :: i, n
   list_del = 0
   list_del_reverse = 0

   n=0
   do i = 1, mo_num
     if (mo_class(i) == 'Deleted')then
       n += 1
       list_del(n) = i
       list_del_reverse(i) = n
     endif
   enddo
   print *,  'Deleted MOs:'
   print *,  list_del(1:n_del_orb)
   
END_PROVIDER
 
 BEGIN_PROVIDER [ integer, list_act        , (dim_list_act_orb) ]
&BEGIN_PROVIDER [ integer, list_act_reverse, (mo_num) ]
   implicit none
   BEGIN_DOC
   ! List of MO indices which are in the active.
   END_DOC
   integer                        :: i, n
   list_act = 0
   list_act_reverse = 0

   n=0
   do i = 1, mo_num
     if (mo_class(i) == 'Active')then
       n += 1
       list_act(n) = i
       list_act_reverse(i) = n
     endif
   enddo
   print *,  'Active MOs:'
   print *,  list_act(1:n_act_orb)
   
END_PROVIDER
 

 
 BEGIN_PROVIDER [ integer, list_core_inact        , (dim_list_core_inact_orb) ]
&BEGIN_PROVIDER [ integer, list_core_inact_reverse, (mo_num) ]
   implicit none
   BEGIN_DOC
   ! List of indices of the core and inactive MOs
   END_DOC
   integer                        :: i,itmp
   call bitstring_to_list(reunion_of_core_inact_bitmask(1,1), list_core_inact, itmp, N_int)
   list_core_inact_reverse = 0
   ASSERT (itmp == n_core_inact_orb)
   do i = 1, n_core_inact_orb
     list_core_inact_reverse(list_core_inact(i)) = i
   enddo
   print *,  'Core and Inactive MOs:'
   print *,  list_core_inact(1:n_core_inact_orb)
END_PROVIDER
 
 
 BEGIN_PROVIDER [ integer, list_core_inact_act        , (n_core_inact_act_orb) ]
&BEGIN_PROVIDER [ integer, list_core_inact_act_reverse, (mo_num) ]
   implicit none
   BEGIN_DOC
   ! List of indices of the core inactive and active MOs
   END_DOC
   integer                        :: i,itmp
   call bitstring_to_list(reunion_of_core_inact_act_bitmask(1,1), list_core_inact_act, itmp, N_int)
   list_core_inact_act_reverse = 0
   ASSERT (itmp == n_core_inact_act_orb)
   do i = 1, n_core_inact_act_orb
     list_core_inact_act_reverse(list_core_inact_act(i)) = i
   enddo
   print *,  'Core, Inactive and Active MOs:'
   print *,  list_core_inact_act(1:n_core_inact_act_orb)
END_PROVIDER
 
 
 BEGIN_PROVIDER [ integer, list_inact_act        , (n_inact_act_orb) ]
&BEGIN_PROVIDER [ integer, list_inact_act_reverse, (mo_num) ]
   implicit none
   BEGIN_DOC
   ! List of indices of the inactive and active MOs
   END_DOC
   integer                        :: i,itmp
   call bitstring_to_list(reunion_of_inact_act_bitmask(1,1), list_inact_act, itmp, N_int)
   list_inact_act_reverse = 0
   ASSERT (itmp == n_inact_act_orb)
   do i = 1, n_inact_act_orb
     list_inact_act_reverse(list_inact_act(i)) = i
   enddo
   print *,  'Inactive and Active MOs:'
   print *,  list_inact_act(1:n_inact_act_orb)
END_PROVIDER
 

BEGIN_PROVIDER [integer, n_all_but_del_orb]
 implicit none
 integer :: i
 n_all_but_del_orb = 0
 do i = 1, mo_num
  if(  trim(mo_class(i))=="Core" & 
  .or. trim(mo_class(i))=="Inactive" & 
  .or. trim(mo_class(i))=="Active"   &
  .or. trim(mo_class(i))=="Virtual" )then
   n_all_but_del_orb +=1 
  endif
 enddo
END_PROVIDER 

BEGIN_PROVIDER [integer, list_all_but_del_orb, (n_all_but_del_orb)]
 implicit none
 integer :: i,j
 j = 0 
 do i = 1, mo_num
  if(  trim(mo_class(i))=="Core" & 
  .or. trim(mo_class(i))=="Inactive" & 
  .or. trim(mo_class(i))=="Active"   &
  .or. trim(mo_class(i))=="Virtual" )then
   j += 1
   list_all_but_del_orb(j) = i
  endif
 enddo

END_PROVIDER
