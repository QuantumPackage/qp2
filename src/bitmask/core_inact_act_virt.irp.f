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
 
!============================================!
!                                            !
!                    kpts                    !
!                                            !
!============================================!
BEGIN_PROVIDER [ integer(bit_kind), kpts_bitmask , (N_int,kpt_num) ]
  implicit none
  BEGIN_DOC
  ! Bitmask identifying each kpt
  END_DOC
  integer :: k,i,di
  integer :: tmp_mo_list(mo_num_per_kpt)
  kpts_bitmask  = 0_bit_kind
  do k=1,kpt_num
    di=(k-1)*mo_num_per_kpt
    do i=1,mo_num_per_kpt
      tmp_mo_list(i) = i+di
    enddo
    call list_to_bitstring( kpts_bitmask(1,k), tmp_mo_list, mo_num_per_kpt, N_int)
  enddo
END_PROVIDER

BEGIN_PROVIDER [ integer, n_core_orb_kpts, (kpt_num)]
  implicit none
  BEGIN_DOC
  ! Number of core MOs
  END_DOC
  integer                        :: i,k,kshift

  do k=1,kpt_num
    n_core_orb_kpts(k) = 0
    kshift = (1-k)*mo_num_per_kpt
    do i = 1, mo_num_per_kpt
      if(mo_class(i+kshift) == 'Core')then
        n_core_orb_kpts(k) += 1
      endif
    enddo
  enddo
  
!  call write_int(6,n_core_orb, 'Number of core     MOs')
   
END_PROVIDER

BEGIN_PROVIDER [ integer, n_inact_orb_kpts, (kpt_num)]
  implicit none
  BEGIN_DOC
  ! Number of inactive MOs
  END_DOC
  integer                        :: i,k,kshift

  do k=1,kpt_num
    n_inact_orb_kpts(k) = 0
    kshift = (1-k)*mo_num_per_kpt
    do i = 1, mo_num_per_kpt
      if(mo_class(i+kshift) == 'Inactive')then
        n_inact_orb_kpts(k) += 1
      endif
    enddo
  enddo
  
!  call write_int(6,n_inact_orb, 'Number of inactive MOs')
   
END_PROVIDER

BEGIN_PROVIDER [ integer, n_act_orb_kpts, (kpt_num)]
  implicit none
  BEGIN_DOC
  ! Number of active MOs
  END_DOC
  integer                        :: i,k,kshift

  do k=1,kpt_num
    n_act_orb_kpts(k) = 0
    kshift = (1-k)*mo_num_per_kpt
    do i = 1, mo_num_per_kpt
      if(mo_class(i+kshift) == 'Active')then
        n_act_orb_kpts(k) += 1
      endif
    enddo
  enddo
  
!  call write_int(6,n_act_orb, 'Number of active   MOs')
   
END_PROVIDER

BEGIN_PROVIDER [ integer, n_virt_orb_kpts, (kpt_num)]
  implicit none
  BEGIN_DOC
  ! Number of virtual MOs
  END_DOC
  integer                        :: i,k,kshift

  do k=1,kpt_num
    n_virt_orb_kpts(k) = 0
    kshift = (1-k)*mo_num_per_kpt
    do i = 1, mo_num_per_kpt
      if(mo_class(i+kshift) == 'Virtual')then
        n_virt_orb_kpts(k) += 1
      endif
    enddo
  enddo
  
!  call write_int(6,n_virt_orb, 'Number of virtual  MOs')
   
END_PROVIDER

BEGIN_PROVIDER [ integer, n_del_orb_kpts, (kpt_num)]
  implicit none
  BEGIN_DOC
  ! Number of deleted MOs
  END_DOC
  integer                        :: i,k,kshift

  do k=1,kpt_num
    n_del_orb_kpts(k) = 0
    kshift = (1-k)*mo_num_per_kpt
    do i = 1, mo_num_per_kpt
      if(mo_class(i+kshift) == 'Deleted')then
        n_del_orb_kpts(k) += 1
      endif
    enddo
  enddo
  
!  call write_int(6,n_del_orb, 'Number of deleted  MOs')
   
END_PROVIDER

BEGIN_PROVIDER [ integer, n_core_inact_orb_kpts, (kpt_num) ]
  !todo: finish implementation for kpts (will need kpts_bitmask)
  implicit none
  BEGIN_DOC
  ! n_core + n_inact
  END_DOC
  integer                        :: i,k
  do k=1,kpt_num
  n_core_inact_orb_kpts(k) = 0
  do i = 1, N_int
    n_core_inact_orb_kpts(k) += popcnt(iand(kpts_bitmask(i,k),reunion_of_core_inact_bitmask(i,1)))
  enddo
  enddo
END_PROVIDER

BEGIN_PROVIDER [integer, n_inact_act_orb_kpts, (kpt_num) ]
   implicit none
  BEGIN_DOC
  ! n_inact + n_act
  END_DOC
  integer :: k
  do k=1,kpt_num
    n_inact_act_orb_kpts(k) = (n_inact_orb_kpts(k)+n_act_orb_kpts(k))
  enddo
END_PROVIDER
 
BEGIN_PROVIDER [integer, dim_list_core_orb_kpts]
  implicit none
  BEGIN_DOC
  ! dimensions for the allocation of list_core.
  ! it is at least 1
  END_DOC
   dim_list_core_orb_kpts = max(maxval(n_core_orb_kpts),1)
END_PROVIDER

BEGIN_PROVIDER [integer, dim_list_inact_orb_kpts]
   implicit none
   BEGIN_DOC
   ! dimensions for the allocation of list_inact.
   ! it is at least 1
   END_DOC
   dim_list_inact_orb_kpts = max(maxval(n_inact_orb_kpts),1)
END_PROVIDER

BEGIN_PROVIDER [integer, dim_list_core_inact_orb_kpts]
  implicit none
  BEGIN_DOC
  ! dimensions for the allocation of list_core.
  ! it is at least 1
  END_DOC
   dim_list_core_inact_orb_kpts = max(maxval(n_core_inact_orb_kpts),1)
END_PROVIDER

BEGIN_PROVIDER [integer, dim_list_act_orb_kpts]
   implicit none
   BEGIN_DOC
   ! dimensions for the allocation of list_act.
   ! it is at least 1
   END_DOC
   dim_list_act_orb_kpts = max(maxval(n_act_orb_kpts),1)
END_PROVIDER

BEGIN_PROVIDER [integer, dim_list_virt_orb_kpts]
   implicit none
   BEGIN_DOC
   ! dimensions for the allocation of list_virt.
   ! it is at least 1
   END_DOC
   dim_list_virt_orb_kpts = max(maxval(n_virt_orb_kpts),1)
END_PROVIDER

BEGIN_PROVIDER [integer, dim_list_del_orb_kpts]
   implicit none
   BEGIN_DOC
   ! dimensions for the allocation of list_del.
   ! it is at least 1
   END_DOC
   dim_list_del_orb_kpts = max(maxval(n_del_orb_kpts),1)
END_PROVIDER

BEGIN_PROVIDER [integer, dim_list_core_inact_act_orb_kpts]
   implicit none
   BEGIN_DOC
   ! dimensions for the allocation of list_core_inact_act.
   ! it is at least 1
   END_DOC
   dim_list_core_inact_act_orb_kpts = max(maxval(n_core_inact_act_orb_kpts),1)
END_PROVIDER

BEGIN_PROVIDER [integer, dim_list_inact_act_orb_kpts]
   implicit none
   BEGIN_DOC
   ! dimensions for the allocation of list_inact_act.
   ! it is at least 1
   END_DOC
   dim_list_inact_act_orb_kpts = max(maxval(n_inact_act_orb_kpts),1)
END_PROVIDER

BEGIN_PROVIDER [integer, n_core_inact_act_orb_kpts, (kpt_num) ]
  implicit none
  BEGIN_DOC
  !  Number of core inactive and active MOs
  END_DOC
  integer :: k
  do k=1,kpt_num
    n_core_inact_act_orb_kpts(k) = (n_core_orb_kpts(k) + n_inact_orb_kpts(k) + n_act_orb_kpts(k))
  enddo
END_PROVIDER
 



BEGIN_PROVIDER [ integer(bit_kind), core_bitmask_kpts , (N_int,2,kpt_num) ]
  implicit none
  BEGIN_DOC
  ! Bitmask identifying the core MOs 
  END_DOC
  integer :: k,i
  core_bitmask_kpts  = 0_bit_kind
  do k=1,kpt_num
    do i=1,N_int
      core_bitmask_kpts(i,1,k) = iand(core_bitmask(i,1),kpts_bitmask(i,k))
      core_bitmask_kpts(i,2,k) = iand(core_bitmask(i,2),kpts_bitmask(i,k))
    enddo
  enddo
END_PROVIDER

BEGIN_PROVIDER [ integer(bit_kind), inact_bitmask_kpts , (N_int,2,kpt_num) ]
  implicit none
  BEGIN_DOC
  ! Bitmask identifying the inactive MOs 
  END_DOC
  integer :: k,i
  inact_bitmask_kpts  = 0_bit_kind
  do k=1,kpt_num
    do i=1,N_int
      inact_bitmask_kpts(i,1,k) = iand(inact_bitmask(i,1),kpts_bitmask(i,k))
      inact_bitmask_kpts(i,2,k) = iand(inact_bitmask(i,2),kpts_bitmask(i,k))
    enddo
  enddo
END_PROVIDER

BEGIN_PROVIDER [ integer(bit_kind), act_bitmask_kpts , (N_int,2,kpt_num) ]
  implicit none
  BEGIN_DOC
  ! Bitmask identifying the active MOs 
  END_DOC
  integer :: k,i
  act_bitmask_kpts  = 0_bit_kind
  do k=1,kpt_num
    do i=1,N_int
      act_bitmask_kpts(i,1,k) = iand(act_bitmask(i,1),kpts_bitmask(i,k))
      act_bitmask_kpts(i,2,k) = iand(act_bitmask(i,2),kpts_bitmask(i,k))
    enddo
  enddo
END_PROVIDER

BEGIN_PROVIDER [ integer(bit_kind), virt_bitmask_kpts , (N_int,2,kpt_num) ]
  implicit none
  BEGIN_DOC
  ! Bitmask identifying the virtual MOs 
  END_DOC
  integer :: k,i
  virt_bitmask_kpts  = 0_bit_kind
  do k=1,kpt_num
    do i=1,N_int
      virt_bitmask_kpts(i,1,k) = iand(virt_bitmask(i,1),kpts_bitmask(i,k))
      virt_bitmask_kpts(i,2,k) = iand(virt_bitmask(i,2),kpts_bitmask(i,k))
    enddo
  enddo
END_PROVIDER

BEGIN_PROVIDER [ integer(bit_kind), del_bitmask_kpts , (N_int,2,kpt_num) ]
  implicit none
  BEGIN_DOC
  ! Bitmask identifying the deleted MOs 
  END_DOC
  integer :: k,i
  del_bitmask_kpts  = 0_bit_kind
  do k=1,kpt_num
    do i=1,N_int
      del_bitmask_kpts(i,1,k) = iand(del_bitmask(i,1),kpts_bitmask(i,k))
      del_bitmask_kpts(i,2,k) = iand(del_bitmask(i,2),kpts_bitmask(i,k))
    enddo
  enddo
END_PROVIDER

 BEGIN_PROVIDER [ integer, list_core_kpts        , (dim_list_core_orb_kpts,kpt_num) ]
&BEGIN_PROVIDER [ integer, list_core_kpts_reverse, (mo_num_per_kpt,kpt_num) ]
   implicit none
   BEGIN_DOC
   ! List of MO indices which are in the core.
   END_DOC
   integer                        :: i, n,k,di
   list_core_kpts = 0
   list_core_kpts_reverse = 0

   do k=1,kpt_num
     n=0
     di = (k-1)*mo_num_per_kpt
     do i = 1, mo_num_per_kpt
       if(mo_class(i+di) == 'Core')then
         n += 1
         list_core_kpts(n,k) = i
         list_core_kpts_reverse(i,k) = n
       endif
     enddo
     print *,  'Core MOs: ',k
     print *,  list_core_kpts(1:n_core_orb_kpts(k),k)
   enddo
   
END_PROVIDER
 
 BEGIN_PROVIDER [ integer, list_inact_kpts        , (dim_list_inact_orb_kpts,kpt_num) ]
&BEGIN_PROVIDER [ integer, list_inact_kpts_reverse, (mo_num_per_kpt,kpt_num) ]
   implicit none
   BEGIN_DOC
   ! List of MO indices which are inactive.
   END_DOC
   integer                        :: i, n,k,di
   list_inact_kpts = 0
   list_inact_kpts_reverse = 0

   do k=1,kpt_num
     n=0
     di = (k-1)*mo_num_per_kpt
     do i = 1, mo_num_per_kpt
       if(mo_class(i+di) == 'Inactive')then
         n += 1
         list_inact_kpts(n,k) = i
         list_inact_kpts_reverse(i,k) = n
       endif
     enddo
     print *,  'Inactive MOs: ',k
     print *,  list_inact_kpts(1:n_inact_orb_kpts(k),k)
   enddo
   
END_PROVIDER
 
 BEGIN_PROVIDER [ integer, list_virt_kpts        , (dim_list_virt_orb_kpts,kpt_num) ]
&BEGIN_PROVIDER [ integer, list_virt_kpts_reverse, (mo_num_per_kpt,kpt_num) ]
   implicit none
   BEGIN_DOC
   ! List of MO indices which are virtual.
   END_DOC
   integer                        :: i, n,k,di
   list_virt_kpts = 0
   list_virt_kpts_reverse = 0

   do k=1,kpt_num
     n=0
     di = (k-1)*mo_num_per_kpt
     do i = 1, mo_num_per_kpt
       if(mo_class(i+di) == 'Virtual')then
         n += 1
         list_virt_kpts(n,k) = i
         list_virt_kpts_reverse(i,k) = n
       endif
     enddo
     print *,  'Virtual MOs: ',k
     print *,  list_virt_kpts(1:n_virt_orb_kpts(k),k)
   enddo
   
END_PROVIDER
 
 BEGIN_PROVIDER [ integer, list_del_kpts        , (dim_list_del_orb_kpts,kpt_num) ]
&BEGIN_PROVIDER [ integer, list_del_kpts_reverse, (mo_num_per_kpt,kpt_num) ]
   implicit none
   BEGIN_DOC
   ! List of MO indices which are deleted.
   END_DOC
   integer                        :: i, n,k,di
   list_del_kpts = 0
   list_del_kpts_reverse = 0

   do k=1,kpt_num
     n=0
     di = (k-1)*mo_num_per_kpt
     do i = 1, mo_num_per_kpt
       if(mo_class(i+di) == 'Deleted')then
         n += 1
         list_del_kpts(n,k) = i
         list_del_kpts_reverse(i,k) = n
       endif
     enddo
     print *,  'Deleted MOs: ',k
     print *,  list_del_kpts(1:n_del_orb_kpts(k),k)
   enddo
   
END_PROVIDER
 
 BEGIN_PROVIDER [ integer, list_act_kpts        , (dim_list_act_orb_kpts,kpt_num) ]
&BEGIN_PROVIDER [ integer, list_act_kpts_reverse, (mo_num_per_kpt,kpt_num) ]
   implicit none
   BEGIN_DOC
   ! List of MO indices which are active.
   END_DOC
   integer                        :: i, n,k,di
   list_act_kpts = 0
   list_act_kpts_reverse = 0

   do k=1,kpt_num
     n=0
     di = (k-1)*mo_num_per_kpt
     do i = 1, mo_num_per_kpt
       if(mo_class(i+di) == 'Active')then
         n += 1
         list_act_kpts(n,k) = i
         list_act_kpts_reverse(i,k) = n
       endif
     enddo
     print *,  'Active MOs: ',k
     print *,  list_act_kpts(1:n_act_orb_kpts(k),k)
   enddo
   
END_PROVIDER

!todo: finish below for kpts
 
 BEGIN_PROVIDER [ integer, list_core_inact_kpts        , (dim_list_core_inact_orb_kpts,kpt_num) ]
&BEGIN_PROVIDER [ integer, list_core_inact_kpts_reverse, (mo_num_per_kpt,kpt_num) ]
  implicit none
  BEGIN_DOC
  ! List of indices of the core and inactive MOs
  END_DOC
  integer                        :: i,itmp,k
  list_core_inact_kpts_reverse = 0
  do k=1,kpt_num
    !call bitstring_to_list(reunion_of_core_inact_bitmask(1,1), list_core_inact, itmp, N_int)
    call bitstring_to_list(reunion_of_core_inact_bitmask_kpts(1,1,k), list_core_inact_kpts(1,k), itmp, N_int)
    ASSERT (itmp == n_core_inact_orb_kpts(k))
    do i = 1, n_core_inact_orb_kpts(k)
      list_core_inact_kpts_reverse(list_core_inact_kpts(i,k),k) = i
    enddo
    print *,  'Core and Inactive MOs: ',k
    print *,  list_core_inact_kpts(1:n_core_inact_orb_kpts(k),k)
  enddo
END_PROVIDER
 
 
 BEGIN_PROVIDER [ integer, list_core_inact_act_kpts        , (dim_list_core_inact_act_orb_kpts,kpt_num) ]
&BEGIN_PROVIDER [ integer, list_core_inact_act_kpts_reverse, (mo_num_per_kpt,kpt_num) ]
  implicit none
  BEGIN_DOC
  ! List of indices of the core inactive and active MOs
  END_DOC
  integer                        :: i,itmp,k
  list_core_inact_act_kpts_reverse = 0
  do k=1,kpt_num
    !call bitstring_to_list(reunion_of_core_inact_act_bitmask(1,1), list_core_inact_act, itmp, N_int)
    call bitstring_to_list(reunion_of_core_inact_act_bitmask_kpts(1,1,k), list_core_inact_act_kpts(1,k), itmp, N_int)
    ASSERT (itmp == n_core_inact_act_orb_kpts(k))
    do i = 1, n_core_inact_act_orb_kpts(k)
      list_core_inact_act_kpts_reverse(list_core_inact_act_kpts(i,k),k) = i
    enddo
    print *,  'Core, Inactive and Active MOs: ',k
    print *,  list_core_inact_act_kpts(1:n_core_inact_act_orb_kpts(k),k)
  enddo
END_PROVIDER
 
 
 BEGIN_PROVIDER [ integer, list_inact_act_kpts        , (dim_list_inact_act_orb_kpts,kpt_num) ]
&BEGIN_PROVIDER [ integer, list_inact_act_kpts_reverse, (mo_num_per_kpt,kpt_num) ]
   implicit none
   BEGIN_DOC
   ! List of indices of the inactive and active MOs
   END_DOC
   integer                        :: i,itmp,k
   list_inact_act_kpts_reverse = 0
   do k=1,kpt_num
     call bitstring_to_list(reunion_of_inact_act_bitmask_kpts(1,1,k), list_inact_act_kpts(1,k), itmp, N_int)
     ASSERT (itmp == n_inact_act_orb_kpts(k))
     do i = 1, n_inact_act_orb_kpts(k)
       list_inact_act_kpts_reverse(list_inact_act_kpts(i,k),k) = i
     enddo
     print *,  'Inactive and Active MOs: ',k
     print *,  list_inact_act_kpts(1:n_inact_act_orb_kpts(k),k)
   enddo
END_PROVIDER
 
