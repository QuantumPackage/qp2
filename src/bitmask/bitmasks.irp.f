use bitmasks

BEGIN_PROVIDER [ integer, N_int ]
  implicit none
  include 'utils/constants.include.F'
  BEGIN_DOC
  ! Number of 64-bit integers needed to represent determinants as binary strings
  END_DOC
  N_int = (mo_num-1)/bit_kind_size + 1
  call write_int(6,N_int, 'N_int')
  if (N_int > N_int_max) then
    stop 'N_int > N_int_max'
  endif
  
END_PROVIDER


BEGIN_PROVIDER [ integer(bit_kind), full_ijkl_bitmask, (N_int) ]
  implicit none
  BEGIN_DOC
  ! Bitmask to include all possible MOs
  END_DOC
  
  integer                        :: i,j,k
  k=0
  do j=1,N_int
    full_ijkl_bitmask(j) = 0_bit_kind
    do i=0,bit_kind_size-1
      k=k+1
      if (mo_class(k) /= 'Deleted') then
        full_ijkl_bitmask(j) = ibset(full_ijkl_bitmask(j),i)
      endif
      if (k == mo_num) exit
    enddo
  enddo
END_PROVIDER

BEGIN_PROVIDER [ integer(bit_kind), full_ijkl_bitmask_4, (N_int,4) ]
  implicit none
  integer                        :: i
  do i=1,N_int
    full_ijkl_bitmask_4(i,1) = full_ijkl_bitmask(i)
    full_ijkl_bitmask_4(i,2) = full_ijkl_bitmask(i)
    full_ijkl_bitmask_4(i,3) = full_ijkl_bitmask(i)
    full_ijkl_bitmask_4(i,4) = full_ijkl_bitmask(i)
  enddo
END_PROVIDER

BEGIN_PROVIDER [ integer(bit_kind), core_inact_act_bitmask_4, (N_int,4) ]
  implicit none
  integer                        :: i
  do i=1,N_int
    core_inact_act_bitmask_4(i,1) = reunion_of_core_inact_act_bitmask(i,1)
    core_inact_act_bitmask_4(i,2) = reunion_of_core_inact_act_bitmask(i,1)
    core_inact_act_bitmask_4(i,3) = reunion_of_core_inact_act_bitmask(i,1)
    core_inact_act_bitmask_4(i,4) = reunion_of_core_inact_act_bitmask(i,1)
  enddo
END_PROVIDER

BEGIN_PROVIDER [ integer(bit_kind), virt_bitmask_4, (N_int,4) ]
  implicit none
  integer                        :: i
  do i=1,N_int
    virt_bitmask_4(i,1) = virt_bitmask(i,1)
    virt_bitmask_4(i,2) = virt_bitmask(i,1)
    virt_bitmask_4(i,3) = virt_bitmask(i,1)
    virt_bitmask_4(i,4) = virt_bitmask(i,1)
  enddo
END_PROVIDER




BEGIN_PROVIDER [ integer(bit_kind), HF_bitmask, (N_int,2)]
  implicit none
  BEGIN_DOC
  ! Hartree Fock bit mask
  END_DOC
  integer                        :: i,j,n
  integer                        :: occ(elec_alpha_num)
  
  HF_bitmask = 0_bit_kind
  !if (is_complex) then
  if (.False.) then
    integer :: kpt,korb
    kpt=1
    korb=1
    do i=1,elec_alpha_num
      occ(i) = korb + (kpt-1) * ao_num_per_kpt
      kpt += 1
      if (kpt > kpt_num) then
        kpt = 1
        korb += 1
      endif
    enddo
  else
    do i=1,elec_alpha_num
      occ(i) = i
    enddo
  endif
  call list_to_bitstring( HF_bitmask(1,1), occ, elec_alpha_num, N_int)
  ! elec_alpha_num <= elec_beta_num, so occ is already OK.
  call list_to_bitstring( HF_bitmask(1,2), occ, elec_beta_num, N_int)
  
END_PROVIDER

BEGIN_PROVIDER [ integer(bit_kind), ref_bitmask, (N_int,2)]
  implicit none
  BEGIN_DOC
  ! Reference bit mask, used in Slater rules, chosen as Hartree-Fock bitmask
  END_DOC
  ref_bitmask = HF_bitmask
END_PROVIDER



BEGIN_PROVIDER [ integer(bit_kind), generators_bitmask, (N_int,2,6) ]
  implicit none
  BEGIN_DOC
  ! Bitmasks for generator determinants.
  ! (N_int, alpha/beta, hole/particle, generator).
  !
  ! 3rd index is :
  !
  ! * 1 : hole     for single exc
  !
  ! * 2 : particle for single exc
  !
  ! * 3 : hole     for 1st exc of double
  !
  ! * 4 : particle for 1st exc of double
  !
  ! * 5 : hole     for 2nd exc of double
  !
  ! * 6 : particle for 2nd exc of double
  !
  END_DOC
  logical                        :: exists
  PROVIDE ezfio_filename full_ijkl_bitmask 
  
  integer                        :: ispin, i
  do ispin=1,2
      do i=1,N_int
        generators_bitmask(i,ispin,s_hole ) = reunion_of_inact_act_bitmask(i,ispin)
        generators_bitmask(i,ispin,s_part ) = reunion_of_act_virt_bitmask(i,ispin)
        generators_bitmask(i,ispin,d_hole1) = reunion_of_inact_act_bitmask(i,ispin)
        generators_bitmask(i,ispin,d_part1) = reunion_of_act_virt_bitmask(i,ispin)
        generators_bitmask(i,ispin,d_hole2) = reunion_of_inact_act_bitmask(i,ispin)
        generators_bitmask(i,ispin,d_part2) = reunion_of_act_virt_bitmask(i,ispin)
      enddo
  enddo
  
END_PROVIDER

BEGIN_PROVIDER [ integer(bit_kind), reunion_of_core_inact_bitmask, (N_int,2)]
  implicit none
  BEGIN_DOC
  ! Reunion of the core and inactive and virtual bitmasks
  END_DOC
  integer                        :: i
  do i = 1, N_int
    reunion_of_core_inact_bitmask(i,1) = ior(core_bitmask(i,1),inact_bitmask(i,1))
    reunion_of_core_inact_bitmask(i,2) = ior(core_bitmask(i,2),inact_bitmask(i,2))
  enddo
END_PROVIDER


BEGIN_PROVIDER [integer(bit_kind), reunion_of_inact_act_bitmask, (N_int,2)]
  implicit none
  BEGIN_DOC
  ! Reunion of the  inactive and active bitmasks
  END_DOC
  integer                        :: i,j
  
  do i = 1, N_int
    reunion_of_inact_act_bitmask(i,1) = ior(inact_bitmask(i,1),act_bitmask(i,1))
    reunion_of_inact_act_bitmask(i,2) = ior(inact_bitmask(i,2),act_bitmask(i,2))
  enddo
END_PROVIDER

BEGIN_PROVIDER [integer(bit_kind), reunion_of_act_virt_bitmask, (N_int,2)]
  implicit none
  BEGIN_DOC
  ! Reunion of the  inactive and active bitmasks
  END_DOC
  integer                        :: i,j
  
  do i = 1, N_int
    reunion_of_act_virt_bitmask(i,1) = ior(virt_bitmask(i,1),act_bitmask(i,1))
    reunion_of_act_virt_bitmask(i,2) = ior(virt_bitmask(i,2),act_bitmask(i,2))
  enddo
END_PROVIDER


BEGIN_PROVIDER [integer(bit_kind), reunion_of_core_inact_act_bitmask, (N_int,2)]
  implicit none
  BEGIN_DOC
  ! Reunion of the core, inactive and active bitmasks
  END_DOC
  integer                        :: i,j
  
  do i = 1, N_int
    reunion_of_core_inact_act_bitmask(i,1) = ior(reunion_of_core_inact_bitmask(i,1),act_bitmask(i,1))
    reunion_of_core_inact_act_bitmask(i,2) = ior(reunion_of_core_inact_bitmask(i,2),act_bitmask(i,2))
  enddo
END_PROVIDER


BEGIN_PROVIDER [ integer(bit_kind), reunion_of_bitmask, (N_int,2)]
  implicit none
  BEGIN_DOC
  ! Reunion of the inactive, active and virtual bitmasks
  END_DOC
  integer                        :: i,j
  do i = 1, N_int
    reunion_of_bitmask(i,1) = ior(ior(act_bitmask(i,1),inact_bitmask(i,1)),virt_bitmask(i,1))
    reunion_of_bitmask(i,2) = ior(ior(act_bitmask(i,2),inact_bitmask(i,2)),virt_bitmask(i,2))
  enddo
END_PROVIDER


 BEGIN_PROVIDER [ integer(bit_kind), inact_virt_bitmask, (N_int,2)]
&BEGIN_PROVIDER [ integer(bit_kind), core_inact_virt_bitmask, (N_int,2)]
  implicit none
  BEGIN_DOC
  ! Reunion of the inactive and virtual bitmasks
  END_DOC
  integer                        :: i,j
  do i = 1, N_int
    inact_virt_bitmask(i,1) = ior(inact_bitmask(i,1),virt_bitmask(i,1))
    inact_virt_bitmask(i,2) = ior(inact_bitmask(i,2),virt_bitmask(i,2))
    core_inact_virt_bitmask(i,1) = ior(core_bitmask(i,1),inact_virt_bitmask(i,1))
    core_inact_virt_bitmask(i,2) = ior(core_bitmask(i,2),inact_virt_bitmask(i,2))
  enddo
END_PROVIDER


BEGIN_PROVIDER [ integer(bit_kind), unpaired_alpha_electrons, (N_int)]
  implicit none
  BEGIN_DOC
  ! Bitmask reprenting the unpaired alpha electrons in the HF_bitmask
  END_DOC
  integer                        :: i
  unpaired_alpha_electrons = 0_bit_kind
  do i = 1, N_int
    unpaired_alpha_electrons(i) = xor(HF_bitmask(i,1),HF_bitmask(i,2))
  enddo
END_PROVIDER

BEGIN_PROVIDER [integer(bit_kind), closed_shell_ref_bitmask, (N_int,2)]
  implicit none
  integer                        :: i,j
  do i = 1, N_int
    closed_shell_ref_bitmask(i,1) = ior(ref_bitmask(i,1),act_bitmask(i,1))
    closed_shell_ref_bitmask(i,2) = ior(ref_bitmask(i,2),act_bitmask(i,2))
  enddo
END_PROVIDER
