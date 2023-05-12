! Dimensions of MOs


BEGIN_PROVIDER [ integer, n_mo_dim ]
  implicit none
  BEGIN_DOC
  ! Number of different pairs (i,j) of MOs we can build,
  ! with i>j
  END_DOC

  n_mo_dim = mo_num*(mo_num-1)/2

END_PROVIDER

BEGIN_PROVIDER [ integer, n_mo_dim_core ]
  implicit none 
  BEGIN_DOC
  ! Number of different pairs (i,j) of core MOs we can build,
  ! with i>j
  END_DOC

  n_mo_dim_core = dim_list_core_orb*(dim_list_core_orb-1)/2

END_PROVIDER

BEGIN_PROVIDER [ integer, n_mo_dim_act ]
  implicit none
  BEGIN_DOC
  ! Number of different pairs (i,j) of active MOs we can build,
  ! with i>j
  END_DOC

  n_mo_dim_act = dim_list_act_orb*(dim_list_act_orb-1)/2

END_PROVIDER

BEGIN_PROVIDER [ integer, n_mo_dim_inact ]
  implicit none 
  BEGIN_DOC
  ! Number of different pairs (i,j) of inactive MOs we can build,
  ! with i>j
  END_DOC

  n_mo_dim_inact = dim_list_inact_orb*(dim_list_inact_orb-1)/2

END_PROVIDER

BEGIN_PROVIDER [ integer, n_mo_dim_virt ]
  implicit none 
  BEGIN_DOC
  ! Number of different pairs (i,j) of virtual MOs we can build,
  ! with i>j
  END_DOC

  n_mo_dim_virt = dim_list_virt_orb*(dim_list_virt_orb-1)/2

END_PROVIDER

! Energies/criterions

BEGIN_PROVIDER [ double precision, my_st_av_energy ]
  implicit none
  BEGIN_DOC
  ! State average CI energy
  END_DOC

  !call update_st_av_ci_energy(my_st_av_energy)
  call state_average_energy(my_st_av_energy)

END_PROVIDER

! With all the MOs

BEGIN_PROVIDER [ double precision, my_gradient_opt, (n_mo_dim) ]
&BEGIN_PROVIDER [ double precision, my_CC1_opt ]
  implicit none
  BEGIN_DOC
  ! - Gradient of the energy with respect to the MO rotations, for all the MOs.
  ! - Maximal element of the gradient in absolute value 
  END_DOC

  double precision :: norm_grad

  PROVIDE mo_two_e_integrals_in_map

  call gradient_opt(n_mo_dim, my_gradient_opt, my_CC1_opt, norm_grad)

END_PROVIDER

BEGIN_PROVIDER [ double precision, my_hessian_opt, (n_mo_dim, n_mo_dim) ]
  implicit none
  BEGIN_DOC
  ! - Gradient of the energy with respect to the MO rotations, for all the MOs.
  ! - Maximal element of the gradient in absolute value 
  END_DOC

  double precision, allocatable :: h_f(:,:,:,:)

  PROVIDE mo_two_e_integrals_in_map

  allocate(h_f(mo_num, mo_num, mo_num, mo_num))

  call hessian_list_opt(n_mo_dim, my_hessian_opt, h_f)

END_PROVIDER

! With the list of active MOs
! Can be generalized to any mo_class by changing the list/dimension

BEGIN_PROVIDER [ double precision, my_gradient_list_opt, (n_mo_dim_act) ]
&BEGIN_PROVIDER [ double precision, my_CC2_opt ]
  implicit none
  BEGIN_DOC
  ! - Gradient of the energy with respect to the MO rotations, only for the active MOs !
  ! - Maximal element of the gradient in absolute value 
  END_DOC

  double precision :: norm_grad

  PROVIDE mo_two_e_integrals_in_map !one_e_dm_mo two_e_dm_mo mo_one_e_integrals 

  call gradient_list_opt(n_mo_dim_act, dim_list_act_orb, list_act, my_gradient_list_opt, my_CC2_opt, norm_grad)

END_PROVIDER

BEGIN_PROVIDER [ double precision, my_hessian_list_opt, (n_mo_dim_act, n_mo_dim_act) ]
  implicit none
  BEGIN_DOC
  ! - Gradient of the energy with respect to the MO rotations, only for the active MOs !
  ! - Maximal element of the gradient in absolute value 
  END_DOC

  double precision, allocatable :: h_f(:,:,:,:)

  PROVIDE mo_two_e_integrals_in_map

  allocate(h_f(dim_list_act_orb, dim_list_act_orb, dim_list_act_orb, dim_list_act_orb))

  call hessian_list_opt(n_mo_dim_act, dim_list_act_orb, list_act, my_hessian_list_opt, h_f)

END_PROVIDER
