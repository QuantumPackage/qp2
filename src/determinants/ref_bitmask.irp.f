 BEGIN_PROVIDER [ double precision, ref_bitmask_energy ]
&BEGIN_PROVIDER [ double precision, ref_bitmask_one_e_energy ]
&BEGIN_PROVIDER [ double precision, ref_bitmask_kinetic_energy ]
&BEGIN_PROVIDER [ double precision, ref_bitmask_n_e_energy ]
&BEGIN_PROVIDER [ double precision, ref_bitmask_two_e_energy ]
&BEGIN_PROVIDER [ double precision, ref_bitmask_energy_ab ]
&BEGIN_PROVIDER [ double precision, ref_bitmask_energy_bb ]
&BEGIN_PROVIDER [ double precision, ref_bitmask_energy_aa ]

  use bitmasks
  implicit none
  BEGIN_DOC
  ! Energy of the reference bitmask used in Slater rules
  END_DOC

  integer                        :: occ(N_int*bit_kind_size,2)
  integer                        :: i,j

  call bitstring_to_list(ref_bitmask(1,1), occ(1,1), i, N_int)
  call bitstring_to_list(ref_bitmask(1,2), occ(1,2), i, N_int)


  ref_bitmask_energy = 0.d0
  ref_bitmask_one_e_energy = 0.d0
  ref_bitmask_kinetic_energy   = 0.d0
  ref_bitmask_n_e_energy = 0.d0
  ref_bitmask_two_e_energy = 0.d0

  do i = 1, elec_beta_num
    ref_bitmask_energy += mo_one_e_integrals_diag(occ(i,1)) + mo_one_e_integrals_diag(occ(i,2))
    ref_bitmask_kinetic_energy += mo_kinetic_integrals_diag(occ(i,1)) + mo_kinetic_integrals_diag(occ(i,2))
    ref_bitmask_n_e_energy += mo_integrals_n_e_diag(occ(i,1)) + mo_integrals_n_e_diag(occ(i,2))
  enddo

  do i = elec_beta_num+1,elec_alpha_num
    ref_bitmask_energy += mo_one_e_integrals_diag(occ(i,1))
    ref_bitmask_kinetic_energy += mo_kinetic_integrals_diag(occ(i,1))
    ref_bitmask_n_e_energy += mo_integrals_n_e_diag(occ(i,1))
  enddo

  do j= 1, elec_alpha_num
    do i = j+1, elec_alpha_num
      ref_bitmask_two_e_energy += mo_two_e_integrals_jj_anti(occ(i,1),occ(j,1))
      ref_bitmask_energy += mo_two_e_integrals_jj_anti(occ(i,1),occ(j,1))
    enddo
  enddo

  do j= 1, elec_beta_num
    do i = j+1, elec_beta_num
      ref_bitmask_two_e_energy += mo_two_e_integrals_jj_anti(occ(i,2),occ(j,2))
      ref_bitmask_energy += mo_two_e_integrals_jj_anti(occ(i,2),occ(j,2))
    enddo
    do i= 1, elec_alpha_num
      ref_bitmask_two_e_energy += mo_two_e_integrals_jj(occ(i,1),occ(j,2))
      ref_bitmask_energy += mo_two_e_integrals_jj(occ(i,1),occ(j,2))
    enddo
  enddo
  ref_bitmask_one_e_energy = ref_bitmask_kinetic_energy +   ref_bitmask_n_e_energy

 ref_bitmask_energy_ab = 0.d0
 do i = 1, elec_alpha_num
  do j = 1, elec_beta_num
   ref_bitmask_energy_ab += mo_two_e_integrals_jj(occ(i,1),occ(j,2))
  enddo
 enddo

 ref_bitmask_energy_aa = 0.d0
 do i = 1, elec_alpha_num
  do j = 1, elec_alpha_num
   ref_bitmask_energy_aa += mo_two_e_integrals_jj_anti(occ(i,1),occ(j,1))
  enddo
 enddo
 ref_bitmask_energy_aa = ref_bitmask_energy_aa * 0.5d0

 ref_bitmask_energy_bb = 0.d0
 do i = 1, elec_beta_num
  do j = 1, elec_beta_num
   ref_bitmask_energy_bb += mo_two_e_integrals_jj_anti(occ(i,2),occ(j,2))
  enddo
 enddo
 ref_bitmask_energy_bb = ref_bitmask_energy_bb * 0.5d0



END_PROVIDER

