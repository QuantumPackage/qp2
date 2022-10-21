
! ---

BEGIN_PROVIDER [double precision, ao_one_e_integrals_tc_tot, (ao_num,ao_num)]

  implicit none
  integer :: i, j

  ao_one_e_integrals_tc_tot = ao_one_e_integrals      

  provide j1b_type

  if( (j1b_type .eq. 1) .or. (j1b_type .eq. 2) ) then

    do i = 1, ao_num
      do j = 1, ao_num
        ao_one_e_integrals_tc_tot(j,i) += ( j1b_gauss_hermI  (j,i) &
                                          + j1b_gauss_hermII (j,i) &
                                          + j1b_gauss_nonherm(j,i) )
      enddo
    enddo

  endif

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, mo_bi_ortho_tc_one_e, (mo_num, mo_num)]

  BEGIN_DOC
  !
  ! mo_bi_ortho_tc_one_e(k,i) = <MO^L_k | h_c | MO^R_i>
  !
  END_DOC

  implicit none
 
  call ao_to_mo_bi_ortho(ao_one_e_integrals_tc_tot, ao_num, mo_bi_ortho_tc_one_e, mo_num)

END_PROVIDER 

! ---

 BEGIN_PROVIDER [double precision, mo_bi_orth_bipole_x , (mo_num,mo_num)]
&BEGIN_PROVIDER [double precision, mo_bi_orth_bipole_y , (mo_num,mo_num)]
&BEGIN_PROVIDER [double precision, mo_bi_orth_bipole_z , (mo_num,mo_num)]
 BEGIN_DOC
 ! array of the integrals of MO_i * x MO_j
 ! array of the integrals of MO_i * y MO_j
 ! array of the integrals of MO_i * z MO_j
 END_DOC
 implicit none

  call ao_to_mo_bi_ortho(                                                     &
      ao_dipole_x,                                                   &
      size(ao_dipole_x,1),                                           &
      mo_bi_orth_bipole_x,                                                   &
      size(mo_bi_orth_bipole_x,1)                                            &
      )
  call ao_to_mo_bi_ortho(                                                     &
      ao_dipole_y,                                                   &
      size(ao_dipole_y,1),                                           &
      mo_bi_orth_bipole_y,                                                   &
      size(mo_bi_orth_bipole_y,1)                                            &
      )
  call ao_to_mo_bi_ortho(                                                     &
      ao_dipole_z,                                                   &
      size(ao_dipole_z,1),                                           &
      mo_bi_orth_bipole_z,                                                   &
      size(mo_bi_orth_bipole_z,1)                                            &
      )

END_PROVIDER

