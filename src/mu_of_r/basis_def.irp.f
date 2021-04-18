

 BEGIN_PROVIDER [integer, n_occ_val_orb_for_hf,(2)]
&BEGIN_PROVIDER [integer, n_max_occ_val_orb_for_hf]
 implicit none
 BEGIN_DOC
 ! Number of OCCUPIED VALENCE ORBITALS for each spin to build the f_{HF}(r_1,r_2) function 
 !
 ! This is typically elec_alpha_num - n_core_orb for alpha electrons and elec_beta_num - n_core_orb for beta electrons
 !
 ! This determines the size of the space \mathcal{A} of Eqs. (15-16) of Phys.Chem.Lett.2019, 10, 2931   2937
 END_DOC
 integer :: i
 n_occ_val_orb_for_hf = 0
 ! You browse the ALPHA ELECTRONS and check if its not a CORE ORBITAL 
 do i = 1, elec_alpha_num
  if(  trim(mo_class(i))=="Inactive" & 
  .or. trim(mo_class(i))=="Active"   &
  .or. trim(mo_class(i))=="Virtual" )then
   n_occ_val_orb_for_hf(1) +=1 
  endif
 enddo

 ! You browse the BETA  ELECTRONS and check if its not a CORE ORBITAL 
 do i = 1, elec_beta_num
  if(  trim(mo_class(i))=="Inactive" & 
  .or. trim(mo_class(i))=="Active"   &
  .or. trim(mo_class(i))=="Virtual" )then
   n_occ_val_orb_for_hf(2) +=1 
  endif
 enddo
 n_max_occ_val_orb_for_hf = maxval(n_occ_val_orb_for_hf)

END_PROVIDER 

BEGIN_PROVIDER [integer, list_valence_orb_for_hf, (n_max_occ_val_orb_for_hf,2)]
 implicit none
 BEGIN_DOC
 ! List of OCCUPIED valence orbitals for each spin to build the f_{HF}(r_1,r_2) function 
 !
 ! This corresponds to ALL OCCUPIED orbitals in the HF wave function, except those defined as "core" 
 !
 ! This determines the space \mathcal{A} of Eqs. (15-16) of Phys.Chem.Lett.2019, 10, 2931   2937
 END_DOC
 integer :: i,j
 j = 0
 ! You browse the ALPHA ELECTRONS and check if its not a CORE ORBITAL 
 do i = 1, elec_alpha_num
  if(  trim(mo_class(i))=="Inactive" & 
  .or. trim(mo_class(i))=="Active"   &
  .or. trim(mo_class(i))=="Virtual" )then
   j +=1 
   list_valence_orb_for_hf(j,1) = i
  endif
 enddo

 j = 0
 ! You browse the BETA  ELECTRONS and check if its not a CORE ORBITAL 
 do i = 1, elec_beta_num
  if(  trim(mo_class(i))=="Inactive" & 
  .or. trim(mo_class(i))=="Active"   &
  .or. trim(mo_class(i))=="Virtual" )then
   j +=1 
   list_valence_orb_for_hf(j,2) = i
  endif
 enddo

END_PROVIDER 

BEGIN_PROVIDER [integer, n_basis_orb]
 implicit none
 BEGIN_DOC
 ! Defines the number of orbitals you will use to explore the basis set 
 !
 ! This determines the size of the space \mathcal{B} of Eqs. (15-16) of Phys.Chem.Lett.2019, 10, 2931   2937
 !
 ! It corresponds to all MOs except those defined as "deleted" 
 END_DOC
 n_basis_orb = n_all_but_del_orb
END_PROVIDER 

BEGIN_PROVIDER [integer, list_basis, (n_basis_orb)]
 implicit none
 BEGIN_DOC
 ! Defines the set of orbitals you will use to explore the basis set 
 !
 ! This determines the space \mathcal{B} of Eqs. (15-16) of Phys.Chem.Lett.2019, 10, 2931   2937
 !
 ! It corresponds to all MOs except those defined as "deleted" 
 END_DOC
 integer :: i
 do i = 1, n_all_but_del_orb
   list_basis(i) = list_all_but_del_orb(i)
 enddo
END_PROVIDER 

BEGIN_PROVIDER [double precision, basis_mos_in_r_array, (n_basis_orb,n_points_final_grid)]
 implicit none
 integer :: ipoint,i,ii
 do ipoint = 1, n_points_final_grid
  do i = 1, n_basis_orb 
   ii = list_basis(i)
   basis_mos_in_r_array(i,ipoint) = mos_in_r_array(ii,ipoint)
  enddo
 enddo
END_PROVIDER 
