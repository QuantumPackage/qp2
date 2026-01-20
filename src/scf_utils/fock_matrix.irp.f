BEGIN_PROVIDER [ double precision, Fock_matrix_param, (2,3) ]
 implicit none
 BEGIN_DOC
 ! Fock matrix parameters.
 ! (1,1) : A cc  (2,1) : B cc
 ! (1,2) : A oo  (2,2) : B oo
 ! (1,3) : A vv  (2,3) : B vv
 END_DOC

 if (trim(rohf_parameters) == 'Roothaan') then
 ! C. C. J. Roothaan, Rev. Mod. Phys. 32, 179 ␁1960.

   Fock_matrix_param(1,1) = -0.5d0
   Fock_matrix_param(2,1) =  1.5d0
   Fock_matrix_param(1,2) =  0.5d0
   Fock_matrix_param(2,2) =  0.5d0
   Fock_matrix_param(1,3) =  1.5d0
   Fock_matrix_param(2,3) = -0.5d0

 else if (trim(rohf_parameters) == 'McWeeny') then
 ! R. McWeeny and G. Diercksen, J. Chem. Phys. 49, 4852 ␁1968

   Fock_matrix_param(1,1) =  1.d0/3.d0
   Fock_matrix_param(2,1) =  2.d0/3.d0
   Fock_matrix_param(1,2) =  1.d0/3.d0
   Fock_matrix_param(2,2) =  1.d0/3.d0
   Fock_matrix_param(1,3) =  2.d0/3.d0
   Fock_matrix_param(2,3) =  1.d0/3.d0

 else if (trim(rohf_parameters) == 'Davidson') then
 ! E. R. Davidson, Chem. Phys. Lett. 21, 565 ␁1973

   Fock_matrix_param(1,1) =  0.5d0
   Fock_matrix_param(2,1) =  0.5d0
   Fock_matrix_param(1,2) =  1.0d0
   Fock_matrix_param(2,2) =  0.0d0
   Fock_matrix_param(1,3) =  1.0d0
   Fock_matrix_param(2,3) =  0.0d0

 else if (trim(rohf_parameters) == 'Guest') then
 ! M. F. Guest and V. R. Saunders, Mol. Phys. 28, 819 ␁1974.

   Fock_matrix_param(1,1) =  0.5d0
   Fock_matrix_param(2,1) =  0.5d0
   Fock_matrix_param(1,2) =  0.5d0
   Fock_matrix_param(2,2) =  0.5d0
   Fock_matrix_param(1,3) =  0.0d0
   Fock_matrix_param(2,3) =  0.5d0

 else if (trim(rohf_parameters) == 'Binkley') then
 ! J. S. Binkley, J. A. Pople, and P. A. Dobosh, Mol. Phys. 28, 1423 ␁1974.

   Fock_matrix_param(1,1) =  0.5d0
   Fock_matrix_param(2,1) =  0.5d0
   Fock_matrix_param(1,2) =  1.0d0
   Fock_matrix_param(2,2) =  0.0d0
   Fock_matrix_param(1,3) =  0.0d0
   Fock_matrix_param(2,3) =  1.0d0

 else if (trim(rohf_parameters) == 'Faegri') then
 ! K. Faegri and R. Manne, Mol. Phys. 31, 1037 ␁1976.

   Fock_matrix_param(1,1) =  0.5d0
   Fock_matrix_param(2,1) =  0.5d0
   Fock_matrix_param(1,2) =  1.0d0
   Fock_matrix_param(2,2) =  0.0d0
   Fock_matrix_param(1,3) =  0.5d0
   Fock_matrix_param(2,3) =  0.5d0

 else if (trim(rohf_parameters) == 'Euler') then
 ! Plakhutin, B. N., et al.  J. Chem. Phys., . 125, 20, p. 204110, doi:10.1063/1.2393223.

   Fock_matrix_param(1,1) =  0.5d0
   Fock_matrix_param(2,1) =  0.5d0
   Fock_matrix_param(1,2) =  0.5d0
   Fock_matrix_param(2,2) =  0.0d0
   Fock_matrix_param(1,3) =  0.5d0
   Fock_matrix_param(2,3) =  0.5d0

 else if (trim(rohf_parameters) == 'Canonical') then
 ! Plakhutin, B. N., et al.  J. Chem. Phys., . 125, 20, p. 204110, doi:10.1063/1.2393223.

   Fock_matrix_param(1,1) =  0.0d0
   Fock_matrix_param(2,1) =  1.0d0
   Fock_matrix_param(1,2) =  1.0d0
   Fock_matrix_param(2,2) =  0.0d0
   Fock_matrix_param(1,3) =  1.0d0
   Fock_matrix_param(2,3) =  0.0d0

 else

  stop 'Unknown set of ROHF parameters'

 endif


END_PROVIDER

 BEGIN_PROVIDER [ double precision, Fock_matrix_mo, (mo_num,mo_num) ]
&BEGIN_PROVIDER [ double precision, Fock_matrix_diag_mo, (mo_num)]
   implicit none
   BEGIN_DOC
   ! Fock matrix on the MO basis.
   ! For open shells, the ROHF Fock Matrix is ::
   !
   !       |  Rcc  |  F^b  |  Fcv  |
   !       |-----------------------|
   !       |  F^b  |  Roo  |  F^a  |
   !       |-----------------------|
   !       |  Fcv  |  F^a  |  Rvv  |
   !
   ! C: Core, O: Open, V: Virtual
   !
   ! Rcc = Acc Fcc^a + Bcc Fcc^b
   ! Roo = Aoo Foo^a + Boo Foo^b
   ! Rvv = Avv Fvv^a + Bvv Fvv^b
   ! Fcv = (F^a + F^b)/2
   !
   ! F^a: Fock matrix alpha (MO), F^b: Fock matrix beta (MO)
   ! A,B: Coupling parameters
   !
   ! J. Chem. Phys. 133, 141102 (2010), https://doi.org/10.1063/1.3503173
   ! Coupling parameters from J. Chem. Phys. 125, 204110 (2006); https://doi.org/10.1063/1.2393223.
   !       cc   oo   vv
   !  A  -0.5  0.5  1.5
   !  B   1.5  0.5 -0.5
   !
   END_DOC
   integer                        :: i,j,n
   if (all_shells_closed) then
     Fock_matrix_mo = Fock_matrix_mo_alpha
   else
     ! Core
     do j = 1, elec_beta_num
       ! Core
       do i = 1, elec_beta_num
         fock_matrix_mo(i,j) = Fock_matrix_param(1,1) * fock_matrix_mo_alpha(i,j) + &
                               Fock_matrix_param(2,1) * fock_matrix_mo_beta(i,j)
       enddo
       ! Open
       do i = elec_beta_num+1, elec_alpha_num
         fock_matrix_mo(i,j) = fock_matrix_mo_beta(i,j)
       enddo
       ! Virtual
       do i = elec_alpha_num+1, mo_num
         fock_matrix_mo(i,j) = 0.5d0 * fock_matrix_mo_alpha(i,j) + &
                               0.5d0 * fock_matrix_mo_beta(i,j)
       enddo
     enddo
     ! Open
     do j = elec_beta_num+1, elec_alpha_num
       ! Core
       do i = 1, elec_beta_num
         fock_matrix_mo(i,j) = fock_matrix_mo_beta(i,j)
       enddo
       ! Open
       do i = elec_beta_num+1, elec_alpha_num
         fock_matrix_mo(i,j) = Fock_matrix_param(1,2) * fock_matrix_mo_alpha(i,j) + &
                               Fock_matrix_param(2,2) * fock_matrix_mo_beta(i,j)
       enddo
       ! Virtual
       do i = elec_alpha_num+1, mo_num
         fock_matrix_mo(i,j) =  fock_matrix_mo_alpha(i,j)
       enddo
     enddo
     ! Virtual
     do j = elec_alpha_num+1, mo_num
       ! Core
       do i = 1, elec_beta_num
         fock_matrix_mo(i,j) =   0.5d0 * fock_matrix_mo_alpha(i,j) &
                               + 0.5d0 * fock_matrix_mo_beta(i,j)
       enddo
       ! Open
       do i = elec_beta_num+1, elec_alpha_num
         fock_matrix_mo(i,j) = fock_matrix_mo_alpha(i,j)
       enddo
       ! Virtual
       do i = elec_alpha_num+1, mo_num
         fock_matrix_mo(i,j) = Fock_matrix_param(1,3) * fock_matrix_mo_alpha(i,j) + &
                               Fock_matrix_param(2,3) * fock_matrix_mo_beta(i,j)
       enddo
     enddo
   endif


   do i = 1, mo_num
     Fock_matrix_diag_mo(i) = Fock_matrix_mo(i,i)
   enddo


   if(frozen_orb_scf)then
     integer                        :: iorb,jorb
     !       active|core|active
     !active |     | 0  |
     !core   |  0  |    |   0
     !active |     | 0  |
     do i = 1, n_core_orb
      iorb = list_core(i)
      do j = 1, n_act_orb
       jorb = list_act(j)
       Fock_matrix_mo(iorb,jorb) = 0.d0
       Fock_matrix_mo(jorb,iorb) = 0.d0
      enddo
     enddo
   endif

   if(no_oa_or_av_opt)then
     do i = 1, n_act_orb
       iorb = list_act(i)
       do j = 1, n_inact_orb
         jorb = list_inact(j)
         Fock_matrix_mo(iorb,jorb) = 0.d0
         Fock_matrix_mo(jorb,iorb) = 0.d0
       enddo
       do j = 1, n_virt_orb
         jorb = list_virt(j)
         Fock_matrix_mo(iorb,jorb) = 0.d0
         Fock_matrix_mo(jorb,iorb) = 0.d0
       enddo
       do j = 1, n_core_orb
         jorb = list_core(j)
         Fock_matrix_mo(iorb,jorb) = 0.d0
         Fock_matrix_mo(jorb,iorb) = 0.d0
       enddo
     enddo
   endif

END_PROVIDER

BEGIN_PROVIDER [ double precision, Fock_matrix_mo_alpha, (mo_num,mo_num) ]
   implicit none
   BEGIN_DOC
   ! Fock matrix on the MO basis
   END_DOC
   call ao_to_mo(Fock_matrix_ao_alpha,size(Fock_matrix_ao_alpha,1), &
                 Fock_matrix_mo_alpha,size(Fock_matrix_mo_alpha,1))
END_PROVIDER


BEGIN_PROVIDER [ double precision, Fock_matrix_mo_beta, (mo_num,mo_num) ]
   implicit none
   BEGIN_DOC
   ! Fock matrix on the MO basis
   END_DOC
   call ao_to_mo(Fock_matrix_ao_beta,size(Fock_matrix_ao_beta,1), &
                 Fock_matrix_mo_beta,size(Fock_matrix_mo_beta,1))
END_PROVIDER


BEGIN_PROVIDER [ double precision, Fock_matrix_ao, (ao_num, ao_num) ]
 implicit none
 BEGIN_DOC
 ! Fock matrix in AO basis set
 END_DOC

 if(frozen_orb_scf)then
   call mo_to_ao(Fock_matrix_mo,size(Fock_matrix_mo,1),              &
       Fock_matrix_ao,size(Fock_matrix_ao,1))
 else
   if (all_shells_closed.and. (level_shift == 0.)) then
     integer                        :: i,j
     do j=1,ao_num
       do i=1,ao_num
         Fock_matrix_ao(i,j) = Fock_matrix_ao_alpha(i,j)
       enddo
     enddo
   else
     call mo_to_ao(Fock_matrix_mo,size(Fock_matrix_mo,1),            &
         Fock_matrix_ao,size(Fock_matrix_ao,1))
   endif
 endif
END_PROVIDER


BEGIN_PROVIDER [ double precision, SCF_energy ]
 implicit none
 BEGIN_DOC
 ! Hartree-Fock energy
 END_DOC
 integer                        :: i,j

 SCF_energy = 0.d0

 do j=1,ao_num
   do i=1,ao_num
     SCF_energy +=  &
         (ao_one_e_integrals(i,j) + Fock_matrix_ao_alpha(i,j) ) *  SCF_density_matrix_ao_alpha(i,j) +&
         (ao_one_e_integrals(i,j) + Fock_matrix_ao_beta (i,j) ) *  SCF_density_matrix_ao_beta (i,j)
   enddo
 enddo
 SCF_energy = 0.5d0 * SCF_energy + extra_e_contrib_density + nuclear_repulsion

END_PROVIDER

