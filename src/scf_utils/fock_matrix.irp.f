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
   if (elec_alpha_num == elec_beta_num) then
     Fock_matrix_mo = Fock_matrix_mo_alpha
   else
     ! Core
     do j = 1, elec_beta_num
       ! Core
       do i = 1, elec_beta_num
         fock_matrix_mo(i,j) = - 0.5d0 * fock_matrix_mo_alpha(i,j) &
                               + 1.5d0 * fock_matrix_mo_beta(i,j)
       enddo
       ! Open
       do i = elec_beta_num+1, elec_alpha_num
         fock_matrix_mo(i,j) = fock_matrix_mo_beta(i,j)
       enddo
       ! Virtual
       do i = elec_alpha_num+1, mo_num
         fock_matrix_mo(i,j) =   0.5d0 * fock_matrix_mo_alpha(i,j) &
                               + 0.5d0 * fock_matrix_mo_beta(i,j)
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
         fock_matrix_mo(i,j) =   0.5d0 * fock_matrix_mo_alpha(i,j) &
                               + 0.5d0 * fock_matrix_mo_beta(i,j)
       enddo
       ! Virtual
       do i = elec_alpha_num+1, mo_num
         fock_matrix_mo(i,j) = fock_matrix_mo_alpha(i,j)
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
         fock_matrix_mo(i,j) =   1.5d0 * fock_matrix_mo_alpha(i,j) &
                               - 0.5d0 * fock_matrix_mo_beta(i,j)
       enddo
     enddo
   endif

   ! Old
   ! BEGIN_DOC
   ! Fock matrix on the MO basis.
   ! For open shells, the ROHF Fock Matrix is ::
   !
   !       |   F-K    |  F + K/2  |    F     |
   !       |---------------------------------|
   !       | F + K/2  |     F     |  F - K/2 |
   !       |---------------------------------|
   !       |    F     |  F - K/2  |  F + K   |
   !
   !
   ! F = 1/2 (Fa + Fb)
   !
   ! K = Fb - Fa
   !
   ! END_DOC
   !integer                        :: i,j,n
   !if (elec_alpha_num == elec_beta_num) then
   !  Fock_matrix_mo = Fock_matrix_mo_alpha
   !else

   !  do j=1,elec_beta_num
   !    ! F-K
   !    do i=1,elec_beta_num !CC
   !      Fock_matrix_mo(i,j) = 0.5d0*(Fock_matrix_mo_alpha(i,j)+Fock_matrix_mo_beta(i,j))&
   !          - (Fock_matrix_mo_beta(i,j) - Fock_matrix_mo_alpha(i,j))
   !    enddo
   !    ! F+K/2
   !    do i=elec_beta_num+1,elec_alpha_num  !CA
   !      Fock_matrix_mo(i,j) = 0.5d0*(Fock_matrix_mo_alpha(i,j)+Fock_matrix_mo_beta(i,j))&
   !          + 0.5d0*(Fock_matrix_mo_beta(i,j) - Fock_matrix_mo_alpha(i,j))
   !    enddo
   !    ! F
   !    do i=elec_alpha_num+1, mo_num !CV
   !      Fock_matrix_mo(i,j) = 0.5d0*(Fock_matrix_mo_alpha(i,j)+Fock_matrix_mo_beta(i,j))
   !    enddo
   !  enddo

   !  do j=elec_beta_num+1,elec_alpha_num
   !    ! F+K/2
   !    do i=1,elec_beta_num !AC
   !      Fock_matrix_mo(i,j) = 0.5d0*(Fock_matrix_mo_alpha(i,j)+Fock_matrix_mo_beta(i,j))&
   !          + 0.5d0*(Fock_matrix_mo_beta(i,j) - Fock_matrix_mo_alpha(i,j))
   !    enddo
   !    ! F
   !    do i=elec_beta_num+1,elec_alpha_num !AA
   !      Fock_matrix_mo(i,j) = 0.5d0*(Fock_matrix_mo_alpha(i,j)+Fock_matrix_mo_beta(i,j))
   !    enddo
   !    ! F-K/2
   !    do i=elec_alpha_num+1, mo_num !AV
   !      Fock_matrix_mo(i,j) = 0.5d0*(Fock_matrix_mo_alpha(i,j)+Fock_matrix_mo_beta(i,j))&
   !          - 0.5d0*(Fock_matrix_mo_beta(i,j) - Fock_matrix_mo_alpha(i,j))
   !    enddo
   !  enddo

   !  do j=elec_alpha_num+1, mo_num
   !    ! F
   !    do i=1,elec_beta_num !VC
   !      Fock_matrix_mo(i,j) = 0.5d0*(Fock_matrix_mo_alpha(i,j)+Fock_matrix_mo_beta(i,j))
   !    enddo
   !    ! F-K/2
   !    do i=elec_beta_num+1,elec_alpha_num !VA
   !      Fock_matrix_mo(i,j) = 0.5d0*(Fock_matrix_mo_alpha(i,j)+Fock_matrix_mo_beta(i,j))&
   !          - 0.5d0*(Fock_matrix_mo_beta(i,j) - Fock_matrix_mo_alpha(i,j))
   !    enddo
   !    ! F+K
   !    do i=elec_alpha_num+1,mo_num !VV
   !      Fock_matrix_mo(i,j) = 0.5d0*(Fock_matrix_mo_alpha(i,j)+Fock_matrix_mo_beta(i,j)) &
   !          + (Fock_matrix_mo_beta(i,j) - Fock_matrix_mo_alpha(i,j))
   !    enddo
   !  enddo

   !endif

   do i = 1, mo_num
     Fock_matrix_diag_mo(i) = Fock_matrix_mo(i,i)
   enddo


   if(frozen_orb_scf)then
     integer                        :: iorb,jorb
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
   if ( (elec_alpha_num == elec_beta_num).and.                       &
         (level_shift == 0.) )                                       &
         then
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
 SCF_energy = nuclear_repulsion

 integer                        :: i,j
 do j=1,ao_num
   do i=1,ao_num
     SCF_energy += 0.5d0 * (                                          &
         (ao_one_e_integrals(i,j) + Fock_matrix_ao_alpha(i,j) ) *  SCF_density_matrix_ao_alpha(i,j) +&
         (ao_one_e_integrals(i,j) + Fock_matrix_ao_beta (i,j) ) *  SCF_density_matrix_ao_beta (i,j) )
   enddo
 enddo
 SCF_energy += extra_e_contrib_density

END_PROVIDER

! ---

