 BEGIN_PROVIDER [ double precision, Fock_matrix_mo, (mo_num,mo_num) ]
&BEGIN_PROVIDER [ double precision, Fock_matrix_diag_mo, (mo_num)]
   implicit none
   BEGIN_DOC
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
   END_DOC
   integer                        :: i,j,n
   if (elec_alpha_num == elec_beta_num) then
     Fock_matrix_mo = Fock_matrix_mo_alpha
   else

     do j=1,elec_beta_num
       ! F-K
       do i=1,elec_beta_num
         Fock_matrix_mo(i,j) = 0.5d0*(Fock_matrix_mo_alpha(i,j)+Fock_matrix_mo_beta(i,j))&
             - (Fock_matrix_mo_beta(i,j) - Fock_matrix_mo_alpha(i,j))
       enddo
       ! F+K/2
       do i=elec_beta_num+1,elec_alpha_num
         Fock_matrix_mo(i,j) = 0.5d0*(Fock_matrix_mo_alpha(i,j)+Fock_matrix_mo_beta(i,j))&
             + 0.5d0*(Fock_matrix_mo_beta(i,j) - Fock_matrix_mo_alpha(i,j))
       enddo
       ! F
       do i=elec_alpha_num+1, mo_num
         Fock_matrix_mo(i,j) = 0.5d0*(Fock_matrix_mo_alpha(i,j)+Fock_matrix_mo_beta(i,j))
       enddo
     enddo

     do j=elec_beta_num+1,elec_alpha_num
       ! F+K/2
       do i=1,elec_beta_num
         Fock_matrix_mo(i,j) = 0.5d0*(Fock_matrix_mo_alpha(i,j)+Fock_matrix_mo_beta(i,j))&
             + 0.5d0*(Fock_matrix_mo_beta(i,j) - Fock_matrix_mo_alpha(i,j))
       enddo
       ! F
       do i=elec_beta_num+1,elec_alpha_num
         Fock_matrix_mo(i,j) = 0.5d0*(Fock_matrix_mo_alpha(i,j)+Fock_matrix_mo_beta(i,j))
       enddo
       ! F-K/2
       do i=elec_alpha_num+1, mo_num
         Fock_matrix_mo(i,j) = 0.5d0*(Fock_matrix_mo_alpha(i,j)+Fock_matrix_mo_beta(i,j))&
             - 0.5d0*(Fock_matrix_mo_beta(i,j) - Fock_matrix_mo_alpha(i,j))
       enddo
     enddo

     do j=elec_alpha_num+1, mo_num
       ! F
       do i=1,elec_beta_num
         Fock_matrix_mo(i,j) = 0.5d0*(Fock_matrix_mo_alpha(i,j)+Fock_matrix_mo_beta(i,j))
       enddo
       ! F-K/2
       do i=elec_beta_num+1,elec_alpha_num
         Fock_matrix_mo(i,j) = 0.5d0*(Fock_matrix_mo_alpha(i,j)+Fock_matrix_mo_beta(i,j))&
             - 0.5d0*(Fock_matrix_mo_beta(i,j) - Fock_matrix_mo_alpha(i,j))
       enddo
       ! F+K
       do i=elec_alpha_num+1,mo_num
         Fock_matrix_mo(i,j) = 0.5d0*(Fock_matrix_mo_alpha(i,j)+Fock_matrix_mo_beta(i,j)) &
             + (Fock_matrix_mo_beta(i,j) - Fock_matrix_mo_alpha(i,j))
       enddo
     enddo

   endif

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
   call mo_to_ao(Fock_matrix_mo,size(Fock_matrix_mo,1), &
      Fock_matrix_ao,size(Fock_matrix_ao,1))
 else
  if ( (elec_alpha_num == elec_beta_num).and. &
       (level_shift == 0.) ) &
  then
    integer                        :: i,j
    do j=1,ao_num
      do i=1,ao_num
        Fock_matrix_ao(i,j) = Fock_matrix_ao_alpha(i,j)
      enddo
    enddo
  else
    call mo_to_ao(Fock_matrix_mo,size(Fock_matrix_mo,1), &
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

