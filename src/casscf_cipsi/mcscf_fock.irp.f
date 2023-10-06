BEGIN_PROVIDER [real*8, Fipq, (mo_num,mo_num) ]
   BEGIN_DOC
   ! the inactive Fock matrix, in molecular orbitals
   END_DOC
   implicit none
   integer                        :: p,q,k,kk,t,tt,u,uu
   
   do q=1,mo_num
     do p=1,mo_num
       Fipq(p,q)=one_ints_no(p,q)
     end do
   end do
   
   ! the inactive Fock matrix
   do k=1,n_core_inact_orb
     kk=list_core_inact(k)
     do q=1,mo_num
       do p=1,mo_num
         Fipq(p,q)+=2.D0*bielec_pqxx_no(p,q,k,k) -bielec_pxxq_no(p,k,k,q)
       end do
     end do
   end do
   
   if (bavard) then
     integer                        :: i
     write(6,*)
     write(6,*) ' the diagonal of the inactive effective Fock matrix '
     write(6,'(5(i3,F12.5))') (i,Fipq(i,i),i=1,mo_num)
     write(6,*)
   end if
   
   
END_PROVIDER
 
 
BEGIN_PROVIDER [real*8, Fapq, (mo_num,mo_num) ]
   BEGIN_DOC
   ! the active active Fock matrix, in molecular orbitals
   ! we create them in MOs, quite expensive
   !
   ! for an implementation in AOs we need first the natural orbitals
   ! for forming an active density matrix in AOs
   !
   END_DOC
   implicit none
   integer                        :: p,q,k,kk,t,tt,u,uu
   
   Fapq = 0.d0
   
   ! the active Fock matrix, D0tu is diagonal
   do t=1,n_act_orb
     tt=list_act(t)
     do q=1,mo_num
       do p=1,mo_num
         Fapq(p,q)+=occnum(tt)                                       &
             *(bielec_pqxx_no(p,q,tt,tt)-0.5D0*bielec_pxxq_no(p,tt,tt,q))
       end do
     end do
   end do
   
   if (bavard) then
     integer                        :: i
     write(6,*)
     write(6,*) ' the effective Fock matrix over MOs'
     write(6,*)
     
     write(6,*)
     write(6,*) ' the diagonal of the inactive effective Fock matrix '
     write(6,'(5(i3,F12.5))') (i,Fipq(i,i),i=1,mo_num)
     write(6,*)
     write(6,*)
     write(6,*) ' the diagonal of the active Fock matrix '
     write(6,'(5(i3,F12.5))') (i,Fapq(i,i),i=1,mo_num)
     write(6,*)
   end if
   
   
END_PROVIDER
 
 BEGIN_PROVIDER [ double precision, mcscf_fock_alpha_ao, (ao_num, ao_num)] 
&BEGIN_PROVIDER [ double precision, mcscf_fock_beta_ao, (ao_num, ao_num)] 
 implicit none
 BEGIN_DOC
  ! mcscf_fock_alpha_ao are set to usual Fock like operator but computed with the MCSCF densities on the AO basis 
 END_DOC
 SCF_density_matrix_ao_alpha = D0tu_alpha_ao
 SCF_density_matrix_ao_beta = D0tu_beta_ao
 soft_touch SCF_density_matrix_ao_alpha SCF_density_matrix_ao_beta 
 mcscf_fock_beta_ao = fock_matrix_ao_beta
 mcscf_fock_alpha_ao = fock_matrix_ao_alpha
END_PROVIDER 


 BEGIN_PROVIDER [ double precision, mcscf_fock_alpha_mo, (mo_num, mo_num)] 
&BEGIN_PROVIDER [ double precision, mcscf_fock_beta_mo, (mo_num, mo_num)] 
 implicit none
 BEGIN_DOC
  ! Mo_mcscf_fock_alpha are set to usual Fock like operator but computed with the MCSCF densities on the MO basis 
 END_DOC

 call ao_to_mo(mcscf_fock_alpha_ao,ao_num,mcscf_fock_alpha_mo,mo_num)
 call ao_to_mo(mcscf_fock_beta_ao,ao_num,mcscf_fock_beta_mo,mo_num)

END_PROVIDER 

 BEGIN_PROVIDER [ double precision, mcscf_fock_mo, (mo_num,mo_num) ]
&BEGIN_PROVIDER [ double precision, mcscf_fock_diag_mo, (mo_num)]
   implicit none
   BEGIN_DOC
   ! MCSF Fock matrix on the MO basis.
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
     mcscf_fock_mo = mcscf_fock_alpha_mo
   else
     ! Core
     do j = 1, elec_beta_num
       ! Core
       do i = 1, elec_beta_num
         mcscf_fock_mo(i,j) = - 0.5d0 * mcscf_fock_alpha_mo(i,j) &
                               + 1.5d0 * mcscf_fock_beta_mo(i,j)
       enddo
       ! Open
       do i = elec_beta_num+1, elec_alpha_num
         mcscf_fock_mo(i,j) = mcscf_fock_beta_mo(i,j)
       enddo
       ! Virtual
       do i = elec_alpha_num+1, mo_num
         mcscf_fock_mo(i,j) =   0.5d0 * mcscf_fock_alpha_mo(i,j) &
                               + 0.5d0 * mcscf_fock_beta_mo(i,j)
       enddo
     enddo
     ! Open
     do j = elec_beta_num+1, elec_alpha_num
       ! Core
       do i = 1, elec_beta_num
         mcscf_fock_mo(i,j) = mcscf_fock_beta_mo(i,j)
       enddo
       ! Open
       do i = elec_beta_num+1, elec_alpha_num
         mcscf_fock_mo(i,j) =   0.5d0 * mcscf_fock_alpha_mo(i,j) &
                               + 0.5d0 * mcscf_fock_beta_mo(i,j)
       enddo
       ! Virtual
       do i = elec_alpha_num+1, mo_num
         mcscf_fock_mo(i,j) = mcscf_fock_alpha_mo(i,j)
       enddo
     enddo
     ! Virtual
     do j = elec_alpha_num+1, mo_num
       ! Core
       do i = 1, elec_beta_num
         mcscf_fock_mo(i,j) =   0.5d0 * mcscf_fock_alpha_mo(i,j) &
                               + 0.5d0 * mcscf_fock_beta_mo(i,j)
       enddo
       ! Open
       do i = elec_beta_num+1, elec_alpha_num
         mcscf_fock_mo(i,j) = mcscf_fock_alpha_mo(i,j)
       enddo
       ! Virtual
       do i = elec_alpha_num+1, mo_num
         mcscf_fock_mo(i,j) =   1.5d0 * mcscf_fock_alpha_mo(i,j) &
                               - 0.5d0 * mcscf_fock_beta_mo(i,j)
       enddo
     enddo
   endif

 do i = 1, mo_num
  mcscf_fock_diag_mo(i) = mcscf_fock_mo(i,i)
 enddo
END_PROVIDER 
