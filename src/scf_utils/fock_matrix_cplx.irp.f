 BEGIN_PROVIDER [ complex*16, Fock_matrix_mo_complex, (mo_num,mo_num) ]
&BEGIN_PROVIDER [ double precision, Fock_matrix_diag_mo_complex, (mo_num)]
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
     Fock_matrix_mo_complex = Fock_matrix_mo_alpha_complex
   else

     do j=1,elec_beta_num
       ! F-K
       do i=1,elec_beta_num !CC
         Fock_matrix_mo_complex(i,j) = 0.5d0*(Fock_matrix_mo_alpha_complex(i,j)+Fock_matrix_mo_beta_complex(i,j))&
             - (Fock_matrix_mo_beta_complex(i,j) - Fock_matrix_mo_alpha_complex(i,j))
       enddo
       ! F+K/2
       do i=elec_beta_num+1,elec_alpha_num  !CA
         Fock_matrix_mo_complex(i,j) = 0.5d0*(Fock_matrix_mo_alpha_complex(i,j)+Fock_matrix_mo_beta_complex(i,j))&
             + 0.5d0*(Fock_matrix_mo_beta_complex(i,j) - Fock_matrix_mo_alpha_complex(i,j))
       enddo
       ! F
       do i=elec_alpha_num+1, mo_num !CV
         Fock_matrix_mo_complex(i,j) = 0.5d0*(Fock_matrix_mo_alpha_complex(i,j)+Fock_matrix_mo_beta_complex(i,j))
       enddo
     enddo

     do j=elec_beta_num+1,elec_alpha_num
       ! F+K/2
       do i=1,elec_beta_num !AC
         Fock_matrix_mo_complex(i,j) = 0.5d0*(Fock_matrix_mo_alpha_complex(i,j)+Fock_matrix_mo_beta_complex(i,j))&
             + 0.5d0*(Fock_matrix_mo_beta_complex(i,j) - Fock_matrix_mo_alpha_complex(i,j))
       enddo
       ! F
       do i=elec_beta_num+1,elec_alpha_num !AA
         Fock_matrix_mo_complex(i,j) = 0.5d0*(Fock_matrix_mo_alpha_complex(i,j)+Fock_matrix_mo_beta_complex(i,j))
       enddo
       ! F-K/2
       do i=elec_alpha_num+1, mo_num !AV
         Fock_matrix_mo_complex(i,j) = 0.5d0*(Fock_matrix_mo_alpha_complex(i,j)+Fock_matrix_mo_beta_complex(i,j))&
             - 0.5d0*(Fock_matrix_mo_beta_complex(i,j) - Fock_matrix_mo_alpha_complex(i,j))
       enddo
     enddo

     do j=elec_alpha_num+1, mo_num
       ! F
       do i=1,elec_beta_num !VC
         Fock_matrix_mo_complex(i,j) = 0.5d0*(Fock_matrix_mo_alpha_complex(i,j)+Fock_matrix_mo_beta_complex(i,j))
       enddo
       ! F-K/2
       do i=elec_beta_num+1,elec_alpha_num !VA
         Fock_matrix_mo_complex(i,j) = 0.5d0*(Fock_matrix_mo_alpha_complex(i,j)+Fock_matrix_mo_beta_complex(i,j))&
             - 0.5d0*(Fock_matrix_mo_beta_complex(i,j) - Fock_matrix_mo_alpha_complex(i,j))
       enddo
       ! F+K
       do i=elec_alpha_num+1,mo_num !VV
         Fock_matrix_mo_complex(i,j) = 0.5d0*(Fock_matrix_mo_alpha_complex(i,j)+Fock_matrix_mo_beta_complex(i,j)) &
             + (Fock_matrix_mo_beta_complex(i,j) - Fock_matrix_mo_alpha_complex(i,j))
       enddo
     enddo

   endif

   do i = 1, mo_num
     Fock_matrix_diag_mo_complex(i) = dble(Fock_matrix_mo_complex(i,i))
     if (dabs(dimag(Fock_matrix_mo_complex(i,i))) .gt. 1.0d-12) then
       !stop 'diagonal elements of Fock matrix should be real'
       print *, 'diagonal elements of Fock matrix should be real',i,Fock_matrix_mo_complex(i,i)
       stop -1
     endif
   enddo


   if(frozen_orb_scf)then
     integer                        :: iorb,jorb
     do i = 1, n_core_orb
      iorb = list_core(i)
      do j = 1, n_act_orb
       jorb = list_act(j)
       Fock_matrix_mo_complex(iorb,jorb) = (0.d0,0.d0)
       Fock_matrix_mo_complex(jorb,iorb) = (0.d0,0.d0)
      enddo
     enddo
   endif

END_PROVIDER



BEGIN_PROVIDER [ complex*16, Fock_matrix_mo_alpha_complex, (mo_num,mo_num) ]
   implicit none
   BEGIN_DOC
   ! Fock matrix on the MO basis
   END_DOC
   call ao_to_mo_complex(Fock_matrix_ao_alpha_complex,size(Fock_matrix_ao_alpha_complex,1), &
                 Fock_matrix_mo_alpha_complex,size(Fock_matrix_mo_alpha_complex,1))
END_PROVIDER

BEGIN_PROVIDER [ complex*16, Fock_matrix_mo_beta_complex, (mo_num,mo_num) ]
   implicit none
   BEGIN_DOC
   ! Fock matrix on the MO basis
   END_DOC
   call ao_to_mo_complex(Fock_matrix_ao_beta_complex,size(Fock_matrix_ao_beta_complex,1), &
                 Fock_matrix_mo_beta_complex,size(Fock_matrix_mo_beta_complex,1))
END_PROVIDER


BEGIN_PROVIDER [ complex*16, Fock_matrix_ao_complex, (ao_num, ao_num) ]
 implicit none
 BEGIN_DOC
 ! Fock matrix in AO basis set
 END_DOC

 if(frozen_orb_scf)then
   call mo_to_ao_complex(Fock_matrix_mo_complex,size(Fock_matrix_mo_complex,1),              &
       Fock_matrix_ao_complex,size(Fock_matrix_ao_complex,1))
 else
   if ( (elec_alpha_num == elec_beta_num).and.                       &
         (level_shift == 0.) )                                       &
         then
     integer                        :: i,j
     do j=1,ao_num
       do i=1,ao_num
         Fock_matrix_ao_complex(i,j) = Fock_matrix_ao_alpha_complex(i,j)
       enddo
     enddo
   else
     call mo_to_ao_complex(Fock_matrix_mo_complex,size(Fock_matrix_mo_complex,1),            &
         Fock_matrix_ao_complex,size(Fock_matrix_ao_complex,1))
   endif
 endif
END_PROVIDER

