       subroutine driver_optorb
         implicit none
         integer :: i,j

         write(6,*)
!        write(6,*) ' <0|H|0>     (qp)      = ',psi_energy_with_nucl_rep(1)
         write(6,*) ' energy improvement    = ',energy_improvement
!        write(6,*) ' new energy            = ',psi_energy_with_nucl_rep(1)+energy_improvement
         write(6,*)

         write(6,*)
         write(6,*) '   creating new orbitals '
         do i=1,mo_num
          write(6,*) ' Orbital No ',i
          write(6,'(5F14.6)') (NewOrbs(j,i),j=1,mo_num)
          write(6,*)
         end do

         mo_label = "Natural"
         do i=1,mo_num
          do j=1,ao_num
           mo_coef(j,i)=NewOrbs(j,i)
          end do
         end do
         call save_mos
         call map_deinit(mo_integrals_map)
         FREE mo_integrals_map mo_coef mo_two_e_integrals_in_map
 
         write(6,*)
         write(6,*) '   ... all done '

       end 
