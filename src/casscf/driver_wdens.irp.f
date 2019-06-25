       subroutine driver_wdens
         implicit none
         integer :: istate,p,q,r,s,indx,i,j


         write(6,*) ' total energy = ',eone+etwo+ecore
         write(6,*) ' generating natural orbitals '
         write(6,*)
         write(6,*)

         write(6,*) ' all data available ! '
         write(6,*) '    writing out files '

         call trf_to_natorb
real*8 :: approx,np,nq,nr,ns
logical :: lpq,lrs,lps,lqr
         
         open(unit=12,form='formatted',status='unknown',file='onetrf.tmp')
         indx=0
         do q=1,mo_num
          do p=q,mo_num
            if (abs(onetrf(p,q)).gt.1.D-12) then
             write(12,'(2i6,E20.12)') p,q,onetrf(p,q)
             indx+=1
            end if
          end do 
         end do 
         write(6,*) ' wrote ',indx,' mono-electronic integrals'
         close(12)
 
         
         write(6,*)
         write(6,*) '   creating new orbitals '
         do i=1,mo_num
          write(6,*) ' Orbital No ',i
          write(6,'(5F14.6)') (NatOrbsFCI(j,i),j=1,mo_num)
          write(6,*)
         end do

         mo_label = "MCSCF"
         mo_label = "Natural"
         do i=1,mo_num
          do j=1,ao_num
           mo_coef(j,i)=NatOrbsFCI(j,i)
          end do
         end do
         call save_mos
 
         write(6,*) '   ... done '

       end

