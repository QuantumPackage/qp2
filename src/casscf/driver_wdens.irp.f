       subroutine driver_wdens
         implicit none
         integer :: istate,p,q,r,s,indx,i,j


         write(6,*) ' total energy = ',eone+etwo+ecore
         write(6,*) ' generating natural orbitals '
         write(6,*)
         write(6,*)
         call trf_to_natorb

         write(6,*) ' all data available ! '
         write(6,*) '    writing out files '

         open(unit=12,file='D0tu.dat',form='formatted',status='unknown')
         do p=1,n_act_orb
          do q=1,n_act_orb
           if (abs(D0tu(p,q)).gt.1.D-12) then
            write(12,'(2i8,E20.12)') p,q,D0tu(p,q)
           end if
          end do
         end do 
         close(12)
         
real*8 :: approx,np,nq,nr,ns
logical :: lpq,lrs,lps,lqr
         open(unit=12,file='P0tuvx.dat',form='formatted',status='unknown')
         do p=1,n_act_orb
          np=D0tu(p,p)
          do q=1,n_act_orb
           lpq=p.eq.q
           nq=D0tu(q,q)
           do r=1,n_act_orb
            lqr=q.eq.r
            nr=D0tu(r,r)
            do s=1,n_act_orb
             lrs=r.eq.s
             lps=p.eq.s
             approx=0.D0
             if (lpq.and.lrs) then
              if (lqr) then
! pppp
               approx=0.5D0*np*(np-1.D0)
              else
! pprr
               approx=0.5D0*np*nr
              end if
             else
              if (lps.and.lqr.and..not.lpq) then
! pqqp
               approx=-0.25D0*np*nq
              end if
             end if
             if (abs(P0tuvx(p,q,r,s)).gt.1.D-12) then
              write(12,'(4i4,2E20.12)') p,q,r,s,P0tuvx(p,q,r,s),approx
             end if
            end do
           end do
          end do 
         end do 
         close(12)
         
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
 
         
         open(unit=12,form='formatted',status='unknown',file='bielec_PQxx.tmp')
         indx=0
         do p=1,mo_num
          do q=p,mo_num
           do r=1,n_core_orb+n_act_orb
            do s=r,n_core_orb+n_act_orb
             if (abs(bielec_PQxxtmp(p,q,r,s)).gt.1.D-12) then
              write(12,'(4i8,E20.12)') p,q,r,s,bielec_PQxxtmp(p,q,r,s)
              indx+=1
             end if
            end do
           end do
          end do 
         end do 
         write(6,*) ' wrote ',indx,' integrals (PQ|xx)'
         close(12)
 
         open(unit=12,form='formatted',status='unknown',file='bielec_PxxQ.tmp')
         indx=0
         do p=1,mo_num
          do q=1,n_core_orb+n_act_orb
           do r=q,n_core_orb+n_act_orb
integer ::s_start
            if (q.eq.r) then
             s_start=p
            else
             s_start=1
            end if
            do s=s_start,mo_num
             if (abs(bielec_PxxQtmp(p,q,r,s)).gt.1.D-12) then
              write(12,'(4i8,E20.12)') p,q,r,s,bielec_PxxQtmp(p,q,r,s)
              indx+=1
             end if
            end do
           end do
          end do 
         end do 
         write(6,*) ' wrote ',indx,' integrals (Px|xQ)'
         close(12)
 
         open(unit=12,form='formatted',status='unknown',file='bielecCI.tmp')
         indx=0
         do p=1,n_act_orb
          do q=p,n_act_orb
           do r=1,n_act_orb
            do s=1,mo_num
             if (abs(bielecCItmp(p,q,r,s)).gt.1.D-12) then
              write(12,'(4i8,E20.12)') p,q,r,s,bielecCItmp(p,q,r,s)
              indx+=1
             end if
            end do
           end do
          end do 
         end do 
         write(6,*) ' wrote ',indx,' integrals (tu|xP)'
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

