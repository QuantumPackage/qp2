! -*- F90 -*- 
 BEGIN_PROVIDER[real*8, bielec_PQxx, (mo_num, mo_num,n_core_orb+n_act_orb,n_core_orb+n_act_orb)]
&BEGIN_PROVIDER[real*8, bielec_PxxQ, (mo_num,n_core_orb+n_act_orb,n_core_orb+n_act_orb, mo_num)]
BEGIN_DOC
! bielec_PQxx : integral (pq|xx) with p,q arbitrary, x core or active
! bielec_PxxQ : integral (px|xq) with p,q arbitrary, x core or active
! indices are unshifted orbital numbers
! all integrals are read from files
END_DOC
        implicit none
          integer :: i,j,p,q,indx,kk
          real*8 :: hhh
          logical :: lread

          do i=1,n_core_orb+n_act_orb
           do j=1,n_core_orb+n_act_orb
            do p=1,mo_num
             do q=1,mo_num
              bielec_PQxx(p,q,i,j)=0.D0
              bielec_PxxQ(p,i,j,q)=0.D0
             end do
            end do
           end do
          end do

          open(unit=12,form='formatted',status='old',file='bielec_PQxx.tmp')
          lread=.true.
          indx=0
          do while (lread)
           read(12,*,iostat=kk) p,q,i,j,hhh
           if (kk.ne.0) then
            lread=.false.
           else
! stored with p.le.q, i.le.j
            bielec_PQxx(p,q,i,j)=hhh
            bielec_PQxx(q,p,i,j)=hhh
            bielec_PQxx(q,p,j,i)=hhh
            bielec_PQxx(p,q,j,i)=hhh
            indx+=1
           end if
          end do
          close(12)
          write(6,*) ' read ',indx,' integrals PQxx into core '
            
          open(unit=12,form='formatted',status='old',file='bielec_PxxQ.tmp')
          lread=.true.
          indx=0
          do while (lread)
           read(12,*,iostat=kk) p,i,j,q,hhh
           if (kk.ne.0) then
            lread=.false.
           else
! stored with (ip).le.(jq)
            bielec_PxxQ(p,i,j,q)=hhh
            bielec_PxxQ(q,j,i,p)=hhh
            indx+=1
           end if
          end do
          write(6,*) ' read ',indx,' integrals PxxQ into core '
          close(12)
          write(6,*) ' provided integrals (PQ|xx) and (Px|xQ) '
END_PROVIDER

BEGIN_PROVIDER[real*8, bielecCI, (n_act_orb,n_act_orb,n_act_orb, mo_num)]
BEGIN_DOC
! bielecCI : integrals (tu|vp) with p arbitrary, tuv active
! index p runs over the whole basis, t,u,v only over the active orbitals
! inegrals read from file
END_DOC
        implicit none
          integer :: i,j,k,p,t,u,v,kk,indx
          real*8 :: hhh
          logical :: lread

          write(6,*) ' reading integrals bielecCI '

          do i=1,n_act_orb
           do j=1,n_act_orb
            do k=1,n_act_orb
             do p=1,mo_num
              bielecCI(i,k,j,p)=0.D0
             end do
            end do
           end do
          end do

          open(unit=12,form='formatted',status='old',file='bielecCI.tmp')
          lread=.true.
          indx=0
          do while (lread)
           read(12,*,iostat=kk) i,j,k,p,hhh
           if (kk.ne.0) then
            lread=.false.
           else
            bielecCI(i,j,k,p)=hhh
            bielecCI(j,i,k,p)=hhh
            indx+=1
           end if
          end do
          write(6,*) ' read ',indx,' integrals (aa|aP) into core '
          close(12)
          write(6,*) ' provided integrals (tu|xP) '
END_PROVIDER

