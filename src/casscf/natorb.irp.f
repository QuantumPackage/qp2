! -*- F90 -*-
! diagonalize D0tu
! save the diagonal somewhere, in inverse order
! 4-index-transform the 2-particle density matrix over active orbitals
! correct the bielectronic integrals
! correct the monoelectronic integrals
! put integrals on file, as well orbitals, and the density matrices
!  
     subroutine trf_to_natorb
      implicit none
      integer :: i,j,k,l,t,u,p,q,pp
      real*8 :: eigval(n_act_orb),natorbsCI(n_act_orb,n_act_orb)
      real*8 :: d(n_act_orb),d1(n_act_orb),d2(n_act_orb)

      call lapack_diag(eigval,natorbsCI,D0tu,n_act_orb,n_act_orb)
      write(6,*) ' found occupation numbers as '
      do i=1,n_act_orb
       write(6,*) i,eigval(i)
      end do

      if (bavard) then
!
       
integer :: nmx
real*8 :: xmx
       do i=1,n_act_orb
! largest element of the eigenvector should be positive
        xmx=0.D0
        nmx=0
        do j=1,n_act_orb
         if (abs(natOrbsCI(j,i)).gt.xmx) then
          nmx=j
          xmx=abs(natOrbsCI(j,i))
         end if
        end do
        xmx=sign(1.D0,natOrbsCI(nmx,i))
        do j=1,n_act_orb
         natOrbsCI(j,i)*=xmx
        end do
        
         
        write(6,*) ' Eigenvector No ',i
        write(6,'(5(I3,F12.5))') (j,natOrbsCI(j,i),j=1,n_act_orb)
       end do
      end if

      do i=1,n_act_orb
       do j=1,n_act_orb
        D0tu(i,j)=0.D0
       end do
! fill occupation numbers in descending order
       D0tu(i,i)=eigval(n_act_orb-i+1)
      end do
!
! 4-index transformation of 2part matrices
!
! index per index
! first quarter
       do j=1,n_act_orb
        do k=1,n_act_orb
         do l=1,n_act_orb
          do p=1,n_act_orb
           d(p)=0.D0
          end do
          do p=1,n_act_orb
           pp=n_act_orb-p+1
           do q=1,n_act_orb
            d(pp)+=P0tuvx(q,j,k,l)*natorbsCI(q,p)
           end do
          end do
          do p=1,n_act_orb
           P0tuvx(p,j,k,l)=d(p)
          end do
         end do
        end do
       end do
! 2nd quarter
       do j=1,n_act_orb
        do k=1,n_act_orb
         do l=1,n_act_orb
          do p=1,n_act_orb
           d(p)=0.D0
          end do
          do p=1,n_act_orb
           pp=n_act_orb-p+1
           do q=1,n_act_orb
            d(pp)+=P0tuvx(j,q,k,l)*natorbsCI(q,p)
           end do
          end do
          do p=1,n_act_orb
           P0tuvx(j,p,k,l)=d(p)
          end do
         end do
        end do
       end do
! 3rd quarter
       do j=1,n_act_orb
        do k=1,n_act_orb
         do l=1,n_act_orb
          do p=1,n_act_orb
           d(p)=0.D0
          end do
          do p=1,n_act_orb
           pp=n_act_orb-p+1
           do q=1,n_act_orb
            d(pp)+=P0tuvx(j,k,q,l)*natorbsCI(q,p)
           end do
          end do
          do p=1,n_act_orb
           P0tuvx(j,k,p,l)=d(p)
          end do
         end do
        end do
       end do
! 4th quarter
       do j=1,n_act_orb
        do k=1,n_act_orb
         do l=1,n_act_orb
          do p=1,n_act_orb
           d(p)=0.D0
          end do
          do p=1,n_act_orb
           pp=n_act_orb-p+1
           do q=1,n_act_orb
            d(pp)+=P0tuvx(j,k,l,q)*natorbsCI(q,p)
           end do
          end do
          do p=1,n_act_orb
           P0tuvx(j,k,l,p)=d(p)
          end do
         end do
        end do
       end do
       write(6,*) ' transformed P0tuvx '
!
! one-electron integrals 
!
       do i=1,mo_num
        do j=1,mo_num
         onetrf(i,j)=mo_one_e_integrals(i,j)
        end do
       end do
! 1st half-trf
       do j=1,mo_num
        do p=1,n_act_orb
         d(p)=0.D0
        end do
        do p=1,n_act_orb
         pp=n_act_orb-p+1
         do q=1,n_act_orb
          d(pp)+=onetrf(list_act(q),j)*natorbsCI(q,p)
         end do
        end do
         do p=1,n_act_orb
          onetrf(list_act(p),j)=d(p)
         end do
        end do
! 2nd half-trf
       do j=1,mo_num
        do p=1,n_act_orb
         d(p)=0.D0
        end do
        do p=1,n_act_orb
         pp=n_act_orb-p+1
         do q=1,n_act_orb
          d(pp)+=onetrf(j,list_act(q))*natorbsCI(q,p)
         end do
        end do
        do p=1,n_act_orb
         onetrf(j,list_act(p))=d(p)
        end do
       end do
       write(6,*) ' transformed onetrf '
! 
! Orbitals
!
       do j=1,ao_num
        do i=1,mo_num
         NatOrbsFCI(j,i)=mo_coef(j,i)
        end do
       end do

       do j=1,ao_num
        do p=1,n_act_orb
         d(p)=0.D0
        end do
        do p=1,n_act_orb
         pp=n_act_orb-p+1
         do q=1,n_act_orb
          d(pp)+=NatOrbsFCI(j,list_act(q))*natorbsCI(q,p)
         end do
        end do
        do p=1,n_act_orb
         NatOrbsFCI(j,list_act(p))=d(p)
        end do
       end do
       write(6,*) ' transformed orbitals '
!
! now the bielectronic integrals
!
!!$       write(6,*) ' before the transformation '
!!$integer :: kk,ll,ii,jj
!!$real*8 :: h1,h2,h3
!!$       do i=1,n_act_orb
!!$          ii=list_act(i)
!!$          do j=1,n_act_orb
!!$             jj=list_act(j)
!!$             do k=1,n_act_orb
!!$                kk=list_act(k)
!!$                do l=1,n_act_orb
!!$                   ll=list_act(l)
!!$                   h1=bielec_PQxxtmp(ii,jj,k+n_core_orb,l+n_core_orb)
!!$                   h2=bielec_PxxQtmp(ii,j+n_core_orb,k+n_core_orb,ll)
!!$                   h3=bielecCItmp(i,j,k,ll)
!!$                   if ((h1.ne.h2).or.(h1.ne.h3)) then
!!$                      write(6,9901) i,j,k,l,h1,h2,h3
!!$9901                  format(' aie   ',4i4,3E20.12)
!!$9902                  format('correct',4i4,3E20.12)
!!$                   else
!!$                      write(6,9902) i,j,k,l,h1,h2,h3
!!$                   end if
!!$                end do
!!$             end do
!!$          end do
!!$       end do

       do j=1,mo_num
        do k=1,n_core_orb+n_act_orb
         do l=1,n_core_orb+n_act_orb
          do p=1,n_act_orb
           d1(p)=0.D0
           d2(p)=0.D0
          end do
          do p=1,n_act_orb
           pp=n_act_orb-p+1
           do q=1,n_act_orb
            d1(pp)+=bielec_PQxxtmp(list_act(q),j,k,l)*natorbsCI(q,p)
            d2(pp)+=bielec_PxxQtmp(list_act(q),k,l,j)*natorbsCI(q,p)
           end do
          end do
          do p=1,n_act_orb
           bielec_PQxxtmp(list_act(p),j,k,l)=d1(p)
           bielec_PxxQtmp(list_act(p),k,l,j)=d2(p)
          end do
         end do
        end do
       end do
! 2nd quarter
       do j=1,mo_num
        do k=1,n_core_orb+n_act_orb
         do l=1,n_core_orb+n_act_orb
          do p=1,n_act_orb
           d1(p)=0.D0
           d2(p)=0.D0
          end do
          do p=1,n_act_orb
           pp=n_act_orb-p+1
           do q=1,n_act_orb
            d1(pp)+=bielec_PQxxtmp(j,list_act(q),k,l)*natorbsCI(q,p)
            d2(pp)+=bielec_PxxQtmp(j,k,l,list_act(q))*natorbsCI(q,p)
           end do
          end do
          do p=1,n_act_orb
           bielec_PQxxtmp(j,list_act(p),k,l)=d1(p)
           bielec_PxxQtmp(j,k,l,list_act(p))=d2(p)
          end do
         end do
        end do
       end do
! 3rd quarter
       do j=1,mo_num
        do k=1,mo_num
         do l=1,n_core_orb+n_act_orb
          do p=1,n_act_orb
           d1(p)=0.D0
           d2(p)=0.D0
          end do
          do p=1,n_act_orb
           pp=n_act_orb-p+1
           do q=1,n_act_orb
            d1(pp)+=bielec_PQxxtmp(j,k,n_core_orb+q,l)*natorbsCI(q,p)
            d2(pp)+=bielec_PxxQtmp(j,n_core_orb+q,l,k)*natorbsCI(q,p)
           end do
          end do
          do p=1,n_act_orb
           bielec_PQxxtmp(j,k,n_core_orb+p,l)=d1(p)
           bielec_PxxQtmp(j,n_core_orb+p,l,k)=d2(p)
          end do
         end do
        end do
       end do
! 4th quarter
       do j=1,mo_num
        do k=1,mo_num
         do l=1,n_core_orb+n_act_orb
          do p=1,n_act_orb
           d1(p)=0.D0
           d2(p)=0.D0
          end do
          do p=1,n_act_orb
           pp=n_act_orb-p+1
           do q=1,n_act_orb
            d1(pp)+=bielec_PQxxtmp(j,k,l,n_core_orb+q)*natorbsCI(q,p)
            d2(pp)+=bielec_PxxQtmp(j,l,n_core_orb+q,k)*natorbsCI(q,p)
           end do
          end do
          do p=1,n_act_orb
           bielec_PQxxtmp(j,k,l,n_core_orb+p)=d1(p)
           bielec_PxxQtmp(j,l,n_core_orb+p,k)=d2(p)
          end do
         end do
        end do
       end do
       write(6,*) ' transformed PQxx and PxxQ '
!
! and finally the bielecCI integrals
!
       do j=1,n_act_orb
        do k=1,n_act_orb
         do l=1,mo_num
          do p=1,n_act_orb
           d(p)=0.D0
          end do
          do p=1,n_act_orb
           pp=n_act_orb-p+1
           do q=1,n_act_orb
            d(pp)+=bielecCItmp(q,j,k,l)*natorbsCI(q,p)
           end do
          end do
          do p=1,n_act_orb
           bielecCItmp(p,j,k,l)=d(p)
          end do
         end do
        end do
       end do
! 2nd quarter
       do j=1,n_act_orb
        do k=1,n_act_orb
         do l=1,mo_num
          do p=1,n_act_orb
           d(p)=0.D0
          end do
          do p=1,n_act_orb
           pp=n_act_orb-p+1
           do q=1,n_act_orb
            d(pp)+=bielecCItmp(j,q,k,l)*natorbsCI(q,p)
           end do
          end do
          do p=1,n_act_orb
           bielecCItmp(j,p,k,l)=d(p)
          end do
         end do
        end do
       end do
! 3rd quarter
       do j=1,n_act_orb
        do k=1,n_act_orb
         do l=1,mo_num
          do p=1,n_act_orb
           d(p)=0.D0
          end do
          do p=1,n_act_orb
           pp=n_act_orb-p+1
           do q=1,n_act_orb
            d(pp)+=bielecCItmp(j,k,q,l)*natorbsCI(q,p)
           end do
          end do
          do p=1,n_act_orb
           bielecCItmp(j,k,p,l)=d(p)
          end do
         end do
        end do
       end do
! 4th quarter
       do j=1,n_act_orb
        do k=1,n_act_orb
         do l=1,n_act_orb
          do p=1,n_act_orb
           d(p)=0.D0
          end do
          do p=1,n_act_orb
           pp=n_act_orb-p+1
           do q=1,n_act_orb
            d(pp)+=bielecCItmp(j,k,l,list_act(q))*natorbsCI(q,p)
           end do
          end do
          do p=1,n_act_orb
           bielecCItmp(j,k,l,list_act(p))=d(p)
          end do
         end do
        end do
       end do
       write(6,*) ' transformed tuvP '
!
! that's all
!
!!$
!!$! test coherence of the bielectronic integals
!!$!   PQxx = PxxQ = tuvP for some of the indices
!!$       write(6,*) ' after the transformation '
!!$       do i=1,n_act_orb
!!$          ii=list_act(i)
!!$          do j=1,n_act_orb
!!$             jj=list_act(j)
!!$             do k=1,n_act_orb
!!$                kk=list_act(k)
!!$                do l=1,n_act_orb
!!$                   ll=list_act(l)
!!$                   h1=bielec_PQxxtmp(ii,jj,k+n_core_orb,l+n_core_orb)
!!$                   h2=bielec_PxxQtmp(ii,j+n_core_orb,k+n_core_orb,ll)
!!$                   h3=bielecCItmp(i,j,k,ll)
!!$                   if ((abs(h1-h2).gt.1.D-14).or.(abs(h1-h3).gt.1.D-14)) then
!!$                      write(6,9901) i,j,k,l,h1,h1-h2,h1-h3
!!$                   else
!!$                      write(6,9902) i,j,k,l,h1,h2,h3
!!$                   end if
!!$                end do
!!$             end do
!!$          end do
!!$       end do

! we recalculate total energies
       write(6,*) 
       write(6,*)  ' recalculating energies after the transformation '
       write(6,*) 
       write(6,*) 
       real*8 :: e_one_all
       real*8 :: e_two_all
       integer :: ii
       integer :: jj
       integer :: t3
       integer :: tt
       integer :: u3
       integer :: uu
       integer :: v
       integer :: v3
       integer :: vv
       integer :: x
       integer :: x3
       integer :: xx

       e_one_all=0.D0
       e_two_all=0.D0
       do i=1,n_core_orb
        ii=list_core(i)
        e_one_all+=2.D0*onetrf(ii,ii)
        do j=1,n_core_orb
         jj=list_core(j)
         e_two_all+=2.D0*bielec_PQxxtmp(ii,ii,j,j)-bielec_PQxxtmp(ii,jj,j,i)
        end do
        do t=1,n_act_orb
         tt=list_act(t)
         t3=t+n_core_orb
         do u=1,n_act_orb
          uu=list_act(u)
          u3=u+n_core_orb
          e_two_all+=D0tu(t,u)*(2.D0*bielec_PQxxtmp(tt,uu,i,i) &
              -bielec_PQxxtmp(tt,ii,i,u3))
         end do
        end do
       end do
       do t=1,n_act_orb
        tt=list_act(t)
        do u=1,n_act_orb
         uu=list_act(u)
         e_one_all+=D0tu(t,u)*onetrf(tt,uu)
         do v=1,n_act_orb
          v3=v+n_core_orb
          do x=1,n_act_orb
           x3=x+n_core_orb
           e_two_all  +=P0tuvx(t,u,v,x)*bielec_PQxxtmp(tt,uu,v3,x3)
          end do
         end do
        end do
       end do
       write(6,*) ' e_one_all = ',e_one_all
       write(6,*) ' e_two_all = ',e_two_all
       ecore    =nuclear_repulsion
       ecore_bis=nuclear_repulsion
       do i=1,n_core_orb
        ii=list_core(i)
        ecore    +=2.D0*onetrf(ii,ii)
        ecore_bis+=2.D0*onetrf(ii,ii)
        do j=1,n_core_orb
         jj=list_core(j)
         ecore    +=2.D0*bielec_PQxxtmp(ii,ii,j,j)-bielec_PQxxtmp(ii,jj,j,i)
         ecore_bis+=2.D0*bielec_PxxQtmp(ii,i,j,jj)-bielec_PxxQtmp(ii,j,j,ii)
        end do
       end do
       eone    =0.D0
       eone_bis=0.D0
       etwo    =0.D0
       etwo_bis=0.D0
       etwo_ter=0.D0
       do t=1,n_act_orb
        tt=list_act(t)
        t3=t+n_core_orb
        do u=1,n_act_orb
         uu=list_act(u)
         u3=u+n_core_orb
         eone    +=D0tu(t,u)*onetrf(tt,uu)
         eone_bis+=D0tu(t,u)*onetrf(tt,uu)
         do i=1,n_core_orb
          ii=list_core(i)
          eone    +=D0tu(t,u)*(2.D0*bielec_PQxxtmp(tt,uu,i,i) &
            -bielec_PQxxtmp(tt,ii,i,u3))
          eone_bis+=D0tu(t,u)*(2.D0*bielec_PxxQtmp(tt,u3,i,ii) &
            -bielec_PxxQtmp(tt,i,i,uu))
         end do
         do v=1,n_act_orb
          vv=list_act(v)
          v3=v+n_core_orb
          do x=1,n_act_orb
           xx=list_act(x)
           x3=x+n_core_orb
real*8 :: h1,h2,h3
           h1=bielec_PQxxtmp(tt,uu,v3,x3)
           h2=bielec_PxxQtmp(tt,u3,v3,xx)
           h3=bielecCItmp(t,u,v,xx)
           etwo    +=P0tuvx(t,u,v,x)*h1
           etwo_bis+=P0tuvx(t,u,v,x)*h2
           etwo_ter+=P0tuvx(t,u,v,x)*h3
           if ((abs(h1-h2).gt.1.D-14).or.(abs(h1-h3).gt.1.D-14)) then
                write(6,9901) t,u,v,x,h1,h2,h3
9901 format('aie: ',4I4,3E20.12)
           end if
          end do
         end do
        end do
       end do

       write(6,*) ' energy contributions '
       write(6,*) '     core energy       = ',ecore,' using PQxx integrals '
       write(6,*) '     core energy (bis) = ',ecore,' using PxxQ integrals '
       write(6,*) '     1el  energy       = ',eone ,' using PQxx integrals '
       write(6,*) '     1el  energy (bis) = ',eone ,' using PxxQ integrals '
       write(6,*) '     2el  energy       = ',etwo    ,' using PQxx integrals '
       write(6,*) '     2el  energy (bis) = ',etwo_bis,' using PxxQ integrals '
       write(6,*) '     2el  energy (ter) = ',etwo_ter,' using tuvP integrals '
       write(6,*) ' ----------------------------------------- '
       write(6,*) '     sum of all        = ',eone+etwo+ecore
       write(6,*)
       
     end subroutine trf_to_natorb

 BEGIN_PROVIDER [real*8, onetrf, (mo_num,mo_num)]
&BEGIN_PROVIDER [real*8, NatOrbsFCI, (ao_num,mo_num)]
END_PROVIDER
