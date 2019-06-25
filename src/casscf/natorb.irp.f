 BEGIN_PROVIDER [real*8, occnum, (mo_num)]
  implicit none
  BEGIN_DOC
  ! MO occupation numbers
  END_DOC
  
  integer :: i
  occnum=0.D0
  do i=1,n_core_orb
    occnum(list_core(i))=2.D0
  end do
  
  do i=1,n_act_orb
    occnum(list_act(i))=occ_act(n_act_orb-i+1)
  end do

  write(6,*) ' occupation numbers '
  do i=1,mo_num
    write(6,*) i,occnum(i)
  end do

END_PROVIDER


 BEGIN_PROVIDER [ real*8, natorbsCI, (n_act_orb,n_act_orb) ]
&BEGIN_PROVIDER [ real*8, occ_act, (n_act_orb) ]
 implicit none
 BEGIN_DOC
 ! Natural orbitals of CI
 END_DOC
 integer                        :: i, j
 
 call lapack_diag(occ_act,natorbsCI,D0tu,n_act_orb,n_act_orb)
 
 write(6,*) ' found occupation numbers as '
 do i=1,n_act_orb
   write(6,*) i,occ_act(i)
 end do
 
 if (bavard) then
   !
   
   integer                        :: nmx
   real*8                         :: xmx
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
  
END_PROVIDER


BEGIN_PROVIDER [real*8, P0tuvx_no, (n_act_orb,n_act_orb,n_act_orb,n_act_orb)]
  implicit none
  BEGIN_DOC
  ! 4-index transformation of 2part matrices
  END_DOC
  integer                        :: i,j,k,l,p,q,pp
  real*8                         :: d(n_act_orb)

  ! index per index
  ! first quarter
  P0tuvx_no(:,:,:,:) = P0tuvx(:,:,:,:)

  do j=1,n_act_orb
    do k=1,n_act_orb
      do l=1,n_act_orb
        do p=1,n_act_orb
          d(p)=0.D0
        end do
        do p=1,n_act_orb
          pp=n_act_orb-p+1
          do q=1,n_act_orb
            d(pp)+=P0tuvx_no(q,j,k,l)*natorbsCI(q,p)
          end do
        end do
        do p=1,n_act_orb
          P0tuvx_no(p,j,k,l)=d(p)
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
            d(pp)+=P0tuvx_no(j,q,k,l)*natorbsCI(q,p)
          end do
        end do
        do p=1,n_act_orb
          P0tuvx_no(j,p,k,l)=d(p)
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
            d(pp)+=P0tuvx_no(j,k,q,l)*natorbsCI(q,p)
          end do
        end do
        do p=1,n_act_orb
          P0tuvx_no(j,k,p,l)=d(p)
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
            d(pp)+=P0tuvx_no(j,k,l,q)*natorbsCI(q,p)
          end do
        end do
        do p=1,n_act_orb
          P0tuvx_no(j,k,l,p)=d(p)
        end do
      end do
    end do
  end do
  write(6,*) ' transformed P0tuvx '

END_PROVIDER



BEGIN_PROVIDER [real*8, onetrf, (mo_num,mo_num)]
  implicit none
  BEGIN_DOC
  ! Transformed one-e integrals
  END_DOC
  integer :: i,j, p, pp, q
  real*8 :: d(n_act_orb)
  onetrf(:,:)=mo_one_e_integrals(:,:)

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
END_PROVIDER


BEGIN_PROVIDER [real*8, NatOrbsFCI, (ao_num,mo_num)]
  implicit none
  BEGIN_DOC
! FCI natural orbitals
  END_DOC
  integer :: i,j, p, pp, q
  real*8 :: d(n_act_orb)

  NatOrbsFCI(:,:)=mo_coef(:,:)
  
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
END_PROVIDER






subroutine trf_to_natorb()
  implicit none
  BEGIN_DOC
  ! save the diagonal somewhere, in inverse order
  ! 4-index-transform the 2-particle density matrix over active orbitals
  ! correct the bielectronic integrals
  ! correct the monoelectronic integrals
  ! put integrals on file, as well orbitals, and the density matrices
  !
  END_DOC
  integer                        :: i,j,k,l,t,u,p,q,pp
  real*8                         :: d(n_act_orb),d1(n_act_orb),d2(n_act_orb)
  
  ! we recalculate total energies
  write(6,*)
  write(6,*)  ' recalculating energies after the transformation '
  write(6,*)
  write(6,*)
  real*8                         :: e_one_all
  real*8                         :: e_two_all
  integer                        :: ii
  integer                        :: jj
  integer                        :: t3
  integer                        :: tt
  integer                        :: u3
  integer                        :: uu
  integer                        :: v
  integer                        :: v3
  integer                        :: vv
  integer                        :: x
  integer                        :: x3
  integer                        :: xx
  
  e_one_all=0.D0
  e_two_all=0.D0
  do i=1,n_core_orb
    ii=list_core(i)
    e_one_all+=2.D0*onetrf(ii,ii)
    do j=1,n_core_orb
      jj=list_core(j)
      e_two_all+=2.D0*bielec_PQxx_no(ii,ii,j,j)-bielec_PQxx_no(ii,jj,j,i)
    end do
    do t=1,n_act_orb
      tt=list_act(t)
      t3=t+n_core_orb
      e_two_all += occnum(list_act(t)) *                            &
          (2.d0*bielec_PQxx_no(tt,tt,i,i) - bielec_PQxx_no(tt,ii,i,t3))
    end do
  end do

    

  do t=1,n_act_orb
    tt=list_act(t)
    e_one_all += occnum(list_act(t))*onetrf(tt,tt)
    do u=1,n_act_orb
      uu=list_act(u)
      do v=1,n_act_orb
        v3=v+n_core_orb
        do x=1,n_act_orb
          x3=x+n_core_orb
          e_two_all  +=P0tuvx_no(t,u,v,x)*bielec_PQxx_no(tt,uu,v3,x3)
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
      ecore    +=2.D0*bielec_PQxx_no(ii,ii,j,j)-bielec_PQxx_no(ii,jj,j,i)
      ecore_bis+=2.D0*bielec_PxxQ_no(ii,i,j,jj)-bielec_PxxQ_no(ii,j,j,ii)
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
    eone     += occnum(list_act(t))*onetrf(tt,tt)
    eone_bis += occnum(list_act(t))*onetrf(tt,tt)
    do i=1,n_core_orb
      ii=list_core(i)
      eone     += occnum(list_act(t)) * &
        (2.D0*bielec_PQxx_no(tt,tt,i,i ) - bielec_PQxx_no(tt,ii,i,t3))
      eone_bis += occnum(list_act(t)) * &
        (2.D0*bielec_PxxQ_no(tt,t3,i,ii) - bielec_PxxQ_no(tt,i ,i,tt))
    end do
    do u=1,n_act_orb
      uu=list_act(u)
      u3=u+n_core_orb
      do v=1,n_act_orb
        vv=list_act(v)
        v3=v+n_core_orb
        do x=1,n_act_orb
          xx=list_act(x)
          x3=x+n_core_orb
          real*8                         :: h1,h2,h3
          h1=bielec_PQxx_no(tt,uu,v3,x3)
          h2=bielec_PxxQ_no(tt,u3,v3,xx)
          h3=bielecCI_no(t,u,v,xx)
          etwo    +=P0tuvx_no(t,u,v,x)*h1
          etwo_bis+=P0tuvx_no(t,u,v,x)*h2
          etwo_ter+=P0tuvx_no(t,u,v,x)*h3
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
  SOFT_TOUCH ecore ecore_bis eone eone_bis etwo etwo_bis etwo_ter 
  
end subroutine trf_to_natorb

