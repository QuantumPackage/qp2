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

  if (bavard) then
    write(6,*) ' occupation numbers '
    do i=1,mo_num
      write(6,*) i,occnum(i)
    end do
  endif

END_PROVIDER


 BEGIN_PROVIDER [ real*8, natorbsCI, (n_act_orb,n_act_orb) ]
&BEGIN_PROVIDER [ real*8, occ_act, (n_act_orb) ]
 implicit none
 BEGIN_DOC
 ! Natural orbitals of CI
 END_DOC
 integer                        :: i, j
 
 call lapack_diag(occ_act,natorbsCI,D0tu,n_act_orb,n_act_orb)
 
 if (bavard) then
   write(6,*) ' found occupation numbers as '
   do i=1,n_act_orb
     write(6,*) i,occ_act(i)
   end do
 
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

END_PROVIDER



BEGIN_PROVIDER [real*8, one_ints_no, (mo_num,mo_num)]
  implicit none
  BEGIN_DOC
  ! Transformed one-e integrals
  END_DOC
  integer :: i,j, p, pp, q
  real*8 :: d(n_act_orb)
  one_ints_no(:,:)=mo_one_e_integrals(:,:)

  ! 1st half-trf
  do j=1,mo_num
    do p=1,n_act_orb
      d(p)=0.D0
    end do
    do p=1,n_act_orb
      pp=n_act_orb-p+1
      do q=1,n_act_orb
        d(pp)+=one_ints_no(list_act(q),j)*natorbsCI(q,p)
      end do
    end do
    do p=1,n_act_orb
      one_ints_no(list_act(p),j)=d(p)
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
        d(pp)+=one_ints_no(j,list_act(q))*natorbsCI(q,p)
      end do
    end do
    do p=1,n_act_orb
      one_ints_no(j,list_act(p))=d(p)
    end do
  end do
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
END_PROVIDER

