 BEGIN_PROVIDER [real*8, occnum, (mo_num)]
  implicit none
  BEGIN_DOC
  ! MO occupation numbers
  END_DOC
  
  integer :: i
  occnum=0.D0
  do i=1,n_core_inact_orb
    occnum(list_core_inact(i))=2.D0
  end do
  
  do i=1,n_act_orb
    occnum(list_act(i))=occ_act(i)
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
 double precision               :: Vt(n_act_orb,n_act_orb)
 
! call lapack_diag(occ_act,natorbsCI,D0tu,n_act_orb,n_act_orb)
 call svd(D0tu, size(D0tu,1), natorbsCI,size(natorbsCI,1), occ_act, Vt, size(Vt,1),n_act_orb,n_act_orb) 
 
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
  integer                        :: i,j,k,l,p,q
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
          do q=1,n_act_orb
            d(p)+=P0tuvx_no(q,j,k,l)*natorbsCI(q,p)
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
          do q=1,n_act_orb
            d(p)+=P0tuvx_no(j,q,k,l)*natorbsCI(q,p)
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
          do q=1,n_act_orb
            d(p)+=P0tuvx_no(j,k,q,l)*natorbsCI(q,p)
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
          do q=1,n_act_orb
            d(p)+=P0tuvx_no(j,k,l,q)*natorbsCI(q,p)
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
  integer :: i,j, p, q
  real*8 :: d(n_act_orb)
  one_ints_no(:,:)=mo_one_e_integrals(:,:)

  ! 1st half-trf
  do j=1,mo_num
    do p=1,n_act_orb
      d(p)=0.D0
    end do
    do p=1,n_act_orb
      do q=1,n_act_orb
        d(p)+=one_ints_no(list_act(q),j)*natorbsCI(q,p)
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
      do q=1,n_act_orb
        d(p)+=one_ints_no(j,list_act(q))*natorbsCI(q,p)
      end do
    end do
    do p=1,n_act_orb
      one_ints_no(j,list_act(p))=d(p)
    end do
  end do
END_PROVIDER


BEGIN_PROVIDER [ double precision, NatOrbsCI_mos, (mo_num, mo_num) ]
  implicit none
  BEGIN_DOC
  ! Rotation matrix from current MOs to the CI natural MOs
  END_DOC
  integer :: p,q

  NatOrbsCI_mos(:,:) = 0.d0

  do q = 1,mo_num
   NatOrbsCI_mos(q,q) = 1.d0
  enddo

  do q = 1,n_act_orb
    do p = 1,n_act_orb
      NatOrbsCI_mos(list_act(p),list_act(q)) = natorbsCI(p,q)
    enddo
  enddo
END_PROVIDER


BEGIN_PROVIDER [real*8, NatOrbsFCI, (ao_num,mo_num)]
  implicit none
  BEGIN_DOC
! FCI natural orbitals
  END_DOC

  call dgemm('N','N', ao_num,mo_num,mo_num,1.d0,                     &
      mo_coef, size(mo_coef,1),                                      &
      NatOrbsCI_mos, size(NatOrbsCI_mos,1), 0.d0,                    &
      NatOrbsFCI, size(NatOrbsFCI,1))
END_PROVIDER

