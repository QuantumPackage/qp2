 BEGIN_PROVIDER [real*8, bielec_PQxx_no, (mo_num, mo_num,n_core_inact_orb+n_act_orb,n_core_inact_orb+n_act_orb)]
  BEGIN_DOC
  ! integral (pq|xx) in the basis of natural MOs
  ! indices are unshifted orbital numbers
  END_DOC
  implicit none
  integer                        :: i,j,k,l,t,u,p,q,pp
  real*8                         :: d(n_act_orb)

  bielec_PQxx_no(:,:,:,:) = bielec_PQxx(:,:,:,:)

  do j=1,mo_num
    do k=1,n_core_inact_orb+n_act_orb
      do l=1,n_core_inact_orb+n_act_orb
        do p=1,n_act_orb
          d(p)=0.D0
        end do
        do p=1,n_act_orb
          pp=n_act_orb-p+1
          do q=1,n_act_orb
            d(pp)+=bielec_PQxx_no(list_act(q),j,k,l)*natorbsCI(q,p)
          end do
        end do
        do p=1,n_act_orb
          bielec_PQxx_no(list_act(p),j,k,l)=d(p)
        end do
      end do
    end do
  end do
  ! 2nd quarter
  do j=1,mo_num
    do k=1,n_core_inact_orb+n_act_orb
      do l=1,n_core_inact_orb+n_act_orb
        do p=1,n_act_orb
          d(p)=0.D0
        end do
        do p=1,n_act_orb
          pp=n_act_orb-p+1
          do q=1,n_act_orb
            d(pp)+=bielec_PQxx_no(j,list_act(q),k,l)*natorbsCI(q,p)
          end do
        end do
        do p=1,n_act_orb
          bielec_PQxx_no(j,list_act(p),k,l)=d(p)
        end do
      end do
    end do
  end do
  ! 3rd quarter
  do j=1,mo_num
    do k=1,mo_num
      do l=1,n_core_inact_orb+n_act_orb
        do p=1,n_act_orb
          d(p)=0.D0
        end do
        do p=1,n_act_orb
          pp=n_act_orb-p+1
          do q=1,n_act_orb
            d(pp)+=bielec_PQxx_no(j,k,n_core_inact_orb+q,l)*natorbsCI(q,p)
          end do
        end do
        do p=1,n_act_orb
          bielec_PQxx_no(j,k,n_core_inact_orb+p,l)=d(p)
        end do
      end do
    end do
  end do
  ! 4th quarter
  do j=1,mo_num
    do k=1,mo_num
      do l=1,n_core_inact_orb+n_act_orb
        do p=1,n_act_orb
          d(p)=0.D0
        end do
        do p=1,n_act_orb
          pp=n_act_orb-p+1
          do q=1,n_act_orb
            d(pp)+=bielec_PQxx_no(j,k,l,n_core_inact_orb+q)*natorbsCI(q,p)
          end do
        end do
        do p=1,n_act_orb
          bielec_PQxx_no(j,k,l,n_core_inact_orb+p)=d(p)
        end do
      end do
    end do
  end do

END_PROVIDER



BEGIN_PROVIDER [real*8, bielec_PxxQ_no, (mo_num,n_core_inact_orb+n_act_orb,n_core_inact_orb+n_act_orb, mo_num)]
  BEGIN_DOC
  ! integral (px|xq) in the basis of natural MOs
  ! indices are unshifted orbital numbers
  END_DOC
  implicit none
  integer                        :: i,j,k,l,t,u,p,q,pp
  real*8                         :: d(n_act_orb)

  bielec_PxxQ_no(:,:,:,:) = bielec_PxxQ(:,:,:,:)

  do j=1,mo_num
    do k=1,n_core_inact_orb+n_act_orb
      do l=1,n_core_inact_orb+n_act_orb
        do p=1,n_act_orb
          d(p)=0.D0             
        end do
        do p=1,n_act_orb
          pp=n_act_orb-p+1
          do q=1,n_act_orb
            d(pp)+=bielec_PxxQ_no(list_act(q),k,l,j)*natorbsCI(q,p)
          end do
        end do
        do p=1,n_act_orb
          bielec_PxxQ_no(list_act(p),k,l,j)=d(p)
        end do
      end do
    end do
  end do
  ! 2nd quarter
  do j=1,mo_num
    do k=1,n_core_inact_orb+n_act_orb
      do l=1,n_core_inact_orb+n_act_orb
        do p=1,n_act_orb
          d(p)=0.D0
        end do
        do p=1,n_act_orb
          pp=n_act_orb-p+1
          do q=1,n_act_orb
            d(pp)+=bielec_PxxQ_no(j,k,l,list_act(q))*natorbsCI(q,p)
          end do
        end do
        do p=1,n_act_orb
          bielec_PxxQ_no(j,k,l,list_act(p))=d(p)
        end do
      end do
    end do
  end do
  ! 3rd quarter
  do j=1,mo_num
    do k=1,mo_num
      do l=1,n_core_inact_orb+n_act_orb
        do p=1,n_act_orb
          d(p)=0.D0
        end do
        do p=1,n_act_orb
          pp=n_act_orb-p+1
          do q=1,n_act_orb
            d(pp)+=bielec_PxxQ_no(j,n_core_inact_orb+q,l,k)*natorbsCI(q,p)
          end do
        end do
        do p=1,n_act_orb
          bielec_PxxQ_no(j,n_core_inact_orb+p,l,k)=d(p)
        end do
      end do
    end do
  end do
  ! 4th quarter
  do j=1,mo_num
    do k=1,mo_num
      do l=1,n_core_inact_orb+n_act_orb
        do p=1,n_act_orb
          d(p)=0.D0
        end do
        do p=1,n_act_orb
          pp=n_act_orb-p+1
          do q=1,n_act_orb
            d(pp)+=bielec_PxxQ_no(j,l,n_core_inact_orb+q,k)*natorbsCI(q,p)
          end do
        end do
        do p=1,n_act_orb
          bielec_PxxQ_no(j,l,n_core_inact_orb+p,k)=d(p)
        end do
      end do
    end do
  end do

END_PROVIDER


BEGIN_PROVIDER [real*8, bielecCI_no, (n_act_orb,n_act_orb,n_act_orb, mo_num)]
  BEGIN_DOC
  ! integrals (tu|vp) in the basis of natural MOs
  ! index p runs over the whole basis, t,u,v only over the active orbitals
  END_DOC
  implicit none
  integer                        :: i,j,k,l,t,u,p,q,pp
  real*8                         :: d(n_act_orb)
  
  bielecCI_no(:,:,:,:) = bielecCI(:,:,:,:)

  do j=1,n_act_orb
    do k=1,n_act_orb
      do l=1,mo_num
        do p=1,n_act_orb
          d(p)=0.D0
        end do
        do p=1,n_act_orb
          pp=n_act_orb-p+1
          do q=1,n_act_orb
            d(pp)+=bielecCI_no(q,j,k,l)*natorbsCI(q,p)
          end do
        end do
        do p=1,n_act_orb
          bielecCI_no(p,j,k,l)=d(p)
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
            d(pp)+=bielecCI_no(j,q,k,l)*natorbsCI(q,p)
          end do
        end do
        do p=1,n_act_orb
          bielecCI_no(j,p,k,l)=d(p)
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
            d(pp)+=bielecCI_no(j,k,q,l)*natorbsCI(q,p)
          end do
        end do
        do p=1,n_act_orb
          bielecCI_no(j,k,p,l)=d(p)
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
            d(pp)+=bielecCI_no(j,k,l,list_act(q))*natorbsCI(q,p)
          end do
        end do
        do p=1,n_act_orb
          bielecCI_no(j,k,l,list_act(p))=d(p)
        end do
      end do 
    end do  
  end do  

END_PROVIDER

