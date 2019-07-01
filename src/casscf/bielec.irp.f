 BEGIN_PROVIDER [real*8, bielec_PQxx, (mo_num, mo_num,n_core_orb+n_act_orb,n_core_orb+n_act_orb)]
  BEGIN_DOC
  ! bielec_PQxx : integral (pq|xx) with p,q arbitrary, x core or active
  ! indices are unshifted orbital numbers
  END_DOC
  implicit none
  integer                        :: i,j,ii,jj,p,q,i3,j3,t3,v3
  double precision, allocatable  :: integrals_array(:,:)
  real*8                         :: mo_two_e_integral
  
  allocate(integrals_array(mo_num,mo_num))
  
  bielec_PQxx = 0.d0
  
  do i=1,n_core_orb
    ii=list_core(i)
    do j=i,n_core_orb
      jj=list_core(j)
      call get_mo_two_e_integrals_i1j1(ii,jj,mo_num,integrals_array,mo_integrals_map)
      do p=1,mo_num
        do q=1,mo_num
          bielec_PQxx(p,q,i,j)=integrals_array(p,q)
          bielec_PQxx(p,q,j,i)=integrals_array(p,q)
        end do
      end do
    end do
    do j=1,n_act_orb
      jj=list_act(j)
      j3=j+n_core_orb
      call get_mo_two_e_integrals_i1j1(ii,jj,mo_num,integrals_array,mo_integrals_map)
      do p=1,mo_num
        do q=1,mo_num
          bielec_PQxx(p,q,i,j3)=integrals_array(p,q)
          bielec_PQxx(p,q,j3,i)=integrals_array(p,q)
        end do
      end do
    end do
  end do


  ! (ij|pq)
  do i=1,n_act_orb
    ii=list_act(i)
    i3=i+n_core_orb
    do j=i,n_act_orb
      jj=list_act(j)
      j3=j+n_core_orb
      call get_mo_two_e_integrals_i1j1(ii,jj,mo_num,integrals_array,mo_integrals_map)
      do p=1,mo_num
        do q=1,mo_num
          bielec_PQxx(p,q,i3,j3)=integrals_array(p,q)
          bielec_PQxx(p,q,j3,i3)=integrals_array(p,q)
        end do
      end do
    end do
  end do

END_PROVIDER



BEGIN_PROVIDER [real*8, bielec_PxxQ, (mo_num,n_core_orb+n_act_orb,n_core_orb+n_act_orb, mo_num)]
  BEGIN_DOC
  ! bielec_PxxQ : integral (px|xq) with p,q arbitrary, x core or active
  ! indices are unshifted orbital numbers
  END_DOC
  implicit none
  integer                        :: i,j,ii,jj,p,q,i3,j3,t3,v3
  double precision, allocatable  :: integrals_array(:,:)
  real*8                         :: mo_two_e_integral
  
  allocate(integrals_array(mo_num,mo_num))
  
  bielec_PxxQ = 0.d0
  
  do i=1,n_core_orb
    ii=list_core(i)
    do j=i,n_core_orb
      jj=list_core(j)
      call get_mo_two_e_integrals_ij  (ii,jj,mo_num,integrals_array,mo_integrals_map)
      do p=1,mo_num
        do q=1,mo_num
          bielec_PxxQ(p,i,j,q)=integrals_array(p,q)
          bielec_PxxQ(p,j,i,q)=integrals_array(q,p)
        end do
      end do
    end do
    do j=1,n_act_orb
      jj=list_act(j)
      j3=j+n_core_orb
      call get_mo_two_e_integrals_ij  (ii,jj,mo_num,integrals_array,mo_integrals_map)
      do p=1,mo_num
        do q=1,mo_num
          bielec_PxxQ(p,i,j3,q)=integrals_array(p,q)
          bielec_PxxQ(p,j3,i,q)=integrals_array(q,p)
        end do
      end do
    end do
  end do


  ! (ip|qj)
  do i=1,n_act_orb
    ii=list_act(i)
    i3=i+n_core_orb
    do j=i,n_act_orb
      jj=list_act(j)
      j3=j+n_core_orb
      call get_mo_two_e_integrals_ij  (ii,jj,mo_num,integrals_array,mo_integrals_map)
      do p=1,mo_num
        do q=1,mo_num
          bielec_PxxQ(p,i3,j3,q)=integrals_array(p,q)
          bielec_PxxQ(p,j3,i3,q)=integrals_array(q,p)
        end do
      end do
    end do
  end do
END_PROVIDER


BEGIN_PROVIDER [real*8, bielecCI, (n_act_orb,n_act_orb,n_act_orb, mo_num)]
  BEGIN_DOC
  ! bielecCI : integrals (tu|vp) with p arbitrary, tuv active
  ! index p runs over the whole basis, t,u,v only over the active orbitals
  END_DOC
  implicit none
  integer                        :: i,j,k,p,t,u,v
  double precision, allocatable  :: integrals_array(:)
  real*8                         :: mo_two_e_integral
  
  allocate(integrals_array(mo_num))
  
  do i=1,n_act_orb
    t=list_act(i)
    do j=1,n_act_orb
      u=list_act(j)
      do k=1,n_act_orb
        v=list_act(k)
        ! (tu|vp)
        call get_mo_two_e_integrals(t,u,v,mo_num,integrals_array,mo_integrals_map)
        do p=1,mo_num
          bielecCI(i,k,j,p)=integrals_array(p)
        end do
      end do
    end do
  end do
END_PROVIDER

