BEGIN_PROVIDER [real*8, bielec_PQxx, (mo_num, mo_num,n_core_inact_act_orb,n_core_inact_act_orb)]
  BEGIN_DOC
  ! bielec_PQxx : integral (pq|xx) with p,q arbitrary, x core or active
  ! indices are unshifted orbital numbers
  END_DOC
  implicit none
  integer                        :: i,j,ii,jj,p,q,i3,j3,t3,v3
  real*8                         :: mo_two_e_integral
  
  bielec_PQxx(:,:,:,:) = 0.d0
  PROVIDE mo_two_e_integrals_in_map
  
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(i,ii,j,jj,i3,j3) &
  !$OMP SHARED(n_core_inact_orb,list_core_inact,mo_num,bielec_PQxx, &
  !$OMP  n_act_orb,mo_integrals_map,list_act)

  !$OMP DO
  do i=1,n_core_inact_orb
    ii=list_core_inact(i)
    do j=i,n_core_inact_orb
      jj=list_core_inact(j)
      call get_mo_two_e_integrals_i1j1(ii,jj,mo_num,bielec_PQxx(1,1,i,j),mo_integrals_map)
      bielec_PQxx(:,:,j,i)=bielec_PQxx(:,:,i,j)
    end do
    do j=1,n_act_orb
      jj=list_act(j)
      j3=j+n_core_inact_orb
      call get_mo_two_e_integrals_i1j1(ii,jj,mo_num,bielec_PQxx(1,1,i,j3),mo_integrals_map)
      bielec_PQxx(:,:,j3,i)=bielec_PQxx(:,:,i,j3)
    end do
  end do
  !$OMP END DO


  !$OMP DO
  do i=1,n_act_orb
    ii=list_act(i)
    i3=i+n_core_inact_orb
    do j=i,n_act_orb
      jj=list_act(j)
      j3=j+n_core_inact_orb
      call get_mo_two_e_integrals_i1j1(ii,jj,mo_num,bielec_PQxx(1,1,i3,j3),mo_integrals_map)
      bielec_PQxx(:,:,j3,i3)=bielec_PQxx(:,:,i3,j3)
    end do
  end do
  !$OMP END DO

  !$OMP END PARALLEL

END_PROVIDER



BEGIN_PROVIDER [real*8, bielec_PxxQ, (mo_num,n_core_inact_act_orb,n_core_inact_act_orb, mo_num)]
  BEGIN_DOC
  ! bielec_PxxQ : integral (px|xq) with p,q arbitrary, x core or active
  ! indices are unshifted orbital numbers
  END_DOC
  implicit none
  integer                        :: i,j,ii,jj,p,q,i3,j3,t3,v3
  double precision, allocatable  :: integrals_array(:,:)
  real*8                         :: mo_two_e_integral
  
  PROVIDE mo_two_e_integrals_in_map
  bielec_PxxQ = 0.d0

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(i,ii,j,jj,i3,j3,integrals_array) &
  !$OMP SHARED(n_core_inact_orb,list_core_inact,mo_num,bielec_PxxQ, &
  !$OMP  n_act_orb,mo_integrals_map,list_act)

  allocate(integrals_array(mo_num,mo_num))
  
  !$OMP DO
  do i=1,n_core_inact_orb
    ii=list_core_inact(i)
    do j=i,n_core_inact_orb
      jj=list_core_inact(j)
      call get_mo_two_e_integrals_ij(ii,jj,mo_num,integrals_array,mo_integrals_map)
      do q=1,mo_num
        do p=1,mo_num
          bielec_PxxQ(p,i,j,q)=integrals_array(p,q)
          bielec_PxxQ(p,j,i,q)=integrals_array(q,p)
        end do
      end do
    end do
    do j=1,n_act_orb
      jj=list_act(j)
      j3=j+n_core_inact_orb
      call get_mo_two_e_integrals_ij(ii,jj,mo_num,integrals_array,mo_integrals_map)
      do q=1,mo_num
        do p=1,mo_num
          bielec_PxxQ(p,i,j3,q)=integrals_array(p,q)
          bielec_PxxQ(p,j3,i,q)=integrals_array(q,p)
        end do
      end do
    end do
  end do
  !$OMP END DO


  ! (ip|qj)
  !$OMP DO
  do i=1,n_act_orb
    ii=list_act(i)
    i3=i+n_core_inact_orb
    do j=i,n_act_orb
      jj=list_act(j)
      j3=j+n_core_inact_orb
      call get_mo_two_e_integrals_ij(ii,jj,mo_num,integrals_array,mo_integrals_map)
      do q=1,mo_num
        do p=1,mo_num
          bielec_PxxQ(p,i3,j3,q)=integrals_array(p,q)
          bielec_PxxQ(p,j3,i3,q)=integrals_array(q,p)
        end do
      end do
    end do
  end do
  !$OMP END DO

  deallocate(integrals_array)
  !$OMP END PARALLEL

END_PROVIDER


BEGIN_PROVIDER [real*8, bielecCI, (n_act_orb,n_act_orb,n_act_orb, mo_num)]
  BEGIN_DOC
  ! bielecCI : integrals (tu|vp) with p arbitrary, tuv active
  ! index p runs over the whole basis, t,u,v only over the active orbitals
  END_DOC
  implicit none
  integer                        :: i,j,k,p,t,u,v
  double precision, external     :: mo_two_e_integral
  PROVIDE mo_two_e_integrals_in_map
  
  !$OMP PARALLEL DO DEFAULT(NONE) &
  !$OMP PRIVATE(i,j,k,p,t,u,v) &
  !$OMP SHARED(mo_num,n_act_orb,list_act,bielecCI)
  do p=1,mo_num
    do j=1,n_act_orb
      u=list_act(j)
      do k=1,n_act_orb
        v=list_act(k)
        do i=1,n_act_orb
          t=list_act(i)
          bielecCI(i,k,j,p) = mo_two_e_integral(t,u,v,p)
        end do
      end do
    end do
  end do
  !$OMP END PARALLEL DO

END_PROVIDER
