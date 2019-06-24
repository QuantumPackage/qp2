! -*- F90 -*- 
 BEGIN_PROVIDER[real*8, bielec_PQxxtmp, (mo_num, mo_num,n_core_orb+n_act_orb,n_core_orb+n_act_orb)]
&BEGIN_PROVIDER[real*8, bielec_PxxQtmp, (mo_num,n_core_orb+n_act_orb,n_core_orb+n_act_orb, mo_num)]
&BEGIN_PROVIDER[integer, num_PQxx_written]
&BEGIN_PROVIDER[integer, num_PxxQ_written]
BEGIN_DOC
! bielec_PQxx : integral (pq|xx) with p,q arbitrary, x core or active
! bielec_PxxQ : integral (px|xq) with p,q arbitrary, x core or active
! indices are unshifted orbital numbers
END_DOC
        implicit none
          integer :: i,j,ii,jj,p,q,i3,j3,t3,v3
          double precision, allocatable :: integrals_array1(:,:)
          double precision, allocatable :: integrals_array2(:,:)
          real*8 :: mo_two_e_integral

          allocate(integrals_array1(mo_num,mo_num))
          allocate(integrals_array2(mo_num,mo_num))

          do i=1,n_core_orb+n_act_orb
           do j=1,n_core_orb+n_act_orb
            do p=1,mo_num
             do q=1,mo_num
              bielec_PQxxtmp(p,q,i,j)=0.D0
              bielec_PxxQtmp(p,i,j,q)=0.D0
             end do
            end do
           end do
          end do

          do i=1,n_core_orb
           ii=list_core(i)
           do j=i,n_core_orb
            jj=list_core(j)
! (ij|pq)
            call get_mo_two_e_integrals_i1j1(ii,jj,mo_num,integrals_array1,mo_integrals_map)
! (ip|qj)
            call get_mo_two_e_integrals_ij  (ii,jj,mo_num,integrals_array2,mo_integrals_map)
            do p=1,mo_num
             do q=1,mo_num
              bielec_PQxxtmp(p,q,i,j)=integrals_array1(p,q)
              bielec_PQxxtmp(p,q,j,i)=integrals_array1(p,q)
              bielec_PxxQtmp(p,i,j,q)=integrals_array2(p,q)
              bielec_PxxQtmp(p,j,i,q)=integrals_array2(q,p)
             end do
            end do
           end do
           do j=1,n_act_orb
            jj=list_act(j)
            j3=j+n_core_orb
! (ij|pq)
            call get_mo_two_e_integrals_i1j1(ii,jj,mo_num,integrals_array1,mo_integrals_map)
! (ip|qj)
            call get_mo_two_e_integrals_ij  (ii,jj,mo_num,integrals_array2,mo_integrals_map)
            do p=1,mo_num
             do q=1,mo_num
              bielec_PQxxtmp(p,q,i,j3)=integrals_array1(p,q)
              bielec_PQxxtmp(p,q,j3,i)=integrals_array1(p,q)
              bielec_PxxQtmp(p,i,j3,q)=integrals_array2(p,q)
              bielec_PxxQtmp(p,j3,i,q)=integrals_array2(q,p)
             end do
            end do
           end do
          end do
          do i=1,n_act_orb
           ii=list_act(i)
           i3=i+n_core_orb
           do j=i,n_act_orb
            jj=list_act(j)
            j3=j+n_core_orb
! (ij|pq)
            call get_mo_two_e_integrals_i1j1(ii,jj,mo_num,integrals_array1,mo_integrals_map)
! (ip|qj)
            call get_mo_two_e_integrals_ij  (ii,jj,mo_num,integrals_array2,mo_integrals_map)
            do p=1,mo_num
             do q=1,mo_num
              bielec_PQxxtmp(p,q,i3,j3)=integrals_array1(p,q)
              bielec_PQxxtmp(p,q,j3,i3)=integrals_array1(p,q)
              bielec_PxxQtmp(p,i3,j3,q)=integrals_array2(p,q)
              bielec_PxxQtmp(p,j3,i3,q)=integrals_array2(q,p)
             end do
            end do
           end do
          end do
          write(6,*) ' provided integrals (PQ|xx) '
          write(6,*) ' provided integrals (Px|xQ) '
!!$          write(6,*) ' 1 1 1 2 = ',bielec_PQxxtmp(2,2,2,3),bielec_PxxQtmp(2,2,2,3)
END_PROVIDER

BEGIN_PROVIDER[real*8, bielecCItmp, (n_act_orb,n_act_orb,n_act_orb, mo_num)]
BEGIN_DOC
! bielecCI : integrals (tu|vp) with p arbitrary, tuv active
! index p runs over the whole basis, t,u,v only over the active orbitals
END_DOC
        implicit none
          integer :: i,j,k,p,t,u,v
          double precision, allocatable :: integrals_array1(:)
          real*8 :: mo_two_e_integral

          allocate(integrals_array1(mo_num))
 
          do i=1,n_act_orb
           t=list_act(i)
           do j=1,n_act_orb
            u=list_act(j)
            do k=1,n_act_orb
             v=list_act(k)
! (tu|vp)
             call get_mo_two_e_integrals(t,u,v,mo_num,integrals_array1,mo_integrals_map)
             do p=1,mo_num
              bielecCItmp(i,k,j,p)=integrals_array1(p)
             end do
            end do
           end do
          end do
          write(6,*) ' provided integrals (tu|xP) '
END_PROVIDER

