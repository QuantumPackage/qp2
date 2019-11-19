BEGIN_PROVIDER [ logical, no_vvvv_integrals  ]
  implicit none
  BEGIN_DOC
! If `True`, computes all integrals except for the integrals having 3 or 4 virtual indices
  END_DOC

  no_vvvv_integrals = .False.
END_PROVIDER

BEGIN_PROVIDER [ double precision, mo_coef_novirt, (ao_num,n_core_inact_act_orb) ]
 implicit none
 BEGIN_DOC
 ! MO coefficients without virtual MOs
 END_DOC
 integer :: j,jj

 do j=1,n_core_inact_act_orb
   jj = list_core_inact_act(j)
   mo_coef_novirt(:,j) = mo_coef(:,jj)
 enddo

END_PROVIDER

subroutine ao_to_mo_novirt(A_ao,LDA_ao,A_mo,LDA_mo)
  implicit none
  BEGIN_DOC
  ! Transform A from the |AO| basis to the |MO| basis excluding virtuals
  !
  ! $C^\dagger.A_{ao}.C$
  END_DOC
  integer, intent(in)            :: LDA_ao,LDA_mo
  double precision, intent(in)   :: A_ao(LDA_ao,ao_num)
  double precision, intent(out)  :: A_mo(LDA_mo,n_core_inact_act_orb)
  double precision, allocatable  :: T(:,:)

  allocate ( T(ao_num,n_core_inact_act_orb) )
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: T

  call dgemm('N','N', ao_num, n_core_inact_act_orb, ao_num,          &
      1.d0, A_ao,LDA_ao,                                             &
      mo_coef_novirt, size(mo_coef_novirt,1),                        &
      0.d0, T, size(T,1))
  
  call dgemm('T','N', n_core_inact_act_orb, n_core_inact_act_orb, ao_num,&
      1.d0, mo_coef_novirt,size(mo_coef_novirt,1),                   &
      T, ao_num,                                                     &
      0.d0, A_mo, size(A_mo,1))
  
  deallocate(T)
end


subroutine four_idx_novvvv
  use map_module
  implicit none
  BEGIN_DOC
  ! Retransform MO integrals for next CAS-SCF step
  END_DOC
  integer                        :: i,j,k,l,n_integrals
  double precision, allocatable  :: f(:,:,:), f2(:,:,:), d(:,:), T(:,:,:,:), T2(:,:,:,:)
  double precision, external     :: get_ao_two_e_integral
  integer(key_kind), allocatable :: idx(:)
  real(integral_kind), allocatable :: values(:)

  integer                        :: p,q,r,s
  double precision               :: c
  allocate( T(n_core_inact_act_orb,n_core_inact_act_orb,ao_num,ao_num) , &
            T2(n_core_inact_act_orb,n_core_inact_act_orb,ao_num,ao_num) )
  
  !$OMP PARALLEL DEFAULT(NONE)                                       &
  !$OMP SHARED(mo_num,ao_num,T,n_core_inact_act_orb, mo_coef_transp, &
  !$OMP   mo_integrals_threshold,mo_coef,mo_integrals_map,           &
  !$OMP   list_core_inact_act,T2,ao_integrals_map)            &
  !$OMP PRIVATE(i,j,k,l,p,q,r,s,idx,values,n_integrals,              &
  !$OMP   f,f2,d,c)
  allocate(f(ao_num,ao_num,ao_num), f2(ao_num,ao_num,ao_num), d(mo_num,mo_num), &
           idx(mo_num*mo_num), values(mo_num*mo_num) )
  
  ! <aa|vv>
  !$OMP DO
  do s=1,ao_num
    do r=1,ao_num
      do q=1,ao_num
        do p=1,r
          f (p,q,r) = get_ao_two_e_integral(p,q,r,s,ao_integrals_map)
          f (r,q,p) = f(p,q,r)
        enddo
      enddo
    enddo
    do r=1,ao_num
      do q=1,ao_num
        do p=1,ao_num
          f2(p,q,r) = f(p,r,q)
        enddo
      enddo
    enddo
    ! f (p,q,r) = <pq|rs>
    ! f2(p,q,r) = <pr|qs>

    do r=1,ao_num
      call ao_to_mo_novirt(f (1,1,r),size(f ,1),T (1,1,r,s),size(T,1))
      call ao_to_mo_novirt(f2(1,1,r),size(f2,1),T2(1,1,r,s),size(T,1))
    enddo
    ! T (i,j,p,q) = <ij|rs>
    ! T2(i,j,p,q) = <ir|js>

  enddo
  !$OMP END DO 

  !$OMP DO
  do j=1,n_core_inact_act_orb
    do i=1,n_core_inact_act_orb
      do s=1,ao_num
        do r=1,ao_num
          f (r,s,1) = T (i,j,r,s)
          f2(r,s,1) = T2(i,j,r,s)
        enddo
      enddo
      call ao_to_mo(f ,size(f ,1),d,size(d,1))
      n_integrals = 0
      do l=1,mo_num
        do k=1,mo_num
          n_integrals+=1
          call two_e_integrals_index(list_core_inact_act(i),list_core_inact_act(j),k,l,idx(n_integrals))
          values(n_integrals) = d(k,l)
        enddo
      enddo
      call map_append(mo_integrals_map, idx, values, n_integrals)

      call ao_to_mo(f2,size(f2,1),d,size(d,1))
      n_integrals = 0
      do l=1,mo_num
        do k=1,mo_num
          n_integrals+=1
          call two_e_integrals_index(list_core_inact_act(i),k,list_core_inact_act(j),l,idx(n_integrals))
          values(n_integrals) = d(k,l)
        enddo
      enddo
      call map_append(mo_integrals_map, idx, values, n_integrals)
    enddo
  enddo
  !$OMP END DO
  deallocate(f,f2,d,idx,values)
  
  !$OMP END PARALLEL

  deallocate(T,T2)
  
  
  call map_sort(mo_integrals_map)
  call map_unique(mo_integrals_map)
  call map_shrink(mo_integrals_map,real(mo_integrals_threshold,integral_kind))
  
end

subroutine four_idx_novvvv2
  use bitmasks
  implicit none
  integer                        :: i
  integer(bit_kind)              :: mask_ijkl(N_int,4)

    print*, '<ix|ix>'
    do i = 1,N_int
      mask_ijkl(i,1) =  core_inact_act_bitmask_4(i,1)
      mask_ijkl(i,2) =  full_ijkl_bitmask_4(i,1)
      mask_ijkl(i,3) =  core_inact_act_bitmask_4(i,1)
      mask_ijkl(i,4) =  full_ijkl_bitmask_4(i,1)
    enddo
    call add_integrals_to_map(mask_ijkl)

    print*, '<ii|vv>'
    do i = 1,N_int
      mask_ijkl(i,1) = core_inact_act_bitmask_4(i,1)
      mask_ijkl(i,2) = core_inact_act_bitmask_4(i,1)
      mask_ijkl(i,3) = virt_bitmask(i,1)
      mask_ijkl(i,4) = virt_bitmask(i,1)
    enddo
    call add_integrals_to_map(mask_ijkl)

end
