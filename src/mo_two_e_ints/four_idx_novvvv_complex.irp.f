BEGIN_PROVIDER [ complex*16, mo_coef_novirt_complex, (ao_num,n_core_inact_act_orb) ]
 implicit none
 BEGIN_DOC
 ! MO coefficients without virtual MOs
 END_DOC
 integer :: j,jj

 do j=1,n_core_inact_act_orb
   jj = list_core_inact_act(j)
   mo_coef_novirt_complex(:,j) = mo_coef_complex(:,jj)
 enddo

END_PROVIDER

subroutine ao_to_mo_novirt_complex(A_ao,LDA_ao,A_mo,LDA_mo)
  implicit none
  BEGIN_DOC
  ! Transform A from the |AO| basis to the |MO| basis excluding virtuals
  !
  ! $C^\dagger.A_{ao}.C$
  END_DOC
  integer, intent(in)            :: LDA_ao,LDA_mo
  complex*16, intent(in)   :: A_ao(LDA_ao,ao_num)
  complex*16, intent(out)  :: A_mo(LDA_mo,n_core_inact_act_orb)
  complex*16, allocatable  :: T(:,:)

  allocate ( T(ao_num,n_core_inact_act_orb) )
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: T

  call zgemm('N','N', ao_num, n_core_inact_act_orb, ao_num,          &
      (1.d0,0.d0), A_ao,LDA_ao,                                             &
      mo_coef_novirt_complex, size(mo_coef_novirt_complex,1),                        &
      (0.d0,0.d0), T, size(T,1))
  
  call zgemm('C','N', n_core_inact_act_orb, n_core_inact_act_orb, ao_num,&
      (1.d0,0.d0), mo_coef_novirt_complex,size(mo_coef_novirt_complex,1),  &
      T, ao_num,                                                     &
      (0.d0,0.d0), A_mo, size(A_mo,1))
  
  deallocate(T)
end

subroutine ao_to_mo_novirt_conjg_complex(A_ao,LDA_ao,A_mo,LDA_mo)
  implicit none
  BEGIN_DOC
  ! Transform A from the |AO| basis to the |MO| basis excluding virtuals
  !
  ! $C^\dagger.A_{ao}.C^*$
  ! half-transformed ints as handled by four_idx_novvvv need to use this
  END_DOC
  integer, intent(in)            :: LDA_ao,LDA_mo
  complex*16, intent(in)   :: A_ao(LDA_ao,ao_num)
  complex*16, intent(out)  :: A_mo(LDA_mo,n_core_inact_act_orb)
  complex*16, allocatable  :: T(:,:)

  allocate ( T(ao_num,n_core_inact_act_orb) )
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: T

  call zgemm('N','N', ao_num, n_core_inact_act_orb, ao_num,          &
      (1.d0,0.d0), A_ao,LDA_ao,                                             &
      dconjg(mo_coef_novirt_complex), size(mo_coef_novirt_complex,1),                        &
      (0.d0,0.d0), T, size(T,1))
  
  call zgemm('C','N', n_core_inact_act_orb, n_core_inact_act_orb, ao_num,&
      (1.d0,0.d0), mo_coef_novirt_complex,size(mo_coef_novirt_complex,1),  &
      T, ao_num,                                                     &
      (0.d0,0.d0), A_mo, size(A_mo,1))
  
  deallocate(T)
end


subroutine four_idx_novvvv_complex
  use map_module
  implicit none
  BEGIN_DOC
  ! Retransform MO integrals for next CAS-SCF step
  END_DOC
  integer                        :: i,j,k,l,n_integrals1,n_integrals2
  logical                        :: use_map1
  complex*16, allocatable        :: f(:,:,:), f2(:,:,:), d(:,:), T(:,:,:,:), T2(:,:,:,:)
  complex*16, external           :: get_ao_two_e_integral_complex
  integer(key_kind), allocatable :: idx1(:),idx2(:)
  complex(integral_kind), allocatable :: values1(:),values2(:)
  double precision               :: sign_tmp
  integer(key_kind)              :: idx_tmp

  integer                        :: p,q,r,s
  allocate( T(n_core_inact_act_orb,n_core_inact_act_orb,ao_num,ao_num) , &
            T2(n_core_inact_act_orb,n_core_inact_act_orb,ao_num,ao_num) )
  
  !$OMP PARALLEL DEFAULT(NONE)                                       &
  !$OMP SHARED(mo_num,ao_num,T,n_core_inact_act_orb, &
  !$OMP   mo_integrals_threshold,mo_integrals_map,           &
  !$OMP   mo_integrals_map_2,ao_integrals_map_2,           &
  !$OMP   list_core_inact_act,T2,ao_integrals_map)            &
  !$OMP PRIVATE(i,j,k,l,p,q,r,s,idx1,idx2,values1,values2,n_integrals1, &
  !$OMP   n_integrals2,use_map1,idx_tmp,sign_tmp,              &
  !$OMP   f,f2,d)
  allocate(f(ao_num,ao_num,ao_num), f2(ao_num,ao_num,ao_num), d(mo_num,mo_num), &
           idx1(2*mo_num*mo_num), values1(2*mo_num*mo_num), &
           idx2(2*mo_num*mo_num), values2(2*mo_num*mo_num) )
  
  ! <aa|vv>
  !$OMP DO
  do s=1,ao_num
    do r=1,ao_num
      do q=1,ao_num
        do p=1,r
          f (p,q,r) = get_ao_two_e_integral_complex(p,q,r,s,ao_integrals_map,ao_integrals_map_2)
          f (r,q,p) = get_ao_two_e_integral_complex(r,q,p,s,ao_integrals_map,ao_integrals_map_2)
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
      call ao_to_mo_novirt_conjg_complex(f (1,1,r),size(f ,1),T (1,1,r,s),size(T,1))
      call ao_to_mo_novirt_complex(f2(1,1,r),size(f2,1),T2(1,1,r,s),size(T,1))
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
      call ao_to_mo_noconjg_complex(f ,size(f ,1),d,size(d,1))
      n_integrals1 = 0
      n_integrals2 = 0
      do l=1,mo_num
        do k=1,mo_num
          call ao_two_e_integral_complex_map_idx_sign(list_core_inact_act(i),list_core_inact_act(j),k,l,use_map1,idx_tmp,sign_tmp)
          if (use_map1) then
            n_integrals1+=1
            values1(n_integrals1) = dble(d(k,l))
            idx1(n_integrals1) = idx_tmp
            if (sign_tmp /= 0.d0) then ! should always be true, but might change in the future
              n_integrals1+=1
              values1(n_integrals1) = sign_tmp*dimag(d(k,l))
              idx1(n_integrals1) = idx_tmp+1
            endif
          else
            n_integrals2+=1
            values2(n_integrals2) = dble(d(k,l))
            idx2(n_integrals2) = idx_tmp
            if (sign_tmp /= 0.d0) then
              n_integrals2+=1
              values2(n_integrals2) = sign_tmp*dimag(d(k,l))
              idx2(n_integrals2) = idx_tmp+1
            endif
          endif
        enddo
      enddo
      call map_append(mo_integrals_map, idx1, values1, n_integrals1)
      call map_append(mo_integrals_map_2, idx2, values2, n_integrals2)

      call ao_to_mo(f2,size(f2,1),d,size(d,1))
      n_integrals1 = 0
      n_integrals2 = 0
      do l=1,mo_num
        do k=1,mo_num
          call ao_two_e_integral_complex_map_idx_sign(list_core_inact_act(i),k,list_core_inact_act(j),l,use_map1,idx_tmp,sign_tmp)
          if (use_map1) then
            n_integrals1+=1
            values1(n_integrals1) = dble(d(k,l))
            idx1(n_integrals1) = idx_tmp
            if (sign_tmp /= 0.d0) then ! should always be true, but might change in the future
              n_integrals1+=1
              values1(n_integrals1) = sign_tmp*dimag(d(k,l))
              idx1(n_integrals1) = idx_tmp+1
            endif
          else
            n_integrals2+=1
            values2(n_integrals2) = dble(d(k,l))
            idx2(n_integrals2) = idx_tmp
            if (sign_tmp /= 0.d0) then
              n_integrals2+=1
              values2(n_integrals2) = sign_tmp*dimag(d(k,l))
              idx2(n_integrals2) = idx_tmp+1
            endif
          endif
        enddo
      enddo
      call map_append(mo_integrals_map, idx1, values1, n_integrals1)
      call map_append(mo_integrals_map_2, idx2, values2, n_integrals2)
    enddo
  enddo
  !$OMP END DO
  deallocate(f,f2,d,idx1,idx2,values1,values2)
  
  !$OMP END PARALLEL

  deallocate(T,T2)
  
  
  call map_sort(mo_integrals_map)
  call map_unique(mo_integrals_map)
  call map_shrink(mo_integrals_map,real(mo_integrals_threshold,integral_kind))

  call map_sort(mo_integrals_map_2)
  call map_unique(mo_integrals_map_2)
  call map_shrink(mo_integrals_map_2,real(mo_integrals_threshold,integral_kind))
  
end

subroutine four_idx_novvvv2_complex
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
    call add_integrals_to_map_complex(mask_ijkl)

    print*, '<ii|vv>'
    do i = 1,N_int
      mask_ijkl(i,1) = core_inact_act_bitmask_4(i,1)
      mask_ijkl(i,2) = core_inact_act_bitmask_4(i,1)
      mask_ijkl(i,3) = virt_bitmask(i,1)
      mask_ijkl(i,4) = virt_bitmask(i,1)
    enddo
    call add_integrals_to_map_complex(mask_ijkl)

end
