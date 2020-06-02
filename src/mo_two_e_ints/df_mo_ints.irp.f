BEGIN_PROVIDER [complex*16, df_mo_integrals_complex, (mo_num_per_kpt,mo_num_per_kpt,df_num,kpt_pair_num)]
  implicit none
  BEGIN_DOC
  !  df MO integrals
  END_DOC
  integer :: i,j,k,l

  if (read_df_mo_integrals) then
    call ezfio_get_mo_two_e_ints_df_mo_integrals_complex(df_mo_integrals_complex)
    print *,  'df MO integrals read from disk'
  else
    call df_mo_from_df_ao(df_mo_integrals_complex,df_ao_integrals_complex,mo_num_per_kpt,ao_num_per_kpt,df_num,kpt_pair_num)
  endif

  if (write_df_mo_integrals) then
    call ezfio_set_mo_two_e_ints_df_mo_integrals_complex(df_mo_integrals_complex)
    print *,  'df MO integrals written to disk'
  endif

END_PROVIDER

subroutine mo_map_fill_from_df_single
  use map_module
  implicit none
  BEGIN_DOC
  ! fill mo bielec integral map using 3-index df integrals
  END_DOC

  integer :: i,k,j,l
  integer :: ki,kk,kj,kl
  integer :: ii,ik,ij,il
  integer :: kikk2,kjkl2,jl2,ik2
  integer :: i_mo,j_mo,i_df

  complex*16,allocatable :: ints_ik(:,:,:), ints_jl(:,:,:)

  complex*16 :: integral
  integer                        :: n_integrals_1, n_integrals_2
  integer                        :: size_buffer
  integer(key_kind),allocatable  :: buffer_i_1(:), buffer_i_2(:)
  real(integral_kind),allocatable :: buffer_values_1(:), buffer_values_2(:)
  double precision               :: tmp_re,tmp_im
  integer                        :: mo_num_kpt_2

  double precision               :: cpu_1, cpu_2, wall_1, wall_2, wall_0
  double precision               :: map_mb

  logical :: use_map1
  integer(key_kind) :: idx_tmp
  double precision :: sign

  mo_num_kpt_2 = mo_num_per_kpt * mo_num_per_kpt

  size_buffer = min(mo_num_per_kpt*mo_num_per_kpt*mo_num_per_kpt,16000000)
  print*, 'Providing the mo_bielec integrals from 3-index df integrals'
  call write_time(6)
!  call ezfio_set_integrals_bielec_disk_access_mo_integrals('Write')
!  TOUCH read_mo_integrals read_ao_integrals write_mo_integrals write_ao_integrals 

  call wall_time(wall_1)
  call cpu_time(cpu_1)

  allocate( ints_jl(mo_num_per_kpt,mo_num_per_kpt,df_num))
  allocate( ints_ik(mo_num_per_kpt,mo_num_per_kpt,df_num))

  wall_0 = wall_1
  do kl=1, kpt_num
    do kj=1, kl
      call idx2_tri_int(kj,kl,kjkl2)
      if (kj < kl) then
        do i_mo=1,mo_num_per_kpt
          do j_mo=1,mo_num_per_kpt
            do i_df=1,df_num
              ints_jl(i_mo,j_mo,i_df) = dconjg(df_mo_integrals_complex(j_mo,i_mo,i_df,kjkl2))
            enddo
          enddo
        enddo
      else
        ints_jl = df_mo_integrals_complex(:,:,:,kjkl2)
      endif

      do kk=1,kl
        ki=kconserv(kl,kk,kj)
        if (ki>kl) cycle
        call idx2_tri_int(ki,kk,kikk2)
        if (ki < kk) then
          do i_mo=1,mo_num_per_kpt
            do j_mo=1,mo_num_per_kpt
              do i_df=1,df_num
                ints_ik(i_mo,j_mo,i_df) = dconjg(df_mo_integrals_complex(j_mo,i_mo,i_df,kikk2))
              enddo
            enddo
          enddo
!          ints_ik = conjg(reshape(df_mo_integral_array(:,:,:,kikk2),(/mo_num_per_kpt,mo_num_per_kpt,df_num/),order=(/2,1,3/)))
        else
          ints_ik = df_mo_integrals_complex(:,:,:,kikk2)
        endif

        !$OMP PARALLEL PRIVATE(i,k,j,l,ii,ik,ij,il,jl2,ik2, &
            !$OMP  n_integrals_1, buffer_i_1, buffer_values_1, &
            !$OMP  n_integrals_2, buffer_i_2, buffer_values_2, &
            !$OMP  idx_tmp, tmp_re, tmp_im, integral,sign,use_map1) &
            !$OMP  DEFAULT(NONE)  &
            !$OMP  SHARED(size_buffer, kpt_num, df_num, mo_num_per_kpt, mo_num_kpt_2, &
            !$OMP  kl,kj,kjkl2,ints_jl, &
            !$OMP  ki,kk,kikk2,ints_ik, &
            !$OMP  kconserv, df_mo_integrals_complex, mo_integrals_threshold, &
            !$OMP  mo_integrals_map, mo_integrals_map_2)

        allocate( &
          buffer_i_1(size_buffer), &
          buffer_i_2(size_buffer), &
          buffer_values_1(size_buffer), &
          buffer_values_2(size_buffer) &
        )
      
        n_integrals_1=0
        n_integrals_2=0
        !$OMP DO SCHEDULE(guided)
        do mu=1,df_num
          do il=1,mo_num_per_kpt
            l=il+(kl-1)*mo_num_per_kpt
            do ij=1,mo_num_per_kpt
              j=ij+(kj-1)*mo_num_per_kpt
              if (j>l) exit
              call idx2_tri_int(j,l,jl2)
              mjl = ints_jl(ij,il,mu)
              if (mjl.eq.(0.d0,0.d0)) cycle
              do ik=1,mo_num_per_kpt
                k=ik+(kk-1)*mo_num_per_kpt
                if (k>l) exit
                do ii=1,mo_num_per_kpt
                  i=ii+(ki-1)*mo_num_per_kpt
                  if ((j==l) .and. (i>k)) exit
                  call idx2_tri_int(i,k,ik2)
                  if (ik2 > jl2) exit
                  mik = ints_ik(ii,ik,mu)
                  integral = mik * dconjg(mjl)
!                  print*,i,k,j,l,real(integral),imag(integral)
                  if (cdabs(integral) < mo_integrals_threshold) then
                    cycle
                  endif
                  call ao_two_e_integral_complex_map_idx_sign(i,j,k,l,use_map1,idx_tmp,sign)
                  tmp_re = dble(integral)
                  tmp_im = dimag(integral)
                  if (use_map1) then
                    n_integrals_1 += 1
                    buffer_i_1(n_integrals_1)=idx_tmp
                    buffer_values_1(n_integrals_1)=tmp_re
                    if (sign.ne.0.d0) then
                      n_integrals_1 += 1
                      buffer_i_1(n_integrals_1)=idx_tmp+1
                      buffer_values_1(n_integrals_1)=tmp_im*sign
                    endif 
                    if (n_integrals_1 >= size(buffer_i_1)-1) then
                      !call map_append(mo_integrals_map, buffer_i_1, buffer_values_1, n_integrals_1)
                      call insert_into_ao_integrals_map(n_integrals_1,buffer_i_1,buffer_values_1)
                      n_integrals_1 = 0
                    endif
                  else
                    n_integrals_2 += 1
                    buffer_i_2(n_integrals_2)=idx_tmp
                    buffer_values_2(n_integrals_2)=tmp_re
                    if (sign.ne.0.d0) then
                      n_integrals_2 += 1
                      buffer_i_2(n_integrals_2)=idx_tmp+1
                      buffer_values_2(n_integrals_2)=tmp_im*sign
                    endif
                    if (n_integrals_2 >= size(buffer_i_2)-1) then
                      !call map_append(mo_integrals_map_2, buffer_i_2, buffer_values_2, n_integrals_2)
                      call insert_into_ao_integrals_map_2(n_integrals_2,buffer_i_2,buffer_values_2)
                      n_integrals_2 = 0
                    endif
                  endif

                enddo !ii
              enddo !ik
            enddo !ij
          enddo !il
        enddo !mu
        !$OMP END DO NOWAIT

        if (n_integrals_1 > 0) then
          !call map_append(mo_integrals_map, buffer_i_1, buffer_values_1, n_integrals_1)
          call insert_into_ao_integrals_map(n_integrals_1,buffer_i_1,buffer_values_1)
        endif
        if (n_integrals_2 > 0) then
          !call map_append(mo_integrals_map_2, buffer_i_2, buffer_values_2, n_integrals_2)
          call insert_into_ao_integrals_map_2(n_integrals_2,buffer_i_2,buffer_values_2)
        endif
        deallocate( &
          buffer_i_1, &
          buffer_i_2, &
          buffer_values_1, &
          buffer_values_2 &
          )
        !$OMP END PARALLEL
      enddo !kk
    enddo !kj
    call wall_time(wall_2)
    if (wall_2 - wall_0 > 1.d0) then
      wall_0 = wall_2
      print*, 100.*float(kl)/float(kpt_num), '% in ',               &
          wall_2-wall_1,'s',map_mb(mo_integrals_map),'+',map_mb(mo_integrals_map_2),'MB'
    endif

  enddo !kl
  deallocate( ints_jl,ints_ik )

  !call map_sort(mo_integrals_map)
  !call map_unique(mo_integrals_map)
  !call map_sort(mo_integrals_map_2)
  !call map_unique(mo_integrals_map_2)
  call map_merge(mo_integrals_map)
  call map_merge(mo_integrals_map_2)
  !!call map_save_to_disk(trim(ezfio_filename)//'/work/mo_ints_complex_1',mo_integrals_map)
  !!call map_save_to_disk(trim(ezfio_filename)//'/work/mo_ints_complex_2',mo_integrals_map_2)
  !!call ezfio_set_mo_two_e_ints_io_mo_two_e_integrals('Read')
  
  call wall_time(wall_2)
  call cpu_time(cpu_2)

  integer*8                      :: get_mo_map_size, mo_map_size
  mo_map_size = get_mo_map_size()
  
  print*,'MO integrals provided:'
  print*,' Size of MO map           ', map_mb(mo_integrals_map),'+',map_mb(mo_integrals_map_2),'MB'
  print*,' Number of MO integrals: ',  mo_map_size
  print*,' cpu  time :',cpu_2 - cpu_1, 's'
  print*,' wall time :',wall_2 - wall_1, 's  ( x ', (cpu_2-cpu_1)/(wall_2-wall_1), ')'
  
end subroutine mo_map_fill_from_df_single

subroutine mo_map_fill_from_df
  use map_module
  implicit none
  BEGIN_DOC
  ! fill mo bielec integral map using 3-index df integrals
  END_DOC

  integer :: i,k,j,l
  integer :: ki,kk,kj,kl
  integer :: ii,ik,ij,il
  integer :: kikk2,kjkl2,jl2,ik2
  integer :: i_mo,j_mo,i_df

  complex*16,allocatable :: ints_ik(:,:,:), ints_jl(:,:,:), ints_ikjl(:,:,:,:)

  complex*16 :: integral
  integer                        :: n_integrals_1, n_integrals_2
  integer                        :: size_buffer
  integer(key_kind),allocatable  :: buffer_i_1(:), buffer_i_2(:)
  real(integral_kind),allocatable :: buffer_values_1(:), buffer_values_2(:)
  double precision               :: tmp_re,tmp_im
  integer                        :: mo_num_kpt_2

  double precision               :: cpu_1, cpu_2, wall_1, wall_2, wall_0
  double precision               :: map_mb

  logical :: use_map1
  integer(keY_kind) :: idx_tmp
  double precision :: sign

  mo_num_kpt_2 = mo_num_per_kpt * mo_num_per_kpt

  size_buffer = min(mo_num_per_kpt*mo_num_per_kpt*mo_num_per_kpt,16000000)
  print*, 'Providing the mo_bielec integrals from 3-index df integrals'
  call write_time(6)
!  call ezfio_set_integrals_bielec_disk_access_mo_integrals('Write')
!  TOUCH read_mo_integrals read_ao_integrals write_mo_integrals write_ao_integrals 

  call wall_time(wall_1)
  call cpu_time(cpu_1)

  allocate( ints_jl(mo_num_per_kpt,mo_num_per_kpt,df_num))

  wall_0 = wall_1
  do kl=1, kpt_num
    do kj=1, kl
      call idx2_tri_int(kj,kl,kjkl2)
      if (kj < kl) then
        do i_mo=1,mo_num_per_kpt
          do j_mo=1,mo_num_per_kpt
            do i_df=1,df_num
              ints_jl(i_mo,j_mo,i_df) = dconjg(df_mo_integrals_complex(j_mo,i_mo,i_df,kjkl2))
            enddo
          enddo
        enddo
      else
        ints_jl = df_mo_integrals_complex(:,:,:,kjkl2)
      endif

  !$OMP PARALLEL PRIVATE(i,k,j,l,ki,kk,ii,ik,ij,il,kikk2,jl2,ik2, &
      !$OMP  ints_ik, ints_ikjl, i_mo, j_mo, i_df, &
      !$OMP  n_integrals_1, buffer_i_1, buffer_values_1, &
      !$OMP  n_integrals_2, buffer_i_2, buffer_values_2, &
      !$OMP  idx_tmp, tmp_re, tmp_im, integral,sign,use_map1) &
      !$OMP  DEFAULT(NONE)  &
      !$OMP  SHARED(size_buffer, kpt_num, df_num, mo_num_per_kpt, mo_num_kpt_2, &
      !$OMP  kl,kj,kjkl2,ints_jl, & 
      !$OMP  kconserv, df_mo_integrals_complex, mo_integrals_threshold, mo_integrals_map, mo_integrals_map_2)
  
  allocate( &
    ints_ik(mo_num_per_kpt,mo_num_per_kpt,df_num), &
    ints_ikjl(mo_num_per_kpt,mo_num_per_kpt,mo_num_per_kpt,mo_num_per_kpt), &
    buffer_i_1(size_buffer), &
    buffer_i_2(size_buffer), &
    buffer_values_1(size_buffer), &
    buffer_values_2(size_buffer) &
  )

  !$OMP DO SCHEDULE(guided)
      do kk=1,kl
        ki=kconserv(kl,kk,kj)
        if (ki>kl) cycle
      !  if ((kl == kj) .and. (ki > kk)) cycle
        call idx2_tri_int(ki,kk,kikk2)
      !  if (kikk2 > kjkl2) cycle
        if (ki < kk) then
          do i_mo=1,mo_num_per_kpt
            do j_mo=1,mo_num_per_kpt
              do i_df=1,df_num
                ints_ik(i_mo,j_mo,i_df) = dconjg(df_mo_integrals_complex(j_mo,i_mo,i_df,kikk2))
              enddo
            enddo
          enddo
!          ints_ik = conjg(reshape(df_mo_integral_array(:,:,:,kikk2),(/mo_num_per_kpt,mo_num_per_kpt,df_num/),order=(/2,1,3/)))
        else
          ints_ik = df_mo_integrals_complex(:,:,:,kikk2)
        endif

        call zgemm('N','T', mo_num_kpt_2, mo_num_kpt_2, df_num, &
               (1.d0,0.d0), ints_ik, mo_num_kpt_2, &
               ints_jl, mo_num_kpt_2, &
               (0.d0,0.d0), ints_ikjl, mo_num_kpt_2)

        n_integrals_1=0
        n_integrals_2=0
        do il=1,mo_num_per_kpt
          l=il+(kl-1)*mo_num_per_kpt
          do ij=1,mo_num_per_kpt
            j=ij+(kj-1)*mo_num_per_kpt
            if (j>l) exit
            call idx2_tri_int(j,l,jl2)
            do ik=1,mo_num_per_kpt
              k=ik+(kk-1)*mo_num_per_kpt
              if (k>l) exit
              do ii=1,mo_num_per_kpt
                i=ii+(ki-1)*mo_num_per_kpt
                if ((j==l) .and. (i>k)) exit
                call idx2_tri_int(i,k,ik2)
                if (ik2 > jl2) exit
                integral = ints_ikjl(ii,ik,ij,il)
!                print*,i,k,j,l,real(integral),imag(integral)
                if (cdabs(integral) < mo_integrals_threshold) then
                  cycle
                endif
                call ao_two_e_integral_complex_map_idx_sign(i,j,k,l,use_map1,idx_tmp,sign)
                tmp_re = dble(integral)
                tmp_im = dimag(integral)
                if (use_map1) then
                  n_integrals_1 += 1
                  buffer_i_1(n_integrals_1)=idx_tmp
                  buffer_values_1(n_integrals_1)=tmp_re
                  if (sign.ne.0.d0) then 
                    n_integrals_1 += 1
                    buffer_i_1(n_integrals_1)=idx_tmp+1
                    buffer_values_1(n_integrals_1)=tmp_im*sign
                  endif 
                  if (n_integrals_1 >= size(buffer_i_1)-1) then
                    call map_append(mo_integrals_map, buffer_i_1, buffer_values_1, n_integrals_1)
                    !call insert_into_ao_integrals_map(n_integrals_1,buffer_i_1,buffer_values_1)
                    n_integrals_1 = 0
                  endif
                else
                  n_integrals_2 += 1
                  buffer_i_2(n_integrals_2)=idx_tmp
                  buffer_values_2(n_integrals_2)=tmp_re
                  if (sign.ne.0.d0) then
                    n_integrals_2 += 1
                    buffer_i_2(n_integrals_2)=idx_tmp+1
                    buffer_values_2(n_integrals_2)=tmp_im*sign
                  endif 
                  if (n_integrals_2 >= size(buffer_i_2)-1) then
                    call map_append(mo_integrals_map_2, buffer_i_2, buffer_values_2, n_integrals_2)
                    !call insert_into_ao_integrals_map_2(n_integrals_2,buffer_i_2,buffer_values_2)
                    n_integrals_2 = 0
                  endif
                endif

              enddo !ii
            enddo !ik
          enddo !ij
        enddo !il

        if (n_integrals_1 > 0) then
          call map_append(mo_integrals_map, buffer_i_1, buffer_values_1, n_integrals_1)
          !call insert_into_ao_integrals_map(n_integrals_1,buffer_i_1,buffer_values_1)
        endif
        if (n_integrals_2 > 0) then
          call map_append(mo_integrals_map_2, buffer_i_2, buffer_values_2, n_integrals_2)
          !call insert_into_ao_integrals_map_2(n_integrals_2,buffer_i_2,buffer_values_2)
        endif
      enddo !kk
  !$OMP END DO NOWAIT
  deallocate( &
    ints_ik, &
    ints_ikjl, &
    buffer_i_1, &
    buffer_i_2, &
    buffer_values_1, &
    buffer_values_2 &
    )
  !$OMP END PARALLEL
    enddo !kj
    call wall_time(wall_2)
    if (wall_2 - wall_0 > 1.d0) then
      wall_0 = wall_2
      print*, 100.*float(kl)/float(kpt_num), '% in ',               &
          wall_2-wall_1,'s',map_mb(mo_integrals_map),'+',map_mb(mo_integrals_map_2),'MB'
    endif

  enddo !kl
  deallocate( ints_jl ) 

  call map_sort(mo_integrals_map)
  call map_unique(mo_integrals_map)
  call map_sort(mo_integrals_map_2)
  call map_unique(mo_integrals_map_2)
  !call map_save_to_disk(trim(ezfio_filename)//'/work/mo_ints_complex_1',mo_integrals_map)
  !call map_save_to_disk(trim(ezfio_filename)//'/work/mo_ints_complex_2',mo_integrals_map_2)
  !call ezfio_set_mo_two_e_ints_io_mo_two_e_integrals('Read')
  
  call wall_time(wall_2)
  call cpu_time(cpu_2)

  integer*8                      :: get_mo_map_size, mo_map_size
  mo_map_size = get_mo_map_size()
  
  print*,'MO integrals provided:'
  print*,' Size of MO map           ', map_mb(mo_integrals_map),'+',map_mb(mo_integrals_map_2),'MB'
  print*,' Number of MO integrals: ',  mo_map_size
  print*,' cpu  time :',cpu_2 - cpu_1, 's'
  print*,' wall time :',wall_2 - wall_1, 's  ( x ', (cpu_2-cpu_1)/(wall_2-wall_1), ')'
  
end subroutine mo_map_fill_from_df

subroutine df_mo_from_df_ao(df_mo,df_ao,n_mo,n_ao,n_df,n_k_pairs)
  use map_module
  implicit none
  BEGIN_DOC
  ! create 3-idx mo ints from 3-idx ao ints
  END_DOC
  integer,intent(in) :: n_mo,n_ao,n_df,n_k_pairs
  complex*16,intent(out) :: df_mo(n_mo,n_mo,n_df,n_k_pairs)
  complex*16,intent(in)  :: df_ao(n_ao,n_ao,n_df,n_k_pairs)
  integer :: kl,kj,kjkl2,mu,p,q
  complex*16,allocatable :: coef_l(:,:), coef_j(:,:), ints_jl(:,:), ints_tmp(:,:)
  double precision :: wall_1,wall_2,cpu_1,cpu_2

  print*,'providing 3-index MO integrals from 3-index AO integrals'

  call wall_time(wall_1)
  call cpu_time(cpu_1)
  allocate( &
              coef_l(n_ao,n_mo),&
              coef_j(n_ao,n_mo),&
             ints_jl(n_ao,n_ao),&
            ints_tmp(n_mo,n_ao)&
          )

  do kl=1, kpt_num
    coef_l = mo_coef_complex_kpts(:,:,kl)
    do kj=1, kl
      coef_j = mo_coef_complex_kpts(:,:,kj)
      kjkl2 = kj+shiftr(kl*kl-kl,1)
      do mu=1, df_num
        ints_jl = df_ao(:,:,mu,kjkl2)
        call zgemm('C','N',n_mo,n_ao,n_ao, &
              (1.d0,0.d0), coef_l, n_ao, &
              ints_jl, n_ao, &
              (0.d0,0.d0), ints_tmp, n_mo)

        call zgemm('N','N',n_mo,n_mo,n_ao, &
              (1.d0,0.d0), ints_tmp, n_mo, &
              coef_j, n_ao, &
              (0.d0,0.d0), df_mo(:,:,mu,kjkl2), n_mo)
      enddo
    enddo
    call wall_time(wall_2)
    print*,100.*float(kl*(kl+1))/(2.*n_k_pairs), '% in ', &
                wall_2-wall_1, 's'
  enddo

  deallocate( &
            coef_l, &
            coef_j, &
            ints_jl, &
            ints_tmp &
          )
  call wall_time(wall_2)
  call cpu_time(cpu_2)
  print*,' 3-idx MO provided'
  print*,'  cpu  time:',cpu_2-cpu_1,'s'
  print*,'  wall time:',wall_2-wall_1,'s  ( x ',(cpu_2-cpu_1)/(wall_2-wall_1),')'

end subroutine df_mo_from_df_ao
