!  BEGIN_PROVIDER [double precision, df_ao_integrals_real, (ao_num_per_kpt,ao_num_per_kpt,df_num,kpt_pair_num)]
! &BEGIN_PROVIDER [double precision, df_ao_integrals_imag, (ao_num_per_kpt,ao_num_per_kpt,df_num,kpt_pair_num)]
! &BEGIN_PROVIDER [complex*16,    df_ao_integrals_complex, (ao_num_per_kpt,ao_num_per_kpt,df_num,kpt_pair_num)]
!   implicit none
!   BEGIN_DOC
!   !  df AO integrals
!   END_DOC
!   integer :: i,j,k,l
! 
!   if (read_df_ao_integrals) then
!     df_ao_integrals_real = 0.d0
!     df_ao_integrals_imag = 0.d0
!     call ezfio_get_ao_two_e_ints_df_ao_integrals_real(df_ao_integrals_real)
!     call ezfio_get_ao_two_e_ints_df_ao_integrals_imag(df_ao_integrals_imag)
!     print *,  'df AO integrals read from disk'
!     do l=1,kpt_pair_num
!       do k=1,df_num
!         do j=1,ao_num_per_kpt
!           do i=1,ao_num_per_kpt
!             df_ao_integrals_complex(i,j,k,l) = dcmplx(df_ao_integrals_real(i,j,k,l), &
!                                                       df_ao_integrals_imag(i,j,k,l))
!           enddo
!         enddo
!       enddo
!     enddo
!   else
!     print*,'df ao integrals must be provided',irp_here
!     stop -1
!   endif
! 
!   if (write_df_ao_integrals) then
!     ! this probably shouldn't happen
!     do l=1,kpt_pair_num
!       do k=1,df_num
!         do j=1,ao_num_per_kpt
!           do i=1,ao_num_per_kpt
!             df_ao_integrals_real(i,j,k,l) = dble(df_ao_integrals_complex(i,j,k,l))
!             df_ao_integrals_imag(i,j,k,l) = dimag(df_ao_integrals_complex(i,j,k,l))
!           enddo
!         enddo
!       enddo
!     enddo
!     call ezfio_set_ao_two_e_ints_df_ao_integrals_real(df_ao_integrals_real)
!     call ezfio_set_ao_two_e_ints_df_ao_integrals_imag(df_ao_integrals_imag)
!     print *,  'df AO integrals written to disk'
!   endif
! 
! END_PROVIDER
 BEGIN_PROVIDER [complex*16, df_ao_integrals_complex, (ao_num_per_kpt,ao_num_per_kpt,df_num,kpt_pair_num)]
  implicit none
  BEGIN_DOC
  !  df AO integrals
  END_DOC
  integer :: i,j,k,l

  if (read_df_ao_integrals) then
    call ezfio_get_ao_two_e_ints_df_ao_integrals_complex(df_ao_integrals_complex)
    print *,  'df AO integrals read from disk'
  else
    print*,'df ao integrals must be provided',irp_here
    stop -1
  endif

  if (write_df_ao_integrals) then
    call ezfio_set_ao_two_e_ints_df_ao_integrals_complex(df_ao_integrals_complex)
    print *,  'df AO integrals written to disk'
  endif

END_PROVIDER


subroutine ao_map_fill_from_df
  use map_module
  implicit none
  BEGIN_DOC
  ! fill ao bielec integral map using 3-index df integrals
  END_DOC

  integer :: i,k,j,l
  integer :: ki,kk,kj,kl
  integer :: ii,ik,ij,il
  integer :: kikk2,kjkl2,jl2,ik2
  integer :: i_ao,j_ao,i_df

  complex*16,allocatable :: ints_ik(:,:,:), ints_jl(:,:,:), ints_ikjl(:,:,:,:)

  complex*16 :: integral
  integer                        :: n_integrals_1, n_integrals_2
  integer                        :: size_buffer
  integer(key_kind),allocatable  :: buffer_i_1(:), buffer_i_2(:)
  real(integral_kind),allocatable :: buffer_values_1(:), buffer_values_2(:)
  double precision               :: tmp_re,tmp_im
  integer                        :: ao_num_kpt_2

  double precision               :: cpu_1, cpu_2, wall_1, wall_2, wall_0
  double precision               :: map_mb

  logical :: use_map1
  integer(keY_kind) :: idx_tmp
  double precision :: sign

  ao_num_kpt_2 = ao_num_per_kpt * ao_num_per_kpt

  size_buffer = min(ao_num_per_kpt*ao_num_per_kpt*ao_num_per_kpt,16000000)
  print*, 'Providing the ao_bielec integrals from 3-index df integrals'
  call write_time(6)
!  call ezfio_set_integrals_bielec_disk_access_mo_integrals('Write')
!  TOUCH read_mo_integrals read_ao_integrals write_mo_integrals write_ao_integrals 

  call wall_time(wall_1)
  call cpu_time(cpu_1)

  allocate( ints_jl(ao_num_per_kpt,ao_num_per_kpt,df_num))

  wall_0 = wall_1
  do kl=1, kpt_num
    do kj=1, kl
      call idx2_tri_int(kj,kl,kjkl2)
      if (kj < kl) then
        do i_ao=1,ao_num_per_kpt
          do j_ao=1,ao_num_per_kpt
            do i_df=1,df_num
              ints_jl(i_ao,j_ao,i_df) = dconjg(df_ao_integrals_complex(j_ao,i_ao,i_df,kjkl2))
            enddo
          enddo
        enddo
      else
        ints_jl = df_ao_integrals_complex(:,:,:,kjkl2)
      endif

  !$OMP PARALLEL PRIVATE(i,k,j,l,ki,kk,ii,ik,ij,il,kikk2,jl2,ik2, &
      !$OMP  ints_ik, ints_ikjl, i_ao, j_ao, i_df, &
      !$OMP  n_integrals_1, buffer_i_1, buffer_values_1, &
      !$OMP  n_integrals_2, buffer_i_2, buffer_values_2, &
      !$OMP  idx_tmp, tmp_re, tmp_im, integral,sign,use_map1) &
      !$OMP  DEFAULT(NONE)  &
      !$OMP  SHARED(size_buffer, kpt_num, df_num, ao_num_per_kpt, ao_num_kpt_2, &
      !$OMP  kl,kj,kjkl2,ints_jl, & 
      !$OMP  kconserv, df_ao_integrals_complex, ao_integrals_threshold, ao_integrals_map, ao_integrals_map_2)
  
  allocate( &
    ints_ik(ao_num_per_kpt,ao_num_per_kpt,df_num), &
    ints_ikjl(ao_num_per_kpt,ao_num_per_kpt,ao_num_per_kpt,ao_num_per_kpt), &
    buffer_i_1(size_buffer), &
    buffer_i_2(size_buffer), &
    buffer_values_1(size_buffer), &
    buffer_values_2(size_buffer) &
  )

  !$OMP DO SCHEDULE(guided)
      do kk=1,kl
        ki=kconserv(kl,kk,kj)
        if ((kl == kj) .and. (ki > kk)) cycle
        call idx2_tri_int(ki,kk,kikk2)
        if (kikk2 > kjkl2) cycle
        if (ki < kk) then
          do i_ao=1,ao_num_per_kpt
            do j_ao=1,ao_num_per_kpt
              do i_df=1,df_num
                ints_ik(i_ao,j_ao,i_df) = dconjg(df_ao_integrals_complex(j_ao,i_ao,i_df,kikk2))
              enddo
            enddo
          enddo
!          ints_ik = conjg(reshape(df_mo_integral_array(:,:,:,kikk2),(/mo_num_per_kpt,mo_num_per_kpt,df_num/),order=(/2,1,3/)))
        else
          ints_ik = df_ao_integrals_complex(:,:,:,kikk2)
        endif

        call zgemm('N','T', ao_num_kpt_2, ao_num_kpt_2, df_num, &
               (1.d0,0.d0), ints_ik, ao_num_kpt_2, &
               ints_jl, ao_num_kpt_2, &
               (0.d0,0.d0), ints_ikjl, ao_num_kpt_2)

        n_integrals_1=0
        n_integrals_2=0
        do il=1,ao_num_per_kpt
          l=il+(kl-1)*ao_num_per_kpt
          do ij=1,ao_num_per_kpt
            j=ij+(kj-1)*ao_num_per_kpt
            if (j>l) exit
            call idx2_tri_int(j,l,jl2)
            do ik=1,ao_num_per_kpt
              k=ik+(kk-1)*ao_num_per_kpt
              if (k>l) exit
              do ii=1,ao_num_per_kpt
                i=ii+(ki-1)*ao_num_per_kpt
                if ((j==l) .and. (i>k)) exit
                call idx2_tri_int(i,k,ik2)
                if (ik2 > jl2) exit
                integral = ints_ikjl(ii,ik,ij,il)
!                print*,i,k,j,l,real(integral),imag(integral)
                if (cdabs(integral) < ao_integrals_threshold) then
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
                    call insert_into_ao_integrals_map_2(n_integrals_2,buffer_i_2,buffer_values_2)
                    n_integrals_2 = 0
                  endif
                endif

              enddo !ii
            enddo !ik
          enddo !ij
        enddo !il

        if (n_integrals_1 > 0) then
          call insert_into_ao_integrals_map(n_integrals_1,buffer_i_1,buffer_values_1)
        endif
        if (n_integrals_2 > 0) then
          call insert_into_ao_integrals_map_2(n_integrals_2,buffer_i_2,buffer_values_2)
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
          wall_2-wall_1,'s',map_mb(ao_integrals_map),'+',map_mb(ao_integrals_map_2),'MB'
    endif

  enddo !kl
  deallocate( ints_jl ) 

  call map_sort(ao_integrals_map)
  call map_unique(ao_integrals_map)
  call map_sort(ao_integrals_map_2)
  call map_unique(ao_integrals_map_2)
  !call map_save_to_disk(trim(ezfio_filename)//'/work/ao_ints_complex_1',ao_integrals_map)
  !call map_save_to_disk(trim(ezfio_filename)//'/work/ao_ints_complex_2',ao_integrals_map_2)
  !call ezfio_set_ao_two_e_ints_io_ao_two_e_integrals('Read')
  
  call wall_time(wall_2)
  call cpu_time(cpu_2)

  integer*8                      :: get_ao_map_size, ao_map_size
  ao_map_size = get_ao_map_size()
  
  print*,'AO integrals provided:'
  print*,' Size of AO map           ', map_mb(ao_integrals_map),'+',map_mb(ao_integrals_map_2),'MB'
  print*,' Number of AO integrals: ',  ao_map_size
  print*,' cpu  time :',cpu_2 - cpu_1, 's'
  print*,' wall time :',wall_2 - wall_1, 's  ( x ', (cpu_2-cpu_1)/(wall_2-wall_1), ')'
  
end subroutine ao_map_fill_from_df

