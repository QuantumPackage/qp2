BEGIN_PROVIDER [complex*16, chol_mo_integrals_complex, (mo_num_per_kpt,mo_num_per_kpt,chol_num_max,kpt_num,unique_kpt_num)]
  implicit none
  BEGIN_DOC
  !  CD MO integrals
  END_DOC
  integer :: i,j,k,l

  if (read_chol_mo_integrals) then
    call ezfio_get_mo_two_e_ints_chol_mo_integrals_complex(chol_mo_integrals_complex)
    print *,  'CD MO integrals read from disk'
  else
    call chol_mo_from_chol_ao(chol_mo_integrals_complex,chol_ao_integrals_complex,mo_num_per_kpt,ao_num_per_kpt, &
                              chol_num_max,kpt_num,unique_kpt_num)
  endif

  if (write_chol_mo_integrals) then
    call ezfio_set_mo_two_e_ints_chol_mo_integrals_complex(chol_mo_integrals_complex)
    print *,  'CD MO integrals written to disk'
  endif

END_PROVIDER

subroutine mo_map_fill_from_chol_dot
  use map_module
  implicit none
  BEGIN_DOC
  ! TODO: verify correct indexing and conj.transp.
  ! fill mo bielec integral map using 3-index cd integrals
  END_DOC

  integer :: i,k,j,l,mu
  integer :: ki,kk,kj,kl
  integer :: ii,ik,ij,il
  integer :: kikk2,kjkl2,jl2,ik2
  integer :: i_mo,j_mo,i_cd
  integer :: kQ, Q_idx

  complex*16,allocatable :: ints_ik(:,:,:), ints_jl(:,:,:)

  complex*16 :: integral,mjl,mik
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
  !complex*16, external :: zdotc
  complex*16, external :: zdotu

  mo_num_kpt_2 = mo_num_per_kpt * mo_num_per_kpt

  size_buffer = min(mo_num_per_kpt*mo_num_per_kpt*mo_num_per_kpt,16000000)
  print*, 'Providing the mo_bielec integrals from 3-index CD integrals'
  call write_time(6)
!  call ezfio_set_integrals_bielec_disk_access_mo_integrals('Write')
!  TOUCH read_mo_integrals read_ao_integrals write_mo_integrals write_ao_integrals 

  call wall_time(wall_1)
  call cpu_time(cpu_1)

  allocate( ints_jl(chol_num_max,mo_num_per_kpt,mo_num_per_kpt))
  allocate( ints_ik(chol_num_max,mo_num_per_kpt,mo_num_per_kpt))

  wall_0 = wall_1
  do kQ = 1, kpt_num
    Q_idx = kpt_sparse_map(kQ)
    do kl = 1, kpt_num
      kj = qktok2(kQ,kl)
      assert(kQ == qktok2(kj,kl))
      if (kj>kl) cycle
      call idx2_tri_int(kj,kl,kjkl2)
      ints_jl = 0.d0
      if (Q_idx > 0) then
        do i_mo=1,mo_num_per_kpt
          do j_mo=1,mo_num_per_kpt
            do i_cd=1,chol_num(abs(Q_idx))
              !ints_jl(i_cd,i_mo,j_mo) = chol_mo_integrals_complex(i_mo,j_mo,i_cd,kl,Q_idx)
              ints_jl(i_cd,i_mo,j_mo) = dconjg(chol_mo_integrals_complex(i_mo,j_mo,i_cd,kl,Q_idx))
            enddo
          enddo
        enddo
      else
        do i_mo=1,mo_num_per_kpt
          do j_mo=1,mo_num_per_kpt
            do i_cd=1,chol_num(abs(Q_idx))
              !ints_jl(i_cd,i_mo,j_mo) = dconjg(chol_mo_integrals_complex(j_mo,i_mo,i_cd,kj,-Q_idx))
              ints_jl(i_cd,i_mo,j_mo) = chol_mo_integrals_complex(j_mo,i_mo,i_cd,kj,-Q_idx)
            enddo
          enddo
        enddo
      endif

      do kk=1,kl
        ki = qktok2(minusk(kk),kQ)
        assert(ki == kconserv(kl,kk,kj))
        if (ki>kl) cycle
        call idx2_tri_int(ki,kk,kikk2)
        ints_ik = 0.d0
        if (Q_idx > 0) then
          do i_mo=1,mo_num_per_kpt
            do j_mo=1,mo_num_per_kpt
              do i_cd=1,chol_num(abs(Q_idx))
                ints_ik(i_cd,i_mo,j_mo) = chol_mo_integrals_complex(i_mo,j_mo,i_cd,ki,Q_idx)
              enddo
            enddo
          enddo
!          ints_ik = conjg(reshape(df_mo_integral_array(:,:,:,kikk2),(/mo_num_per_kpt,mo_num_per_kpt,df_num/),order=(/2,1,3/)))
        else
          do i_mo=1,mo_num_per_kpt
            do j_mo=1,mo_num_per_kpt
              do i_cd=1,chol_num(abs(Q_idx))
                ints_ik(i_cd,i_mo,j_mo) = dconjg(chol_mo_integrals_complex(j_mo,i_mo,i_cd,kk,-Q_idx))
              enddo
            enddo
          enddo
        endif

        !$OMP PARALLEL PRIVATE(i,k,j,l,ii,ik,ij,il,jl2,ik2, &
            !$OMP  mu, mik, mjl, &
            !$OMP  n_integrals_1, buffer_i_1, buffer_values_1, &
            !$OMP  n_integrals_2, buffer_i_2, buffer_values_2, &
            !$OMP  idx_tmp, tmp_re, tmp_im, integral,sign,use_map1) &
            !$OMP  DEFAULT(NONE)  &
            !$OMP  SHARED(size_buffer, kpt_num, mo_num_per_kpt, mo_num_kpt_2, &
            !$OMP  kl,kj,kjkl2,ints_jl, &
            !$OMP  ki,kk,kikk2,ints_ik, &
            !$OMP  kQ, Q_idx, chol_num, &
            !$OMP  kconserv, chol_mo_integrals_complex, mo_integrals_threshold, &
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
                !integral = zdotc(df_num,ints_jl(1,ij,il),1,ints_ik(1,ii,ik),1)
                !integral = zdotu(chol_num(abs(Q_idx)),ints_jl(1,ij,il),1,ints_ik(1,ii,ik),1)
                integral = zdotu(chol_num(abs(Q_idx)),ints_jl(1,il,ij),1,ints_ik(1,ii,ik),1)
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
                    !call insert_into_mo_integrals_map(n_integrals_1,buffer_i_1,buffer_values_1,mo_integrals_threshold)
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
                    !call insert_into_mo_integrals_map_2(n_integrals_2,buffer_i_2,buffer_values_2,mo_integrals_threshold)
                    n_integrals_2 = 0
                  endif
                endif

              enddo !ii
            enddo !ik
          enddo !ij
        enddo !il
        !$OMP END DO NOWAIT

        if (n_integrals_1 > 0) then
          call map_append(mo_integrals_map, buffer_i_1, buffer_values_1, n_integrals_1)
          !call insert_into_mo_integrals_map(n_integrals_1,buffer_i_1,buffer_values_1,mo_integrals_threshold)
        endif
        if (n_integrals_2 > 0) then
          call map_append(mo_integrals_map_2, buffer_i_2, buffer_values_2, n_integrals_2)
          !call insert_into_mo_integrals_map_2(n_integrals_2,buffer_i_2,buffer_values_2,mo_integrals_threshold)
        endif
        deallocate( &
          buffer_i_1, &
          buffer_i_2, &
          buffer_values_1, &
          buffer_values_2 &
          )
        !$OMP END PARALLEL
      enddo !kk
    enddo !kl
    call wall_time(wall_2)
    if (wall_2 - wall_0 > 1.d0) then
      wall_0 = wall_2
      print*, 100.*float(kQ)/float(kpt_num), '% in ',               &
          wall_2-wall_1,'s',map_mb(mo_integrals_map),'+',map_mb(mo_integrals_map_2),'MB'
    endif

  enddo !kQ
  deallocate( ints_jl,ints_ik )

  call map_sort(mo_integrals_map)
  call map_unique(mo_integrals_map)
  call map_sort(mo_integrals_map_2)
  call map_unique(mo_integrals_map_2)
  !call map_merge(mo_integrals_map)
  !call map_merge(mo_integrals_map_2)

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
  
end subroutine mo_map_fill_from_chol_dot


subroutine chol_mo_from_chol_ao(cd_mo,cd_ao,n_mo,n_ao,n_cd,n_k,n_unique_k)
  use map_module
  implicit none
  BEGIN_DOC
  ! create 3-idx mo ints from 3-idx ao ints
  END_DOC
  integer,intent(in) :: n_mo,n_ao,n_cd,n_k,n_unique_k
  complex*16,intent(out) :: cd_mo(n_mo,n_mo,n_cd,n_k,n_unique_k)
  complex*16,intent(in)  :: cd_ao(n_ao,n_ao,n_cd,n_k,n_unique_k)
  integer :: ki,kk,mu,kQ,Q_idx
  complex*16,allocatable :: coef_i(:,:), coef_k(:,:), ints_ik(:,:), ints_tmp(:,:)
  double precision :: wall_1,wall_2,cpu_1,cpu_2

  print*,'providing 3-index CD MO integrals from 3-index CD AO integrals'
  
  cd_mo = 0.d0

  call wall_time(wall_1)
  call cpu_time(cpu_1)
  allocate( &
              coef_i(n_ao,n_mo),&
              coef_k(n_ao,n_mo),&
             ints_ik(n_ao,n_ao),&
            ints_tmp(n_mo,n_ao)&
          )

  do ki=1, kpt_num
    coef_i = mo_coef_complex_kpts(:,:,ki)
    do kk=1, kpt_num
      coef_k = mo_coef_complex_kpts(:,:,kk)
      kQ = qktok2(kk,ki)
      Q_idx = kpt_sparse_map(kQ)
      if (Q_idx < 0) cycle

      do mu=1, chol_num(abs(Q_idx))
        ints_ik = cd_ao(:,:,mu,ki,Q_idx)
        call zgemm('C','N',n_mo,n_ao,n_ao, &
              (1.d0,0.d0), coef_i, n_ao, &
              ints_ik, n_ao, &
              (0.d0,0.d0), ints_tmp, n_mo)

        call zgemm('N','N',n_mo,n_mo,n_ao, &
              (1.d0,0.d0), ints_tmp, n_mo, &
              coef_k, n_ao, &
              (0.d0,0.d0), cd_mo(:,:,mu,ki,Q_idx), n_mo)
      enddo
    enddo
    call wall_time(wall_2)
    print*,100.*float(ki)/kpt_num, '% in ', &
                wall_2-wall_1, 's'
  enddo

  deallocate( &
            coef_i, &
            coef_k, &
            ints_ik, &
            ints_tmp &
          )
  call wall_time(wall_2)
  call cpu_time(cpu_2)
  print*,' 3-idx CD MO provided'
  print*,'  cpu  time:',cpu_2-cpu_1,'s'
  print*,'  wall time:',wall_2-wall_1,'s  ( x ',(cpu_2-cpu_1)/(wall_2-wall_1),')'

end subroutine chol_mo_from_chol_ao
