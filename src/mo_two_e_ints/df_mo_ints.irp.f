 BEGIN_PROVIDER [double precision, df_mo_integrals_real, (mo_num_per_kpt,mo_num_per_kpt,df_num,kpt_pair_num)]
&BEGIN_PROVIDER [double precision, df_mo_integrals_imag, (mo_num_per_kpt,mo_num_per_kpt,df_num,kpt_pair_num)]
&BEGIN_PROVIDER [complex*16,    df_mo_integrals_complex, (mo_num_per_kpt,mo_num_per_kpt,df_num,kpt_pair_num)]
  implicit none
  BEGIN_DOC
  !  df AO integrals
  END_DOC
  integer :: i,j,k,l

  if (read_df_mo_integrals) then
    df_mo_integrals_real = 0.d0
    df_mo_integrals_imag = 0.d0
    call ezfio_get_mo_two_e_ints_df_mo_integrals_real(df_mo_integrals_real)
    call ezfio_get_mo_two_e_ints_df_mo_integrals_imag(df_mo_integrals_imag)
    print *,  'df AO integrals read from disk'
    do l=1,kpt_pair_num
      do k=1,df_num
        do j=1,mo_num_per_kpt
          do i=1,mo_num_per_kpt
            df_mo_integrals_complex(i,j,k,l) = dcmplx(df_mo_integrals_real(i,j,k,l), &
                                                      df_mo_integrals_imag(i,j,k,l))
          enddo
        enddo
      enddo
    enddo
  else
    call df_mo_from_df_ao(df_mo_integrals_complex,df_ao_integrals_complex,mo_num_per_kpt,ao_num_per_kpt,df_num,kpt_pair_num)
  endif

  if (write_df_mo_integrals) then
    do l=1,kpt_pair_num
      do k=1,df_num
        do j=1,mo_num_per_kpt
          do i=1,mo_num_per_kpt
            df_mo_integrals_real(i,j,k,l) = dble(df_mo_integrals_complex(i,j,k,l))
            df_mo_integrals_imag(i,j,k,l) = dimag(df_mo_integrals_complex(i,j,k,l))
          enddo
        enddo
      enddo
    enddo
    call ezfio_set_mo_two_e_ints_df_mo_integrals_real(df_mo_integrals_real)
    call ezfio_set_mo_two_e_ints_df_mo_integrals_imag(df_mo_integrals_imag)
    print *,  'df AO integrals written to disk'
  endif

END_PROVIDER

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
    coef_l = mo_coef_kpts(:,:,kl)
    do kj=1, kl
      coef_j = mo_coef_kpts(:,:,kj)
      kjkl2 = kj+shiftr(kl*kl-kl,1)
      do mu=1, df_num
        ints_jl = df_ao(:,:,mu,kjkl2)
        call zgemm('C','N',n_mo,n_ao,n_ao, &
              (1.d0,0.d0), coef_j, n_ao, &
              ints_jl, n_ao, &
              (0.d0,0.d0), ints_tmp, n_mo)

        call zgemm('N','N',n_mo,n_mo,n_ao, &
              (1.d0,0.d0), ints_tmp, n_mo, &
              coef_l, n_ao, &
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
