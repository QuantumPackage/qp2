 BEGIN_PROVIDER [double precision, df_ao_integrals_real, (ao_kpt_num,ao_kpt_num,df_num,kpt_pair_num)]
&BEGIN_PROVIDER [double precision, df_ao_integrals_imag, (ao_kpt_num,ao_kpt_num,df_num,kpt_pair_num)]
&BEGIN_PROVIDER [complex*16,    df_ao_integrals_complex, (ao_kpt_num,ao_kpt_num,df_num,kpt_pair_num)]
  implicit none
  BEGIN_DOC
  !  df AO integrals
  END_DOC
  integer :: i,j,k,l

  if (read_df_ao_integrals) then
    df_ao_integrals_real = 0.d0
    df_ao_integrals_imag = 0.d0
    call ezfio_get_ao_two_e_ints_df_ao_integrals_real(df_ao_integrals_real)
    call ezfio_get_ao_two_e_ints_df_ao_integrals_imag(df_ao_integrals_imag)
    print *,  'df AO integrals read from disk'
    do l=1,kpt_pair_num
      do k=1,df_num
        do j=1,ao_kpt_num
          do i=1,ao_kpt_num
            df_ao_integrals_complex(i,j,k,l) = dcmplx(df_ao_integrals_real(i,j,k,l), &
                                                      df_ao_integrals_imag(i,j,k,l))
          enddo
        enddo
      enddo
    enddo
  else
    print*,'df ao integrals must be provided',irp_here
    stop -1
  endif

  if (write_df_ao_integrals) then
    ! this probably shouldn't happen
    do l=1,kpt_pair_num
      do k=1,df_num
        do j=1,ao_kpt_num
          do i=1,ao_kpt_num
            df_ao_integrals_real(i,j,k,l) = dble(df_ao_integrals_complex(i,j,k,l))
            df_ao_integrals_imag(i,j,k,l) = dimag(df_ao_integrals_complex(i,j,k,l))
          enddo
        enddo
      enddo
    enddo
    call ezfio_set_ao_two_e_ints_df_ao_integrals_real(df_ao_integrals_real)
    call ezfio_set_ao_two_e_ints_df_ao_integrals_imag(df_ao_integrals_imag)
    print *,  'df AO integrals written to disk'
  endif

END_PROVIDER

