program dump_ao_2e_from_df
  call run_ao_dump
end

subroutine run_ao_dump
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

  complex*16 :: integral,intmap, get_ao_two_e_integral_complex
  double precision               :: tmp_re,tmp_im
  integer                        :: ao_num_kpt_2

  logical :: use_map1
  integer(keY_kind) :: idx_tmp
  double precision :: sign

  ao_num_kpt_2 = ao_num_per_kpt * ao_num_per_kpt


  allocate( ints_jl(ao_num_per_kpt,ao_num_per_kpt,df_num))

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
  
  allocate( &
    ints_ik(ao_num_per_kpt,ao_num_per_kpt,df_num), &
    ints_ikjl(ao_num_per_kpt,ao_num_per_kpt,ao_num_per_kpt,ao_num_per_kpt) &
  )

      do kk=1,kl
        ki=kconserv(kl,kk,kj)
        if (ki > kl) cycle
        !if ((kl == kj) .and. (ki > kk)) cycle
        call idx2_tri_int(ki,kk,kikk2)
        !if (kikk2 > kjkl2) cycle
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
        print'((A),4(I4))','IJKL',ki,kj,kk,kl
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
                intmap = get_ao_two_e_integral_complex(i,j,k,l,ao_integrals_map,ao_integrals_map_2)
!                print*,i,k,j,l,real(integral),imag(integral)
                if ((cdabs(integral) + cdabs(intmap)) < ao_integrals_threshold) then
                  cycle
                endif
                if (cdabs(integral-intmap) < 1.d-8) then
                  print'(4(I4),4(E15.7))',i,j,k,l,integral,intmap
                else
                  print'(4(I4),4(E15.7),(A))',i,j,k,l,integral,intmap,'***'
                endif
              enddo !ii
            enddo !ik
          enddo !ij
        enddo !il
      enddo !kk
  deallocate( &
    ints_ik, &
    ints_ikjl &
    )
    enddo !kj
  enddo !kl
  deallocate( ints_jl ) 

  
end

