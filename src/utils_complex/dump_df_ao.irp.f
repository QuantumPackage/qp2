program dump_df_ao
  call run
end

subroutine run
  use map_module
  implicit none

  integer ::i,j,k,mu
  complex*16 :: integral

  provide df_ao_integrals_complex
  do k=1,kpt_pair_num
    do mu=1,df_num
      do i=1,ao_num_per_kpt
        do j=1,ao_num_per_kpt
          integral = df_ao_integrals_complex(i,j,mu,k)
          if (cdabs(integral).gt.1.d-12) then
            print'(4(I4),4(E15.7))',i,j,mu,k,integral,dble(integral),dimag(integral)
          endif
        enddo
      enddo
    enddo
  enddo

end
