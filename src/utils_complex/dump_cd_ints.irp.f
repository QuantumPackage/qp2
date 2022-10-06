program dump_cd_ksym
  call run
end

subroutine run
  use map_module
  implicit none

  integer ::q,k,n,i,j
  double precision :: vr, vi
  complex*16 :: v
  print*,"chol_ao_integrals_complex q,k,n,i,j"
  provide chol_ao_integrals_complex
  do q = 1, unique_kpt_num
    do k = 1, kpt_num
      do n = 1, chol_num_max
        do i = 1, ao_num_per_kpt
          do j = 1, ao_num_per_kpt
            v = chol_ao_integrals_complex(i,j,n,k,q)
            vr = dble(v)
            vi = dimag(v)
            print '(5(I6,X),2(E25.15,X))', q, k, n, i, j, vr, vi
          enddo
        enddo
      enddo
    enddo
  enddo

end
