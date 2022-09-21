program dump_cd_ksym
  call run
end

subroutine run
  use map_module
  implicit none

  integer ::i,j,k,l
  integer(key_kind) :: idx
  logical :: use_map1
  double precision :: sign
  do i=1,5
    do j=1,5
      do k=1,5
        do l=1,5
          call ao_two_e_integral_complex_map_idx_sign(i,j,k,l,use_map1,idx,sign)
          print'(4(I4,X),(L6),(I8),(F10.1))',i,j,k,l,use_map1,idx,sign
        enddo
      enddo
    enddo
  enddo

  provide qktok2 minusk kconserv
  print*,'minusk'
  do i=1,kpt_num
    j = minusk(i)
    print'(2(I4))',i,j
  enddo
  print*,'qktok2'
  do i=1,kpt_num
    do j=1,kpt_num
      k = qktok2(i,j)
      print'(3(I4))',i,j,k
    enddo
  enddo
  print*,'kconserv'
  do i=1,kpt_num
    do j=1,kpt_num
      do k=1,kpt_num
        l = kconserv(i,j,k)
        print'(4(I4))',i,j,k,l
      enddo
    enddo
  enddo

end
