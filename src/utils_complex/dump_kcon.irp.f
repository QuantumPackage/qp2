program dump_kcon
  call run
end

subroutine run
  use map_module
  implicit none

  integer ::i,j,k,l

  provide kconserv
  do i=1,kpt_num
    do j=1,kpt_num
      do k=1,kpt_num
        l = kconserv(i,j,k)
        print'(4(I4))',i,j,k,l
      enddo
    enddo
  enddo

end
