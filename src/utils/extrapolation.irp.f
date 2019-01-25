subroutine extrapolate_data(N_data, data, pt2, output)
  implicit none
  BEGIN_DOC
  ! Extrapolate the data to the FCI limit
  END_DOC
  integer, intent(in)            :: N_data
  double precision, intent(in)   :: data(N_data)
  double precision, intent(in)   :: pt2(N_data)
  double precision, intent(out)  :: output(N_data)

  double precision :: data_rev(N_data), pt2_rev(N_data)
  integer :: ifit, i,j
  double precision :: ab(2), x(N_data,2), x_inv(2,N_data), y(N_data)

  do i=1,N_data
    data_rev(N_data+1-i) = data(i)
    pt2_rev(N_data+1-i) = pt2(i)
  enddo

  do i=1,N_data
    y(i) = data_rev(i)
    x(i,1) = 1.d0
    x(i,2) = pt2_rev(i)
  enddo
  do ifit=2,N_data
    call get_pseudo_inverse(x,size(x,1),ifit,2,x_inv,size(x_inv,1))
    ab = matmul(x_inv(1:2,1:ifit),y(1:ifit))
    output(ifit) = ab(1)
  enddo

end
