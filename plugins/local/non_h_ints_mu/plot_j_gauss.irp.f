program plot_j_gauss
 implicit none
 double precision :: xmin, xmax, x, dx
 double precision :: mu_min, mu_max, mu, d_mu
 double precision :: pot_j_gauss,j_mu_simple,j_gauss_simple,pot_j_mu
 double precision, allocatable :: mu_tab(:),j_mu(:),j_mu_gauss(:)
 double precision, allocatable :: w_mu(:), w_mu_gauss(:)

 character*(128) :: output
 integer :: getUnitAndOpen
 integer :: i_unit_output_wee_gauss,i_unit_output_wee_mu
 integer :: i_unit_output_j_gauss,i_unit_output_j_mu
 output=trim(ezfio_filename)//'.w_ee_mu_gauss'
 i_unit_output_wee_gauss = getUnitAndOpen(output,'w')
 output=trim(ezfio_filename)//'.w_ee_mu'
 i_unit_output_wee_mu = getUnitAndOpen(output,'w')
 output=trim(ezfio_filename)//'.j_mu_gauss'
 i_unit_output_j_gauss = getUnitAndOpen(output,'w')
 output=trim(ezfio_filename)//'.j_mu'
 i_unit_output_j_mu = getUnitAndOpen(output,'w')

 integer :: npt, i, j, n_mu
 n_mu = 3
 allocate(mu_tab(n_mu),j_mu(n_mu),j_mu_gauss(n_mu),w_mu(n_mu), w_mu_gauss(n_mu))
 mu_min = 0.5d0
 mu_max = 2.d0
 d_mu = (mu_max - mu_min)/dble(n_mu)
 mu = mu_min
 do i = 1, n_mu
  mu_tab(i) = mu
  print*,'mu = ',mu
  mu += d_mu
 enddo
 mu_tab(1) = 0.9d0
 mu_tab(2) = 0.95d0
 mu_tab(3) = 1.d0

 xmin = 0.01d0
 xmax = 10.d0
 npt = 1000
 dx = (xmax - xmin)/dble(npt)
 x = xmin
 do i = 1, npt
  do j = 1, n_mu
   mu = mu_tab(j)
   w_mu_gauss(j) = pot_j_gauss(x,mu)
   w_mu(j) = pot_j_mu(x,mu)
   j_mu(j) = j_mu_simple(x,mu)
   j_mu_gauss(j) = j_gauss_simple(x,mu) + j_mu(j)
  enddo
  write(i_unit_output_wee_gauss,'(100(F16.10,X))')x,w_mu_gauss(:)
  write(i_unit_output_wee_mu,'(100(F16.10,X))')x,w_mu(:)
  write(i_unit_output_j_gauss,'(100(F16.10,X))')x,j_mu_gauss(:)
  write(i_unit_output_j_mu,'(100(F16.10,X))')x,j_mu(:)
  x += dx
 enddo


end
