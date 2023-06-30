
! ---

BEGIN_PROVIDER [integer, List_all_comb_b2_size]

  implicit none

  PROVIDE j1b_type

  if((j1b_type .eq. 3) .or. (j1b_type .eq. 103)) then

    List_all_comb_b2_size = 2**nucl_num

  elseif((j1b_type .eq. 4) .or. (j1b_type .eq. 104)) then

    List_all_comb_b2_size = nucl_num + 1

  else 

    print *, 'j1b_type = ', j1b_type, 'is not implemented'
    stop

  endif

  print *, ' nb of linear terms in the envelope is ', List_all_comb_b2_size

END_PROVIDER

! ---

BEGIN_PROVIDER [integer, List_all_comb_b2, (nucl_num, List_all_comb_b2_size)]

  implicit none
  integer :: i, j

  if(nucl_num .gt. 32) then
    print *, ' nucl_num = ', nucl_num, '> 32'
    stop
  endif

  List_all_comb_b2 = 0

  do i = 0, List_all_comb_b2_size-1
    do j = 0, nucl_num-1
      if (btest(i,j)) then
        List_all_comb_b2(j+1,i+1) = 1
      endif
    enddo
  enddo

END_PROVIDER

! ---

 BEGIN_PROVIDER [ double precision, List_all_comb_b2_coef, (   List_all_comb_b2_size)]
&BEGIN_PROVIDER [ double precision, List_all_comb_b2_expo, (   List_all_comb_b2_size)]
&BEGIN_PROVIDER [ double precision, List_all_comb_b2_cent, (3, List_all_comb_b2_size)]

  implicit none
  integer          :: i, j, k, phase
  double precision :: tmp_alphaj, tmp_alphak
  double precision :: tmp_cent_x, tmp_cent_y, tmp_cent_z

  provide j1b_pen
  provide j1b_pen_coef

  List_all_comb_b2_coef = 0.d0
  List_all_comb_b2_expo = 0.d0
  List_all_comb_b2_cent = 0.d0

  if((j1b_type .eq. 3) .or. (j1b_type .eq. 103)) then

    do i = 1, List_all_comb_b2_size

      tmp_cent_x = 0.d0
      tmp_cent_y = 0.d0
      tmp_cent_z = 0.d0
      do j = 1, nucl_num
        tmp_alphaj = dble(List_all_comb_b2(j,i)) * j1b_pen(j)
        List_all_comb_b2_expo(i) += tmp_alphaj
        tmp_cent_x               += tmp_alphaj * nucl_coord(j,1)
        tmp_cent_y               += tmp_alphaj * nucl_coord(j,2)
        tmp_cent_z               += tmp_alphaj * nucl_coord(j,3)
      enddo

      if(List_all_comb_b2_expo(i) .lt. 1d-10) cycle

      List_all_comb_b2_cent(1,i) = tmp_cent_x / List_all_comb_b2_expo(i) 
      List_all_comb_b2_cent(2,i) = tmp_cent_y / List_all_comb_b2_expo(i)
      List_all_comb_b2_cent(3,i) = tmp_cent_z / List_all_comb_b2_expo(i)
    enddo

    ! ---

    do i = 1, List_all_comb_b2_size

      do j = 2, nucl_num, 1
        tmp_alphaj = dble(List_all_comb_b2(j,i)) * j1b_pen(j)
        do k = 1, j-1, 1
          tmp_alphak = dble(List_all_comb_b2(k,i)) * j1b_pen(k)

          List_all_comb_b2_coef(i) += tmp_alphaj * tmp_alphak * ( (nucl_coord(j,1) - nucl_coord(k,1)) * (nucl_coord(j,1) - nucl_coord(k,1)) &
                                                                + (nucl_coord(j,2) - nucl_coord(k,2)) * (nucl_coord(j,2) - nucl_coord(k,2)) &
                                                                + (nucl_coord(j,3) - nucl_coord(k,3)) * (nucl_coord(j,3) - nucl_coord(k,3)) )
        enddo
      enddo

      if(List_all_comb_b2_expo(i) .lt. 1d-10) cycle

      List_all_comb_b2_coef(i) = List_all_comb_b2_coef(i) / List_all_comb_b2_expo(i)
    enddo

    ! ---

    do i = 1, List_all_comb_b2_size

      phase = 0
      do j = 1, nucl_num
        phase += List_all_comb_b2(j,i)
      enddo

      List_all_comb_b2_coef(i) = (-1.d0)**dble(phase) * dexp(-List_all_comb_b2_coef(i))
    enddo

  elseif((j1b_type .eq. 4) .or. (j1b_type .eq. 104)) then

    List_all_comb_b2_coef(    1) = 1.d0
    List_all_comb_b2_expo(    1) = 0.d0
    List_all_comb_b2_cent(1:3,1) = 0.d0
    do i = 1, nucl_num
      List_all_comb_b2_coef(  i+1) = -1.d0 * j1b_pen_coef(i)
      List_all_comb_b2_expo(  i+1) = j1b_pen(i)
      List_all_comb_b2_cent(1,i+1) = nucl_coord(i,1)
      List_all_comb_b2_cent(2,i+1) = nucl_coord(i,2)
      List_all_comb_b2_cent(3,i+1) = nucl_coord(i,3)
    enddo

  else

    print *, 'j1b_type = ', j1b_type, 'is not implemented'
    stop

  endif

  !print *, ' coeff, expo & cent of list b2'
  !do i = 1, List_all_comb_b2_size
  !  print*, i, List_all_comb_b2_coef(i), List_all_comb_b2_expo(i)
  !  print*, List_all_comb_b2_cent(1,i), List_all_comb_b2_cent(2,i), List_all_comb_b2_cent(3,i)
  !enddo

END_PROVIDER

! ---

BEGIN_PROVIDER [ integer, List_all_comb_b3_size]

  implicit none
  double precision :: tmp

  if((j1b_type .eq. 3) .or. (j1b_type .eq. 103)) then

    List_all_comb_b3_size = 3**nucl_num

  elseif((j1b_type .eq. 4) .or. (j1b_type .eq. 104)) then

    tmp                   = 0.5d0 * dble(nucl_num) * (dble(nucl_num) + 3.d0)
    List_all_comb_b3_size = int(tmp) + 1

  else

    print *, 'j1b_type = ', j1b_type, 'is not implemented'
    stop

  endif

  print *, ' nb of linear terms in the square of the envelope is ', List_all_comb_b3_size

END_PROVIDER

! ---

BEGIN_PROVIDER [integer, List_all_comb_b3, (nucl_num, List_all_comb_b3_size)]

  implicit none
  integer              :: i, j, ii, jj
  integer, allocatable :: M(:,:), p(:)

  if(nucl_num .gt. 32) then
    print *, ' nucl_num = ', nucl_num, '> 32'
    stop
  endif

  List_all_comb_b3(:,:)                     = 0
  List_all_comb_b3(:,List_all_comb_b3_size) = 2

  allocate(p(nucl_num))
  p = 0

  do i = 2, List_all_comb_b3_size-1
    do j = 1, nucl_num

      ii = 0
      do jj = 1, j-1, 1
        ii = ii + p(jj) * 3**(jj-1)
      enddo
      p(j) = modulo(i-1-ii, 3**j) / 3**(j-1)

      List_all_comb_b3(j,i) = p(j)
    enddo
  enddo

END_PROVIDER

! ---

 BEGIN_PROVIDER [ double precision, List_all_comb_b3_coef, (   List_all_comb_b3_size)]
&BEGIN_PROVIDER [ double precision, List_all_comb_b3_expo, (   List_all_comb_b3_size)]
&BEGIN_PROVIDER [ double precision, List_all_comb_b3_cent, (3, List_all_comb_b3_size)]

  implicit none
  integer          :: i, j, k, phase
  integer          :: ii
  double precision :: tmp_alphaj, tmp_alphak, facto
  double precision :: tmp1, tmp2, tmp3, tmp4
  double precision :: xi, yi, zi, xj, yj, zj
  double precision :: dx, dy, dz, r2

  provide j1b_pen
  provide j1b_pen_coef

  List_all_comb_b3_coef = 0.d0
  List_all_comb_b3_expo = 0.d0
  List_all_comb_b3_cent = 0.d0

  if((j1b_type .eq. 3) .or. (j1b_type .eq. 103)) then

    do i = 1, List_all_comb_b3_size

      do j = 1, nucl_num
        tmp_alphaj = dble(List_all_comb_b3(j,i)) * j1b_pen(j)
        List_all_comb_b3_expo(i)   += tmp_alphaj
        List_all_comb_b3_cent(1,i) += tmp_alphaj * nucl_coord(j,1)
        List_all_comb_b3_cent(2,i) += tmp_alphaj * nucl_coord(j,2)
        List_all_comb_b3_cent(3,i) += tmp_alphaj * nucl_coord(j,3)

      enddo

      if(List_all_comb_b3_expo(i) .lt. 1d-10) cycle
      ASSERT(List_all_comb_b3_expo(i) .gt. 0d0)

      List_all_comb_b3_cent(1,i) = List_all_comb_b3_cent(1,i) / List_all_comb_b3_expo(i) 
      List_all_comb_b3_cent(2,i) = List_all_comb_b3_cent(2,i) / List_all_comb_b3_expo(i)
      List_all_comb_b3_cent(3,i) = List_all_comb_b3_cent(3,i) / List_all_comb_b3_expo(i)
    enddo

    ! ---

    do i = 1, List_all_comb_b3_size

      do j = 2, nucl_num, 1
        tmp_alphaj = dble(List_all_comb_b3(j,i)) * j1b_pen(j)
        do k = 1, j-1, 1
          tmp_alphak = dble(List_all_comb_b3(k,i)) * j1b_pen(k)

          List_all_comb_b3_coef(i) += tmp_alphaj * tmp_alphak * ( (nucl_coord(j,1) - nucl_coord(k,1)) * (nucl_coord(j,1) - nucl_coord(k,1)) &
                                                                + (nucl_coord(j,2) - nucl_coord(k,2)) * (nucl_coord(j,2) - nucl_coord(k,2)) &
                                                                + (nucl_coord(j,3) - nucl_coord(k,3)) * (nucl_coord(j,3) - nucl_coord(k,3)) )
        enddo
      enddo

      if(List_all_comb_b3_expo(i) .lt. 1d-10) cycle

      List_all_comb_b3_coef(i) = List_all_comb_b3_coef(i) / List_all_comb_b3_expo(i)
    enddo

    ! ---

    do i = 1, List_all_comb_b3_size

      facto = 1.d0
      phase = 0
      do j = 1, nucl_num
        tmp_alphaj = dble(List_all_comb_b3(j,i)) 

        facto *= 2.d0 / (gamma(tmp_alphaj+1.d0) * gamma(3.d0-tmp_alphaj))
        phase += List_all_comb_b3(j,i)
      enddo

      List_all_comb_b3_coef(i) = (-1.d0)**dble(phase) * facto * dexp(-List_all_comb_b3_coef(i))
    enddo

  elseif((j1b_type .eq. 4) .or. (j1b_type .eq. 104)) then

    ii = 1
    List_all_comb_b3_coef(    ii) = 1.d0
    List_all_comb_b3_expo(    ii) = 0.d0
    List_all_comb_b3_cent(1:3,ii) = 0.d0

    do i = 1, nucl_num
      ii = ii + 1
      List_all_comb_b3_coef(  ii) = -2.d0 * j1b_pen_coef(i)
      List_all_comb_b3_expo(  ii) = j1b_pen(i)
      List_all_comb_b3_cent(1,ii) = nucl_coord(i,1)
      List_all_comb_b3_cent(2,ii) = nucl_coord(i,2)
      List_all_comb_b3_cent(3,ii) = nucl_coord(i,3)
    enddo

    do i = 1, nucl_num
      ii = ii + 1
      List_all_comb_b3_coef(  ii) = 1.d0 * j1b_pen_coef(i) * j1b_pen_coef(i)
      List_all_comb_b3_expo(  ii) = 2.d0 * j1b_pen(i)
      List_all_comb_b3_cent(1,ii) = nucl_coord(i,1)
      List_all_comb_b3_cent(2,ii) = nucl_coord(i,2)
      List_all_comb_b3_cent(3,ii) = nucl_coord(i,3)
    enddo

    do i = 1, nucl_num-1

      tmp1 = j1b_pen(i)

      xi = nucl_coord(i,1)
      yi = nucl_coord(i,2)
      zi = nucl_coord(i,3)

      do j = i+1, nucl_num

        tmp2 = j1b_pen(j)
        tmp3 = tmp1 + tmp2
        tmp4 = 1.d0 / tmp3

        xj = nucl_coord(j,1)
        yj = nucl_coord(j,2)
        zj = nucl_coord(j,3)

        dx = xi - xj
        dy = yi - yj
        dz = zi - zj
        r2 = dx*dx + dy*dy + dz*dz
        
        ii = ii + 1
        ! x 2 to avoid doing integrals twice
        List_all_comb_b3_coef(  ii) = 2.d0 * dexp(-tmp1*tmp2*tmp4*r2) * j1b_pen_coef(i) * j1b_pen_coef(j)
        List_all_comb_b3_expo(  ii) = tmp3
        List_all_comb_b3_cent(1,ii) = tmp4 * (tmp1 * xi + tmp2 * xj)
        List_all_comb_b3_cent(2,ii) = tmp4 * (tmp1 * yi + tmp2 * yj)
        List_all_comb_b3_cent(3,ii) = tmp4 * (tmp1 * zi + tmp2 * zj)
      enddo
    enddo

  else

    print *, 'j1b_type = ', j1b_type, 'is not implemented'
    stop

  endif

  !print *, ' coeff, expo & cent of list b3'
  !do i = 1, List_all_comb_b3_size
  !  print*, i, List_all_comb_b3_coef(i), List_all_comb_b3_expo(i)
  !  print*, List_all_comb_b3_cent(1,i), List_all_comb_b3_cent(2,i), List_all_comb_b3_cent(3,i)
  !enddo

END_PROVIDER

! ---

