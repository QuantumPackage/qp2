
! ---

BEGIN_PROVIDER [integer, List_env1s_size]

  implicit none

  PROVIDE env_type

  if(env_type .eq. "None") then

    List_env1s_size = 1

  elseif(env_type .eq. "Prod_Gauss") then

    List_env1s_size = 2**nucl_num

  elseif(env_type .eq. "Sum_Gauss") then

    List_env1s_size = nucl_num + 1

  else 

    print *, ' Error in List_env1s_size: Unknown env_type = ', env_type
    stop

  endif

  print *, ' nb of 1s-Gaussian in the envelope = ', List_env1s_size

END_PROVIDER

! ---

BEGIN_PROVIDER [integer, List_env1s, (nucl_num, List_env1s_size)]

  implicit none
  integer :: i, j

  if(nucl_num .gt. 32) then
    print *, ' nucl_num = ', nucl_num, '> 32'
    stop
  endif

  List_env1s = 0

  do i = 0, List_env1s_size-1
    do j = 0, nucl_num-1
      if (btest(i,j)) then
        List_env1s(j+1,i+1) = 1
      endif
    enddo
  enddo

END_PROVIDER

! ---

 BEGIN_PROVIDER [ double precision, List_env1s_coef, (   List_env1s_size)]
&BEGIN_PROVIDER [ double precision, List_env1s_expo, (   List_env1s_size)]
&BEGIN_PROVIDER [ double precision, List_env1s_cent, (3, List_env1s_size)]

  implicit none
  integer          :: i, j, k, phase
  double precision :: tmp_alphaj, tmp_alphak
  double precision :: tmp_cent_x, tmp_cent_y, tmp_cent_z

  provide env_type env_expo env_coef

  if(env_type .eq. "None") then

    List_env1s_coef(    1) = 1.d0
    List_env1s_expo(    1) = 0.d0
    List_env1s_cent(1:3,1) = 0.d0

  elseif(env_type .eq. "Prod_Gauss") then

    List_env1s_coef = 0.d0
    List_env1s_expo = 0.d0
    List_env1s_cent = 0.d0

    do i = 1, List_env1s_size

      tmp_cent_x = 0.d0
      tmp_cent_y = 0.d0
      tmp_cent_z = 0.d0
      do j = 1, nucl_num
        tmp_alphaj = dble(List_env1s(j,i)) * env_expo(j)
        List_env1s_expo(i) += tmp_alphaj
        tmp_cent_x               += tmp_alphaj * nucl_coord(j,1)
        tmp_cent_y               += tmp_alphaj * nucl_coord(j,2)
        tmp_cent_z               += tmp_alphaj * nucl_coord(j,3)
      enddo

      if(List_env1s_expo(i) .lt. 1d-10) cycle

      List_env1s_cent(1,i) = tmp_cent_x / List_env1s_expo(i) 
      List_env1s_cent(2,i) = tmp_cent_y / List_env1s_expo(i)
      List_env1s_cent(3,i) = tmp_cent_z / List_env1s_expo(i)
    enddo

    ! ---

    do i = 1, List_env1s_size

      do j = 2, nucl_num, 1
        tmp_alphaj = dble(List_env1s(j,i)) * env_expo(j)
        do k = 1, j-1, 1
          tmp_alphak = dble(List_env1s(k,i)) * env_expo(k)

          List_env1s_coef(i) += tmp_alphaj * tmp_alphak * ( (nucl_coord(j,1) - nucl_coord(k,1)) * (nucl_coord(j,1) - nucl_coord(k,1)) &
                                                                + (nucl_coord(j,2) - nucl_coord(k,2)) * (nucl_coord(j,2) - nucl_coord(k,2)) &
                                                                + (nucl_coord(j,3) - nucl_coord(k,3)) * (nucl_coord(j,3) - nucl_coord(k,3)) )
        enddo
      enddo

      if(List_env1s_expo(i) .lt. 1d-10) cycle

      List_env1s_coef(i) = List_env1s_coef(i) / List_env1s_expo(i)
    enddo

    ! ---

    do i = 1, List_env1s_size

      phase = 0
      do j = 1, nucl_num
        phase += List_env1s(j,i)
      enddo

      List_env1s_coef(i) = (-1.d0)**dble(phase) * dexp(-List_env1s_coef(i))
    enddo

  elseif(env_type .eq. "Sum_Gauss") then

    List_env1s_coef(    1) = 1.d0
    List_env1s_expo(    1) = 0.d0
    List_env1s_cent(1:3,1) = 0.d0
    do i = 1, nucl_num
      List_env1s_coef(  i+1) = -1.d0 * env_coef(i)
      List_env1s_expo(  i+1) = env_expo(i)
      List_env1s_cent(1,i+1) = nucl_coord(i,1)
      List_env1s_cent(2,i+1) = nucl_coord(i,2)
      List_env1s_cent(3,i+1) = nucl_coord(i,3)
    enddo

  else

    print *, ' Error in List_env1s: Unknown env_type = ', env_type
    stop

  endif

END_PROVIDER

! ---

BEGIN_PROVIDER [integer, List_env1s_square_size]

  implicit none
  double precision :: tmp

  if(env_type .eq. "None") then

    List_env1s_square_size = 1

  elseif(env_type .eq. "Prod_Gauss") then

    List_env1s_square_size = 3**nucl_num

  elseif(env_type .eq. "Sum_Gauss") then

    tmp                   = 0.5d0 * dble(nucl_num) * (dble(nucl_num) + 3.d0)
    List_env1s_square_size = int(tmp) + 1

  else

    print *, ' Error in List_env1s_square_size: Unknown env_type = ', env_type
    stop

  endif

  print *, ' nb of 1s-Gaussian in the square of envelope = ', List_env1s_square_size

END_PROVIDER

! ---

BEGIN_PROVIDER [integer, List_env1s_square, (nucl_num, List_env1s_square_size)]

  implicit none
  integer              :: i, j, ii, jj
  integer, allocatable :: M(:,:), p(:)

  if(nucl_num .gt. 32) then
    print *, ' nucl_num = ', nucl_num, '> 32'
    stop
  endif

  List_env1s_square(:,:)                      = 0
  List_env1s_square(:,List_env1s_square_size) = 2

  allocate(p(nucl_num))
  p = 0

  do i = 2, List_env1s_square_size-1
    do j = 1, nucl_num

      ii = 0
      do jj = 1, j-1, 1
        ii = ii + p(jj) * 3**(jj-1)
      enddo
      p(j) = modulo(i-1-ii, 3**j) / 3**(j-1)

      List_env1s_square(j,i) = p(j)
    enddo
  enddo

END_PROVIDER

! ---

 BEGIN_PROVIDER [ double precision, List_env1s_square_coef, (   List_env1s_square_size)]
&BEGIN_PROVIDER [ double precision, List_env1s_square_expo, (   List_env1s_square_size)]
&BEGIN_PROVIDER [ double precision, List_env1s_square_cent, (3, List_env1s_square_size)]

  implicit none
  integer          :: i, j, k, phase
  integer          :: ii
  double precision :: tmp_alphaj, tmp_alphak, facto
  double precision :: tmp1, tmp2, tmp3, tmp4
  double precision :: xi, yi, zi, xj, yj, zj
  double precision :: dx, dy, dz, r2

  provide env_type env_expo env_coef

  if(env_type .eq. "None") then

    List_env1s_square_coef(    1) = 1.d0
    List_env1s_square_expo(    1) = 0.d0
    List_env1s_square_cent(1:3,1) = 0.d0

  elseif(env_type .eq. "Prod_Gauss") then

    List_env1s_square_coef = 0.d0
    List_env1s_square_expo = 0.d0
    List_env1s_square_cent = 0.d0

    do i = 1, List_env1s_square_size

      do j = 1, nucl_num
        tmp_alphaj = dble(List_env1s_square(j,i)) * env_expo(j)
        List_env1s_square_expo(i)   += tmp_alphaj
        List_env1s_square_cent(1,i) += tmp_alphaj * nucl_coord(j,1)
        List_env1s_square_cent(2,i) += tmp_alphaj * nucl_coord(j,2)
        List_env1s_square_cent(3,i) += tmp_alphaj * nucl_coord(j,3)

      enddo

      if(List_env1s_square_expo(i) .lt. 1d-10) cycle

      List_env1s_square_cent(1,i) = List_env1s_square_cent(1,i) / List_env1s_square_expo(i) 
      List_env1s_square_cent(2,i) = List_env1s_square_cent(2,i) / List_env1s_square_expo(i)
      List_env1s_square_cent(3,i) = List_env1s_square_cent(3,i) / List_env1s_square_expo(i)
    enddo

    ! ---

    do i = 1, List_env1s_square_size

      do j = 2, nucl_num, 1
        tmp_alphaj = dble(List_env1s_square(j,i)) * env_expo(j)
        do k = 1, j-1, 1
          tmp_alphak = dble(List_env1s_square(k,i)) * env_expo(k)

          List_env1s_square_coef(i) += tmp_alphaj * tmp_alphak * ( (nucl_coord(j,1) - nucl_coord(k,1)) * (nucl_coord(j,1) - nucl_coord(k,1)) &
                                                                + (nucl_coord(j,2) - nucl_coord(k,2)) * (nucl_coord(j,2) - nucl_coord(k,2)) &
                                                                + (nucl_coord(j,3) - nucl_coord(k,3)) * (nucl_coord(j,3) - nucl_coord(k,3)) )
        enddo
      enddo

      if(List_env1s_square_expo(i) .lt. 1d-10) cycle

      List_env1s_square_coef(i) = List_env1s_square_coef(i) / List_env1s_square_expo(i)
    enddo

    ! ---

    do i = 1, List_env1s_square_size

      facto = 1.d0
      phase = 0
      do j = 1, nucl_num
        tmp_alphaj = dble(List_env1s_square(j,i)) 

        facto *= 2.d0 / (gamma(tmp_alphaj+1.d0) * gamma(3.d0-tmp_alphaj))
        phase += List_env1s_square(j,i)
      enddo

      List_env1s_square_coef(i) = (-1.d0)**dble(phase) * facto * dexp(-List_env1s_square_coef(i))
    enddo

  elseif(env_type .eq. "Sum_Gauss") then

    ii = 1
    List_env1s_square_coef(    ii) = 1.d0
    List_env1s_square_expo(    ii) = 0.d0
    List_env1s_square_cent(1:3,ii) = 0.d0

    do i = 1, nucl_num
      ii = ii + 1
      List_env1s_square_coef(  ii) = -2.d0 * env_coef(i)
      List_env1s_square_expo(  ii) = env_expo(i)
      List_env1s_square_cent(1,ii) = nucl_coord(i,1)
      List_env1s_square_cent(2,ii) = nucl_coord(i,2)
      List_env1s_square_cent(3,ii) = nucl_coord(i,3)
    enddo

    do i = 1, nucl_num
      ii = ii + 1
      List_env1s_square_coef(  ii) = 1.d0 * env_coef(i) * env_coef(i)
      List_env1s_square_expo(  ii) = 2.d0 * env_expo(i)
      List_env1s_square_cent(1,ii) = nucl_coord(i,1)
      List_env1s_square_cent(2,ii) = nucl_coord(i,2)
      List_env1s_square_cent(3,ii) = nucl_coord(i,3)
    enddo

    do i = 1, nucl_num-1

      tmp1 = env_expo(i)

      xi = nucl_coord(i,1)
      yi = nucl_coord(i,2)
      zi = nucl_coord(i,3)

      do j = i+1, nucl_num

        tmp2 = env_expo(j)
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
        List_env1s_square_coef(  ii) = 2.d0 * dexp(-tmp1*tmp2*tmp4*r2) * env_coef(i) * env_coef(j)
        List_env1s_square_expo(  ii) = tmp3
        List_env1s_square_cent(1,ii) = tmp4 * (tmp1 * xi + tmp2 * xj)
        List_env1s_square_cent(2,ii) = tmp4 * (tmp1 * yi + tmp2 * yj)
        List_env1s_square_cent(3,ii) = tmp4 * (tmp1 * zi + tmp2 * zj)
      enddo
    enddo

  else

    print *, ' Error in List_env1s_square: Unknown env_type = ', env_type
    stop

  endif

END_PROVIDER

! ---

