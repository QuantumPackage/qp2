double precision function step_function_becke(x)
  implicit none
  BEGIN_DOC
! Step function of the Becke paper (1988, JCP,88(4))
  END_DOC
  double precision, intent(in)   :: x
  double precision               :: f_function_becke
  integer                        :: i,n_max_becke

  step_function_becke = f_function_becke(x)
  do i = 1, 4
    step_function_becke = f_function_becke(step_function_becke)
  enddo
  step_function_becke = 0.5d0*(1.d0 - step_function_becke)
end

double precision function f_function_becke(x)
  implicit none
  double precision, intent(in)   :: x
  f_function_becke = 1.5d0 * x - 0.5d0 * x*x*x
end

! ---

double precision function cell_function_becke(r, atom_number)

  BEGIN_DOC
  ! atom_number :: atom on which the cell function of Becke (1988, JCP,88(4))
  ! r(1:3)                       :: x,y,z coordinantes of the current point
  END_DOC

  implicit none
  double precision, intent(in) :: r(3)
  integer,          intent(in) :: atom_number
  integer                      :: j
  double precision             :: mu_ij, nu_ij
  double precision             :: distance_i, distance_j, step_function_becke

  distance_i  = (r(1) - nucl_coord_transp(1,atom_number) ) * (r(1) - nucl_coord_transp(1,atom_number))
  distance_i += (r(2) - nucl_coord_transp(2,atom_number) ) * (r(2) - nucl_coord_transp(2,atom_number))
  distance_i += (r(3) - nucl_coord_transp(3,atom_number) ) * (r(3) - nucl_coord_transp(3,atom_number))
  distance_i  = dsqrt(distance_i)

  cell_function_becke = 1.d0
  do j = 1, nucl_num
    if(j==atom_number) cycle

    distance_j  = (r(1) - nucl_coord_transp(1,j) ) * (r(1) - nucl_coord_transp(1,j))
    distance_j += (r(2) - nucl_coord_transp(2,j) ) * (r(2) - nucl_coord_transp(2,j))
    distance_j += (r(3) - nucl_coord_transp(3,j) ) * (r(3) - nucl_coord_transp(3,j))
    distance_j  = dsqrt(distance_j)

    mu_ij = (distance_i - distance_j) * nucl_dist_inv(atom_number,j)
    nu_ij = mu_ij + slater_bragg_type_inter_distance_ua(atom_number,j) * (1.d0 - mu_ij*mu_ij)

    cell_function_becke *= step_function_becke(nu_ij)
  enddo

  return
end

