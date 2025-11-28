BEGIN_PROVIDER [ integer, ao_extra_prim_num_max ]
 implicit none
 BEGIN_DOC
 ! Max number of primitives.
 END_DOC
 ao_extra_prim_num_max = maxval(ao_extra_prim_num)
END_PROVIDER

BEGIN_PROVIDER [ integer, ao_extra_shell, (ao_extra_num) ]
 implicit none
 BEGIN_DOC
 ! Index of the shell to which the ao_extra corresponds
 END_DOC
 integer :: i, j, k, n
 k=0
 do i=1,shell_num
   n = shell_ang_mom(i)+1
   do j=1,(n*(n+1))/2
     k = k+1
     ao_extra_shell(k) = i
   enddo
 enddo

END_PROVIDER

BEGIN_PROVIDER [ integer, ao_extra_first_of_shell, (shell_num) ]
 implicit none
 BEGIN_DOC
 ! Index of the shell to which the ao_extra corresponds
 END_DOC
 integer :: i, j, k, n
 k=1
 do i=1,shell_num
   ao_extra_first_of_shell(i) = k
   n = shell_ang_mom(i)+1
   k = k+(n*(n+1))/2
 enddo

END_PROVIDER

 BEGIN_PROVIDER [ double precision, ao_extra_coef_normalized, (ao_extra_num,ao_extra_prim_num_max) ]
&BEGIN_PROVIDER [ double precision, ao_extra_coef_normalization_factor, (ao_extra_num) ]
  implicit none
  BEGIN_DOC
  ! Coefficients including the |ao_extra| normalization
  END_DOC
  double precision               :: norm,overlap_x,overlap_y,overlap_z,C_A(3), c
  integer                        :: l, powA(3), nz
  integer                        :: i,j,k
  nz=100
  C_A(1) = 0.d0
  C_A(2) = 0.d0
  C_A(3) = 0.d0
  ao_extra_coef_normalized = 0.d0

  do i=1,ao_extra_num

!    powA(1) = ao_extra_power(i,1) +  ao_extra_power(i,2) +  ao_extra_power(i,3)
!    powA(2) = 0
!    powA(3) = 0
    powA(1) = ao_extra_power(i,1)
    powA(2) = ao_extra_power(i,2)
    powA(3) = ao_extra_power(i,3)

    ! Normalization of the primitives
    if (primitives_normalized_extra) then
      do j=1,ao_extra_prim_num(i)
        call overlap_gaussian_xyz(C_A,C_A,ao_extra_expo(i,j),ao_extra_expo(i,j), &
           powA,powA,overlap_x,overlap_y,overlap_z,norm,nz)
        ao_extra_coef_normalized(i,j) = ao_extra_coef(i,j)/dsqrt(norm)
      enddo
    else
      do j=1,ao_extra_prim_num(i)
        ao_extra_coef_normalized(i,j) = ao_extra_coef(i,j)
      enddo
    endif

    ! Normalization of the contracted basis functions
    if (ao_extra_normalized) then
      norm = 0.d0
      do j=1,ao_extra_prim_num(i)
        do k=1,ao_extra_prim_num(i)
          call overlap_gaussian_xyz(C_A,C_A,ao_extra_expo(i,j),ao_extra_expo(i,k),powA,powA,overlap_x,overlap_y,overlap_z,c,nz)
          norm = norm+c*ao_extra_coef_normalized(i,j)*ao_extra_coef_normalized(i,k)
        enddo
      enddo
      ao_extra_coef_normalization_factor(i) = 1.d0/dsqrt(norm)
    else
      ao_extra_coef_normalization_factor(i) = 1.d0
    endif
  enddo

END_PROVIDER

 BEGIN_PROVIDER [ double precision, ao_extra_coef_normalized_ordered, (ao_extra_num,ao_extra_prim_num_max) ]
&BEGIN_PROVIDER [ double precision, ao_extra_expo_ordered, (ao_extra_num,ao_extra_prim_num_max) ]
  implicit none
  BEGIN_DOC
  ! Sorted primitives to accelerate 4 index |MO| transformation
  END_DOC

  integer                        :: iorder(ao_extra_prim_num_max)
  double precision               :: d(ao_extra_prim_num_max,2)
  integer                        :: i,j
  do i=1,ao_extra_num
    do j=1,ao_extra_prim_num(i)
      iorder(j) = j
      d(j,1) = ao_extra_expo(i,j)
      d(j,2) = ao_extra_coef_normalized(i,j)
    enddo
    call dsort(d(1,1),iorder,ao_extra_prim_num(i))
    call dset_order(d(1,2),iorder,ao_extra_prim_num(i))
    do j=1,ao_extra_prim_num(i)
      ao_extra_expo_ordered(i,j) = d(j,1)
      ao_extra_coef_normalized_ordered(i,j) = d(j,2)
    enddo
  enddo
END_PROVIDER


BEGIN_PROVIDER [ double precision, ao_extra_coef_normalized_ordered_transp, (ao_extra_prim_num_max,ao_extra_num) ]
  implicit none
  BEGIN_DOC
  ! Transposed :c:data:`ao_extra_coef_normalized_ordered`
  END_DOC
  integer                        :: i,j
  do j=1, ao_extra_num
    do i=1, ao_extra_prim_num_max
      ao_extra_coef_normalized_ordered_transp(i,j) = ao_extra_coef_normalized_ordered(j,i)
    enddo
  enddo

END_PROVIDER

BEGIN_PROVIDER [ double precision, ao_extra_expo_ordered_transp, (ao_extra_prim_num_max,ao_extra_num) ]
  implicit none
  BEGIN_DOC
  ! Transposed :c:data:`ao_extra_expo_ordered`
  END_DOC
  integer                        :: i,j
  do j=1, ao_extra_num
    do i=1, ao_extra_prim_num_max
      ao_extra_expo_ordered_transp(i,j) = ao_extra_expo_ordered(j,i)
    enddo
  enddo

END_PROVIDER

 BEGIN_PROVIDER [ integer, ao_extra_l, (ao_extra_num) ]
&BEGIN_PROVIDER [ integer, ao_extra_l_max  ]
&BEGIN_PROVIDER [ character*(128), ao_extra_l_char, (ao_extra_num) ]
 implicit none
 BEGIN_DOC
! :math:`l` value of the |ao_extra|: :math`a+b+c` in :math:`x^a y^b z^c`
 END_DOC
 integer :: i
 do i=1,ao_extra_num
   ao_extra_l(i) = ao_extra_power(i,1) + ao_extra_power(i,2) + ao_extra_power(i,3)
   ao_extra_l_char(i) = l_to_character(ao_extra_l(i))
 enddo
 ao_extra_l_max = maxval(ao_extra_l)
END_PROVIDER

integer function ao_extra_power_index(nx,ny,nz)
  implicit none
  integer, intent(in)            :: nx, ny, nz
  BEGIN_DOC
  ! Unique index given to a triplet of powers:
  !
  ! :math:`\frac{1}{2} (l-n_x) (l-n_x+1) + n_z + 1`
  END_DOC
  integer                        :: l
  l = nx + ny + nz
  ao_extra_power_index = ((l-nx)*(l-nx+1))/2 + nz + 1
end


 BEGIN_PROVIDER [ integer, Nucl_N_ao_extras, (extra_nucl_num)]
&BEGIN_PROVIDER [ integer, N_ao_extras_max ]
 implicit none
 BEGIN_DOC
 ! Number of |ao_extras| per atom
 END_DOC
 integer :: i
 Nucl_N_ao_extras = 0
 do i = 1, ao_extra_num
  Nucl_N_ao_extras(ao_extra_nucl(i)) +=1
 enddo
 N_ao_extras_max = maxval(Nucl_N_ao_extras)
END_PROVIDER

 BEGIN_PROVIDER [ integer, Nucl_ao_extras, (extra_nucl_num,N_ao_extras_max)]
 implicit none
 BEGIN_DOC
 ! List of |ao_extras| centered on each atom
 END_DOC
 integer :: i
 integer, allocatable :: nucl_tmp(:)
 allocate(nucl_tmp(nucl_num))
 nucl_tmp = 0
 Nucl_ao_extras = 0
 do i = 1, ao_extra_num
  nucl_tmp(ao_extra_nucl(i))+=1
  Nucl_ao_extras(ao_extra_nucl(i),nucl_tmp(ao_extra_nucl(i))) = i
 enddo
 deallocate(nucl_tmp)
END_PROVIDER


 BEGIN_PROVIDER [ integer, Nucl_list_shell_ao_extras, (extra_nucl_num,N_ao_extras_max)]
&BEGIN_PROVIDER [ integer, Nucl_num_shell_ao_extras, (nucl_num)]
 implicit none
 integer :: i,j,k
 BEGIN_DOC
 ! Index of the shell type |ao_extras| and of the corresponding |ao_extras|
 ! By convention, for p,d,f and g |ao_extras|, we take the index
 ! of the |ao_extra| with the the corresponding power in the x axis
 END_DOC
 do i = 1, extra_nucl_num
  Nucl_num_shell_ao_extras(i) = 0
  do j = 1, Nucl_N_ao_extras(i)
    if (ao_extra_power(Nucl_ao_extras(i,j),1) == ao_extra_l(Nucl_ao_extras(i,j))) then
     Nucl_num_shell_ao_extras(i)+=1
     Nucl_list_shell_ao_extras(i,Nucl_num_shell_ao_extras(i))=Nucl_ao_extras(i,j)
    endif
  enddo
 enddo

END_PROVIDER


BEGIN_PROVIDER [ character*(4), ao_extra_l_char_space, (ao_extra_num) ]
 implicit none
 BEGIN_DOC
! Converts an l value to a string
 END_DOC
 integer :: i
 character*(4) :: give_ao_extra_character_space
 do i=1,ao_extra_num

  if(ao_extra_l(i)==0)then
  ! S type ao_extra
   give_ao_extra_character_space = 'S   '
  elseif(ao_extra_l(i) == 1)then
  ! P type ao_extra
   if(ao_extra_power(i,1)==1)then
    give_ao_extra_character_space = 'X   '
   elseif(ao_extra_power(i,2) == 1)then
    give_ao_extra_character_space = 'Y   '
   else
    give_ao_extra_character_space = 'Z   '
   endif
  elseif(ao_extra_l(i) == 2)then
  ! D type ao_extra
   if(ao_extra_power(i,1)==2)then
    give_ao_extra_character_space = 'XX  '
   elseif(ao_extra_power(i,2) == 2)then
    give_ao_extra_character_space = 'YY  '
   elseif(ao_extra_power(i,3) == 2)then
    give_ao_extra_character_space = 'ZZ  '
   elseif(ao_extra_power(i,1) == 1 .and. ao_extra_power(i,2) == 1)then
    give_ao_extra_character_space = 'XY  '
   elseif(ao_extra_power(i,1) == 1 .and. ao_extra_power(i,3) == 1)then
    give_ao_extra_character_space = 'XZ  '
   else
    give_ao_extra_character_space = 'YZ  '
   endif
  elseif(ao_extra_l(i) == 3)then
  ! F type ao_extra
   if(ao_extra_power(i,1)==3)then
    give_ao_extra_character_space = 'XXX '
   elseif(ao_extra_power(i,2) == 3)then
    give_ao_extra_character_space = 'YYY '
   elseif(ao_extra_power(i,3) == 3)then
    give_ao_extra_character_space = 'ZZZ '
   elseif(ao_extra_power(i,1) == 2 .and. ao_extra_power(i,2) == 1)then
    give_ao_extra_character_space = 'XXY '
   elseif(ao_extra_power(i,1) == 2 .and. ao_extra_power(i,3) == 1)then
    give_ao_extra_character_space = 'XXZ '
   elseif(ao_extra_power(i,2) == 2 .and. ao_extra_power(i,1) == 1)then
    give_ao_extra_character_space = 'YYX '
   elseif(ao_extra_power(i,2) == 2 .and. ao_extra_power(i,3) == 1)then
    give_ao_extra_character_space = 'YYZ '
   elseif(ao_extra_power(i,3) == 2 .and. ao_extra_power(i,1) == 1)then
    give_ao_extra_character_space = 'ZZX '
   elseif(ao_extra_power(i,3) == 2 .and. ao_extra_power(i,2) == 1)then
    give_ao_extra_character_space = 'ZZY '
   elseif(ao_extra_power(i,3) == 1 .and. ao_extra_power(i,2) == 1 .and. ao_extra_power(i,3) == 1)then
    give_ao_extra_character_space = 'XYZ '
   endif
  elseif(ao_extra_l(i) == 4)then
  ! G type ao_extra
   if(ao_extra_power(i,1)==4)then
    give_ao_extra_character_space = 'XXXX'
   elseif(ao_extra_power(i,2) == 4)then
    give_ao_extra_character_space = 'YYYY'
   elseif(ao_extra_power(i,3) == 4)then
    give_ao_extra_character_space = 'ZZZZ'
   elseif(ao_extra_power(i,1) == 3 .and. ao_extra_power(i,2) == 1)then
    give_ao_extra_character_space = 'XXXY'
   elseif(ao_extra_power(i,1) == 3 .and. ao_extra_power(i,3) == 1)then
    give_ao_extra_character_space = 'XXXZ'
   elseif(ao_extra_power(i,2) == 3 .and. ao_extra_power(i,1) == 1)then
    give_ao_extra_character_space = 'YYYX'
   elseif(ao_extra_power(i,2) == 3 .and. ao_extra_power(i,3) == 1)then
    give_ao_extra_character_space = 'YYYZ'
   elseif(ao_extra_power(i,3) == 3 .and. ao_extra_power(i,1) == 1)then
    give_ao_extra_character_space = 'ZZZX'
   elseif(ao_extra_power(i,3) == 3 .and. ao_extra_power(i,2) == 1)then
    give_ao_extra_character_space = 'ZZZY'
   elseif(ao_extra_power(i,1) == 2 .and. ao_extra_power(i,2) == 2)then
    give_ao_extra_character_space = 'XXYY'
   elseif(ao_extra_power(i,2) == 2 .and. ao_extra_power(i,3) == 2)then
    give_ao_extra_character_space = 'YYZZ'
   elseif(ao_extra_power(i,1) == 2 .and. ao_extra_power(i,2) == 1 .and. ao_extra_power(i,3) == 1)then
    give_ao_extra_character_space = 'XXYZ'
   elseif(ao_extra_power(i,2) == 2 .and. ao_extra_power(i,1) == 1 .and. ao_extra_power(i,3) == 1)then
    give_ao_extra_character_space = 'YYXZ'
   elseif(ao_extra_power(i,3) == 2 .and. ao_extra_power(i,1) == 1 .and. ao_extra_power(i,2) == 1)then
    give_ao_extra_character_space = 'ZZXY'
   endif
  endif
   ao_extra_l_char_space(i) = give_ao_extra_character_space
 enddo
END_PROVIDER


 ---

double precision function ao_extra_value(i, r)

  BEGIN_DOC
  ! Returns the value of the i-th ao belonging the EXTRA BASIS at point $\textbf{r}$
  END_DOC

  implicit none
  integer,          intent(in) :: i
  double precision, intent(in) :: r(3)

  integer          :: m, num_ao
  integer          :: power_ao(3)
  double precision :: center_ao(3)
  double precision :: beta
  double precision :: accu, dx, dy, dz, r2

  num_ao = ao_extra_nucl(i)
  power_ao(1:3) = ao_extra_power(i,1:3)
  center_ao(1:3) = extra_nucl_coord(num_ao,1:3)
  dx = r(1) - center_ao(1)
  dy = r(2) - center_ao(2)
  dz = r(3) - center_ao(3)
  r2 = dx*dx + dy*dy + dz*dz
  dx = dx**power_ao(1)
  dy = dy**power_ao(2)
  dz = dz**power_ao(3)
 
  accu = 0.d0
  do m = 1, ao_extra_prim_num(i)
    beta = ao_extra_expo_ordered_transp(m,i)
    accu += ao_extra_coef_normalized_ordered_transp(m,i) * dexp(-beta*r2)
  enddo
  ao_extra_value = accu * dx * dy * dz

end


subroutine give_all_aos_extra_at_r(r, tmp_array)
 implicit none
  BEGIN_dOC
  !
  ! input  : r == r(1) = x and so on
  !
  ! output : tmp_array(i) = EXTRA aos(i) evaluated in $\textbf{r}$
  !
  END_DOC
  double precision, intent(in)  :: r(3)
  double precision, intent(out) :: tmp_array(ao_extra_num)
  integer :: i
 double precision :: ao_extra_value
 do i = 1, ao_extra_num
  tmp_array(i) = ao_extra_value(i, r)
 enddo
end

double precision function extra_density_at_r(r)
 implicit none
  BEGIN_dOC
  !
  ! input  : r == r(1) = x and so on
  !
  ! output : density corresponding to the extra system at r
  !
  END_DOC
 double precision, intent(in)  :: r(3)
 integer :: mu,nu
 double precision, allocatable :: tmp_array(:)
 allocate(tmp_array(ao_extra_num))
 call give_all_aos_extra_at_r(r, tmp_array) 
 extra_density_at_r = 0.d0
 do nu = 1, ao_extra_num
  do mu = 1, ao_extra_num
   extra_density_at_r += ao_extra_one_e_dm(mu,nu,1) * tmp_array(mu) * tmp_array(nu)
  enddo
 enddo
end

BEGIN_PROVIDER [ double precision, ao_extra_one_e_dm_at_extra_r, (n_points_extra_final_grid)]
 implicit none
 BEGIN_DOC
! ao_extra_one_e_dm_at_extra_r(i) = extra density on the extra grid points
 END_DOC
 integer :: i
 double precision :: extra_density_at_r, r(3)
 do i = 1, n_points_extra_final_grid
  r(1) = final_grid_points_extra(1,i)
  r(2) = final_grid_points_extra(2,i)
  r(3) = final_grid_points_extra(3,i)
  ao_extra_one_e_dm_at_extra_r(i) = extra_density_at_r(r)
 enddo
END_PROVIDER 
