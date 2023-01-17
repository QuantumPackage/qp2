BEGIN_PROVIDER [ integer, ao_shell, (ao_num) ]
 implicit none
 BEGIN_DOC
 ! Index of the shell to which the AO corresponds
 END_DOC
 integer :: i, j, k, n
 k=0
 do i=1,shell_num
   n = shell_ang_mom(i)+1
   do j=1,(n*(n+1))/2
     k = k+1
     ao_shell(k) = i
   enddo
 enddo

END_PROVIDER


 BEGIN_PROVIDER [ double precision, ao_coef , (ao_num,ao_prim_num_max) ]
&BEGIN_PROVIDER [ double precision, ao_expo , (ao_num,ao_prim_num_max) ]
 implicit none
 BEGIN_DOC
! Primitive coefficients and exponents for each atomic orbital. Copied from shell info.
 END_DOC

 integer :: i, l
 do i=1,ao_num
   l = ao_shell(i)
   ao_coef(i,:) = shell_coef(l,:)
   ao_expo(i,:) = shell_expo(l,:)
 end do

END_PROVIDER


BEGIN_PROVIDER [ integer, ao_prim_num_max ]
 implicit none
 BEGIN_DOC
 ! Max number of primitives.
 END_DOC
 ao_prim_num_max = shell_prim_num_max
END_PROVIDER

BEGIN_PROVIDER [ integer, ao_first_of_shell, (shell_num) ]
 implicit none
 BEGIN_DOC
 ! Index of the shell to which the AO corresponds
 END_DOC
 integer :: i, j, k, n
 k=1
 do i=1,shell_num
   ao_first_of_shell(i) = k
   n = shell_ang_mom(i)+1
   k = k+(n*(n+1))/2
 enddo

END_PROVIDER

 BEGIN_PROVIDER [ double precision, ao_coef_normalized, (ao_num,ao_prim_num_max) ]
&BEGIN_PROVIDER [ double precision, ao_coef_normalization_factor, (ao_num) ]
  implicit none
  BEGIN_DOC
  ! Coefficients including the |AO| normalization
  END_DOC


  double precision               :: norm,overlap_x,overlap_y,overlap_z,C_A(3), c
  integer                        :: l, powA(3)
  integer, parameter             :: nz=100
  integer                        :: i,j,k

   ao_coef_normalized(:,:) = ao_coef(:,:)

  C_A = 0.d0

  do i=1,ao_num

    powA(1) = ao_power(i,1)
    powA(2) = ao_power(i,2)
    powA(3) = ao_power(i,3)

    ! GAMESS-type normalization of the primitives
    if (primitives_normalized) then
      do j=1,ao_prim_num(i)
        call overlap_gaussian_xyz(C_A,C_A,ao_expo(i,j),ao_expo(i,j), &
           powA,powA,overlap_x,overlap_y,overlap_z,norm,nz)
        ao_coef_normalized(i,j) = ao_coef_normalized(i,j)/dsqrt(norm)
      enddo
    endif
    ! Normalization of the contracted basis functions
    if (ao_normalized) then
      norm = 0.d0
      l = ao_shell(i)
      do j=1,ao_prim_num(i)
        do k=1,ao_prim_num(i)
          call overlap_gaussian_xyz(C_A,C_A,ao_expo(i,j),ao_expo(i,k),powA,powA,overlap_x,overlap_y,overlap_z,c,nz)
          norm = norm+c*ao_coef_normalized(i,j)*ao_coef_normalized(i,k)
        enddo
      enddo
      ao_coef_normalization_factor(i) = 1.d0/dsqrt(norm)
      ao_coef_normalized(i,:) *= ao_coef_normalization_factor(i)
    else
      ao_coef_normalization_factor(i) = 1.d0
    endif
  enddo

END_PROVIDER

 BEGIN_PROVIDER [ double precision, ao_coef_normalized_ordered, (ao_num,ao_prim_num_max) ]
&BEGIN_PROVIDER [ double precision, ao_expo_ordered, (ao_num,ao_prim_num_max) ]
  implicit none
  BEGIN_DOC
  ! Sorted primitives to accelerate 4 index |MO| transformation
  END_DOC

  integer                        :: iorder(ao_prim_num_max)
  double precision               :: d(ao_prim_num_max,2)
  integer                        :: i,j
  do i=1,ao_num
    do j=1,ao_prim_num(i)
      iorder(j) = j
      d(j,1) = ao_expo(i,j)
      d(j,2) = ao_coef_normalized(i,j)
    enddo
    call dsort(d(1,1),iorder,ao_prim_num(i))
    call dset_order(d(1,2),iorder,ao_prim_num(i))
    do j=1,ao_prim_num(i)
      ao_expo_ordered(i,j) = d(j,1)
      ao_coef_normalized_ordered(i,j) = d(j,2)
    enddo
  enddo
END_PROVIDER


BEGIN_PROVIDER [ double precision, ao_coef_normalized_ordered_transp, (ao_prim_num_max,ao_num) ]
  implicit none
  BEGIN_DOC
  ! Transposed :c:data:`ao_coef_normalized_ordered`
  END_DOC
  integer                        :: i,j
  do j=1, ao_num
    do i=1, ao_prim_num_max
      ao_coef_normalized_ordered_transp(i,j) = ao_coef_normalized_ordered(j,i)
    enddo
  enddo

END_PROVIDER

BEGIN_PROVIDER [ double precision, ao_expo_ordered_transp, (ao_prim_num_max,ao_num) ]
  implicit none
  BEGIN_DOC
  ! Transposed :c:data:`ao_expo_ordered`
  END_DOC
  integer                        :: i,j
  do j=1, ao_num
    do i=1, ao_prim_num_max
      ao_expo_ordered_transp(i,j) = ao_expo_ordered(j,i)
    enddo
  enddo

END_PROVIDER

 BEGIN_PROVIDER [ integer, ao_l, (ao_num) ]
&BEGIN_PROVIDER [ integer, ao_l_max  ]
&BEGIN_PROVIDER [ character*(128), ao_l_char, (ao_num) ]
 implicit none
 BEGIN_DOC
! :math:`l` value of the |AO|: :math`a+b+c` in :math:`x^a y^b z^c`
 END_DOC
 integer :: i
 do i=1,ao_num
   ao_l(i) = ao_power(i,1) + ao_power(i,2) + ao_power(i,3)
   ao_l_char(i) = l_to_character(ao_l(i))
 enddo
 ao_l_max = maxval(ao_l)
END_PROVIDER

integer function ao_power_index(nx,ny,nz)
  implicit none
  integer, intent(in)            :: nx, ny, nz
  BEGIN_DOC
  ! Unique index given to a triplet of powers:
  !
  ! :math:`\frac{1}{2} (l-n_x) (l-n_x+1) + n_z + 1`
  END_DOC
  integer                        :: l
  l = nx + ny + nz
  ao_power_index = ((l-nx)*(l-nx+1))/2 + nz + 1
end


BEGIN_PROVIDER [ character*(128), l_to_character, (0:7)]
 BEGIN_DOC
 ! Character corresponding to the "l" value of an |AO|
 END_DOC
 implicit none
 l_to_character(0)='s'
 l_to_character(1)='p'
 l_to_character(2)='d'
 l_to_character(3)='f'
 l_to_character(4)='g'
 l_to_character(5)='h'
 l_to_character(6)='i'
 l_to_character(7)='j'
END_PROVIDER


 BEGIN_PROVIDER [ integer, Nucl_N_Aos, (nucl_num)]
&BEGIN_PROVIDER [ integer, N_AOs_max ]
 implicit none
 BEGIN_DOC
 ! Number of |AOs| per atom
 END_DOC
 integer :: i
 Nucl_N_Aos = 0
 do i = 1, ao_num
  Nucl_N_Aos(ao_nucl(i)) +=1
 enddo
 N_AOs_max = maxval(Nucl_N_Aos)
END_PROVIDER

 BEGIN_PROVIDER [ integer, Nucl_Aos, (nucl_num,N_AOs_max)]
 implicit none
 BEGIN_DOC
 ! List of |AOs| centered on each atom
 END_DOC
 integer :: i
 integer, allocatable :: nucl_tmp(:)
 allocate(nucl_tmp(nucl_num))
 nucl_tmp = 0
 Nucl_Aos = 0
 do i = 1, ao_num
  nucl_tmp(ao_nucl(i))+=1
  Nucl_Aos(ao_nucl(i),nucl_tmp(ao_nucl(i))) = i
 enddo
 deallocate(nucl_tmp)
END_PROVIDER


 BEGIN_PROVIDER [ integer, Nucl_list_shell_Aos, (nucl_num,N_AOs_max)]
&BEGIN_PROVIDER [ integer, Nucl_num_shell_Aos, (nucl_num)]
 implicit none
 integer :: i,j,k
 BEGIN_DOC
 ! Index of the shell type |AOs| and of the corresponding |AOs|
 ! By convention, for p,d,f and g |AOs|, we take the index
 ! of the |AO| with the the corresponding power in the x axis
 END_DOC
 do i = 1, nucl_num
  Nucl_num_shell_Aos(i) = 0
  do j = 1, Nucl_N_Aos(i)
    if (ao_power(Nucl_Aos(i,j),1) == ao_l(Nucl_Aos(i,j))) then
     Nucl_num_shell_Aos(i)+=1
     Nucl_list_shell_Aos(i,Nucl_num_shell_Aos(i))=Nucl_Aos(i,j)
    endif
  enddo
 enddo

END_PROVIDER


BEGIN_PROVIDER [ character*(4), ao_l_char_space, (ao_num) ]
 implicit none
 BEGIN_DOC
! Converts an l value to a string
 END_DOC
 integer :: i
 character*(4) :: give_ao_character_space
 do i=1,ao_num

  if(ao_l(i)==0)then
  ! S type AO
   give_ao_character_space = 'S   '
  elseif(ao_l(i) == 1)then
  ! P type AO
   if(ao_power(i,1)==1)then
    give_ao_character_space = 'X   '
   elseif(ao_power(i,2) == 1)then
    give_ao_character_space = 'Y   '
   else
    give_ao_character_space = 'Z   '
   endif
  elseif(ao_l(i) == 2)then
  ! D type AO
   if(ao_power(i,1)==2)then
    give_ao_character_space = 'XX  '
   elseif(ao_power(i,2) == 2)then
    give_ao_character_space = 'YY  '
   elseif(ao_power(i,3) == 2)then
    give_ao_character_space = 'ZZ  '
   elseif(ao_power(i,1) == 1 .and. ao_power(i,2) == 1)then
    give_ao_character_space = 'XY  '
   elseif(ao_power(i,1) == 1 .and. ao_power(i,3) == 1)then
    give_ao_character_space = 'XZ  '
   else
    give_ao_character_space = 'YZ  '
   endif
  elseif(ao_l(i) == 3)then
  ! F type AO
   if(ao_power(i,1)==3)then
    give_ao_character_space = 'XXX '
   elseif(ao_power(i,2) == 3)then
    give_ao_character_space = 'YYY '
   elseif(ao_power(i,3) == 3)then
    give_ao_character_space = 'ZZZ '
   elseif(ao_power(i,1) == 2 .and. ao_power(i,2) == 1)then
    give_ao_character_space = 'XXY '
   elseif(ao_power(i,1) == 2 .and. ao_power(i,3) == 1)then
    give_ao_character_space = 'XXZ '
   elseif(ao_power(i,2) == 2 .and. ao_power(i,1) == 1)then
    give_ao_character_space = 'YYX '
   elseif(ao_power(i,2) == 2 .and. ao_power(i,3) == 1)then
    give_ao_character_space = 'YYZ '
   elseif(ao_power(i,3) == 2 .and. ao_power(i,1) == 1)then
    give_ao_character_space = 'ZZX '
   elseif(ao_power(i,3) == 2 .and. ao_power(i,2) == 1)then
    give_ao_character_space = 'ZZY '
   elseif(ao_power(i,3) == 1 .and. ao_power(i,2) == 1 .and. ao_power(i,3) == 1)then
    give_ao_character_space = 'XYZ '
   endif
  elseif(ao_l(i) == 4)then
  ! G type AO
   if(ao_power(i,1)==4)then
    give_ao_character_space = 'XXXX'
   elseif(ao_power(i,2) == 4)then
    give_ao_character_space = 'YYYY'
   elseif(ao_power(i,3) == 4)then
    give_ao_character_space = 'ZZZZ'
   elseif(ao_power(i,1) == 3 .and. ao_power(i,2) == 1)then
    give_ao_character_space = 'XXXY'
   elseif(ao_power(i,1) == 3 .and. ao_power(i,3) == 1)then
    give_ao_character_space = 'XXXZ'
   elseif(ao_power(i,2) == 3 .and. ao_power(i,1) == 1)then
    give_ao_character_space = 'YYYX'
   elseif(ao_power(i,2) == 3 .and. ao_power(i,3) == 1)then
    give_ao_character_space = 'YYYZ'
   elseif(ao_power(i,3) == 3 .and. ao_power(i,1) == 1)then
    give_ao_character_space = 'ZZZX'
   elseif(ao_power(i,3) == 3 .and. ao_power(i,2) == 1)then
    give_ao_character_space = 'ZZZY'
   elseif(ao_power(i,1) == 2 .and. ao_power(i,2) == 2)then
    give_ao_character_space = 'XXYY'
   elseif(ao_power(i,2) == 2 .and. ao_power(i,3) == 2)then
    give_ao_character_space = 'YYZZ'
   elseif(ao_power(i,1) == 2 .and. ao_power(i,2) == 1 .and. ao_power(i,3) == 1)then
    give_ao_character_space = 'XXYZ'
   elseif(ao_power(i,2) == 2 .and. ao_power(i,1) == 1 .and. ao_power(i,3) == 1)then
    give_ao_character_space = 'YYXZ'
   elseif(ao_power(i,3) == 2 .and. ao_power(i,1) == 1 .and. ao_power(i,2) == 1)then
    give_ao_character_space = 'ZZXY'
   endif
  endif
   ao_l_char_space(i) = give_ao_character_space
 enddo
END_PROVIDER
