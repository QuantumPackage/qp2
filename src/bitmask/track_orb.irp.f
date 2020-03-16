BEGIN_PROVIDER [ double precision, mo_coef_begin_iteration, (ao_num,mo_num) ]
   implicit none
   BEGIN_DOC
   ! Void provider to store the coefficients of the |MO| basis at the beginning of the SCF iteration
   !
   ! Useful to track some orbitals
   END_DOC
END_PROVIDER

BEGIN_PROVIDER [ complex*16, mo_coef_begin_iteration_complex, (ao_num,mo_num) ]
   implicit none
   BEGIN_DOC
   ! Void provider to store the coefficients of the |MO| basis at the beginning of the SCF iteration
   !
   ! Useful to track some orbitals
   END_DOC
END_PROVIDER

subroutine initialize_mo_coef_begin_iteration
 implicit none
 BEGIN_DOC
 !
 ! Initialize :c:data:`mo_coef_begin_iteration` to the current :c:data:`mo_coef`
 END_DOC
 if (is_complex) then
   mo_coef_begin_iteration_complex = mo_coef_complex
 else
   mo_coef_begin_iteration = mo_coef
 endif
end

subroutine reorder_core_orb
  implicit none
  BEGIN_DOC
  ! TODO: test for complex
 ! routines that takes the current :c:data:`mo_coef` and reorder the core orbitals (see :c:data:`list_core` and :c:data:`n_core_orb`) according to the overlap with :c:data:`mo_coef_begin_iteration`
  END_DOC
  integer :: i,j,iorb
  integer :: k,l
  integer, allocatable :: index_core_orb(:),iorder(:)
  double precision, allocatable :: accu(:)
  integer :: i1,i2
  if (is_complex) then
    complex*16, allocatable :: accu_c(:)
    allocate(accu(mo_num),accu_c(mo_num),index_core_orb(n_core_orb),iorder(mo_num))
    do i = 1, n_core_orb
      iorb = list_core(i)
      do j = 1, mo_num
        accu(j) = 0.d0
        accu_c(j) = (0.d0,0.d0)
        iorder(j) = j
        do k = 1, ao_num
          do l = 1, ao_num
            accu_c(j) += dconjg(mo_coef_begin_iteration_complex(k,iorb)) * &
                         mo_coef_complex(l,j) * ao_overlap_complex(k,l)
          enddo
        enddo
        accu(j) = -cdabs(accu_c(j))
      enddo
      call dsort(accu,iorder,mo_num)
      index_core_orb(i) = iorder(1)
    enddo

    complex*16 :: x_c
    do j = 1, n_core_orb
      i1 = list_core(j)
      i2 = index_core_orb(j)
      do i=1,ao_num
         x_c = mo_coef_complex(i,i1)
         mo_coef_complex(i,i1) = mo_coef_complex(i,i2)
         mo_coef_complex(i,i2) = x_c
      enddo
    enddo
    !call loc_cele_routine
   
    deallocate(accu,accu_c,index_core_orb, iorder)
  else
    allocate(accu(mo_num),index_core_orb(n_core_orb),iorder(mo_num))
   
    do i = 1, n_core_orb
      iorb = list_core(i)
      do j = 1, mo_num
        accu(j) = 0.d0
        iorder(j) = j
        do k = 1, ao_num
          do l = 1, ao_num
            accu(j) += mo_coef_begin_iteration(k,iorb) * mo_coef(l,j) * ao_overlap(k,l)
          enddo
        enddo
        accu(j) = -dabs(accu(j))
      enddo
      call dsort(accu,iorder,mo_num)
      index_core_orb(i) = iorder(1)
    enddo
   
    double precision :: x
    do j = 1, n_core_orb
      i1 = list_core(j)
      i2 = index_core_orb(j)
      do i=1,ao_num
         x = mo_coef(i,i1)
         mo_coef(i,i1) = mo_coef(i,i2)
         mo_coef(i,i2) = x
      enddo
    enddo
    !call loc_cele_routine
   
    deallocate(accu,index_core_orb, iorder)
  endif
end
