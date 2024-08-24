! First gradient

subroutine first_gradient_opt(n,v_grad)

  include 'constants.h'

  implicit none

  !===================================================================
  ! Compute the gradient of energy with respects to orbital rotations
  !===================================================================

  ! Check if read_wf = true, else :
  ! qp set determinant read_wf true

  END_DOC

  ! in
  integer, intent(in) :: n
  ! n : integer, n = mo_num*(mo_num-1)/2
  
  ! out
  double precision, intent(out) :: v_grad(n)
  ! v_grad : double precision vector of length n containeing the gradient

  ! internal
  double precision, allocatable :: grad(:,:),A(:,:)
  double precision :: norm
  integer :: i,p,q,r,s,t
  integer :: istate
  ! grad : double precision matrix containing the gradient before the permutation
  ! A : double precision matrix containing the gradient after the permutation
  ! norm : double precision number, the norm of the vector gradient
  ! i,p,q,r,s,t : integer, indexes 
  ! istate : integer, the electronic state

  ! Function
  double precision :: get_two_e_integral, norm2
  ! get_two_e_integral :  double precision function that gives the two e integrals
  ! norm2 : double precision function that gives the norm of a vector
 
  ! Provided :
  ! mo_one_e_integrals : mono e- integrals
  ! get_two_e_integral : two e- integrals
  ! one_e_dm_mo : one body density matrix (state average)
  ! two_e_dm_mo : two body density matrix (state average)

  !============
  ! Allocation
  !============

  allocate(grad(mo_num,mo_num),A(mo_num,mo_num))

  !=============
  ! Calculation
  !=============

  if (debug) then
    print*,'---first_gradient---'
  endif

  v_grad = 0d0

  do p = 1, mo_num
    do q = 1, mo_num
      grad(p,q) = 0d0
      do r = 1, mo_num
        grad(p,q) = grad(p,q) + mo_one_e_integrals(p,r) * one_e_dm_mo(r,q) &
                               - mo_one_e_integrals(r,q) * one_e_dm_mo(p,r)

      enddo

      do r = 1, mo_num
        do s = 1, mo_num
          do t= 1, mo_num

            grad(p,q) = grad(p,q) &
                + get_two_e_integral(p,t,r,s,mo_integrals_map) * two_e_dm_mo(r,s,q,t) &
                - get_two_e_integral(r,s,q,t,mo_integrals_map) * two_e_dm_mo(p,t,r,s)
          enddo
        enddo
      enddo
    enddo
  enddo

  ! Conversion mo_num*mo_num matrix to mo_num(mo_num-1)/2 vector
  do i=1,n
    call vec_to_mat_index(i,p,q)
    v_grad(i)=(grad(p,q) - grad(q,p))
  enddo  

  ! Display, vector containing the gradient elements 
  if (debug) then  
    print*,'Vector containing the gradient :'
    write(*,'(100(F10.5))') v_grad(1:n)
  endif  

  ! Norm of the vector
  norm = norm2(v_grad)
  print*, 'Norm : ', norm

  ! Matrix gradient
  A = 0d0
  do q=1,mo_num
    do p=1,mo_num
      A(p,q) = grad(p,q) - grad(q,p)
    enddo
  enddo

  ! Display, matrix containting the gradient elements
  if (debug) then    
    print*,'Matrix containing the gradient :'
    do i = 1, mo_num
      write(*,'(100(ES12.5))') A(i,1:mo_num)
    enddo
  endif

  !==============
  ! Deallocation
  !==============

  deallocate(grad,A)

  if (debug) then
    print*,'---End first_gradient---'
  endif

end subroutine
