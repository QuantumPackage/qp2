! Vector to matrix indexes
  
! *Compute the indexes p,q of a matrix element with the vector index i*

! Vector (i) -> lower diagonal matrix (p,q), p > q

! If a matrix is antisymmetric it can be reshaped as a vector. And the
! vector can be reshaped as an antisymmetric matrix

! \begin{align*}
! \begin{pmatrix}
! 0 & -1 & -2 & -4 \\
! 1 & 0  & -3 & -5 \\
! 2 & 3 & 0  & -6  \\
! 4 & 5 & 6 & 0
! \end{pmatrix}
! \Leftrightarrow
! \begin{pmatrix}
! 1 & 2 & 3 & 4 & 5 & 6
! \end{pmatrix}
! \end{align*}

! !!! Here the algorithm only work for the lower diagonal !!!

! Input:
! | i | integer | index in the vector |

! Ouput:
! | p,q | integer | corresponding indexes in the lower diagonal of a matrix |
! |     |         | p > q,                                                  |
! |     |         | p -> row,                                               |
! |     |         | q -> column                                             |


subroutine vec_to_mat_index(i,p,q)

  include 'pi.h'

  BEGIN_DOC
  ! Compute the indexes (p,q) of the element in the lower diagonal matrix knowing
  ! its index i a vector
  END_DOC

  implicit none

  ! Variables

  ! in
  integer,intent(in)   :: i
  
  ! out
  integer, intent(out) :: p,q
  
  ! internal 
  integer              :: a,b
  double precision     :: da

  da = 0.5d0*(1+ sqrt(1d0+8d0*DBLE(i)))
  a = INT(da) 
  if ((a*(a-1))/2==i) then
    p = a-1
  else
    p = a
  endif
  b = p*(p-1)/2
 
  ! Matrix element indexes
  p = p + 1
  q = i - b 

end subroutine
