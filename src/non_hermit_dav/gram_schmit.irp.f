subroutine bi_ortho_gram_schmidt(wi,vi,n,ni,wk,wk_schmidt)
 implicit none
 BEGIN_DOC
! you enter with a set of "ni" BI-ORTHONORMAL vectors of length "n" 
!
! vi(j,i) = <j|vi>, wi(j,i) = <j|wi>, <vi|wj> = delta_{ij} S_ii, S_ii =<vi|wi>
!
! and a vector vk(j) = <j|vk>
!
! you go out with a vector vk_schmidt(j) = <j|vk_schmidt> 
!
! which is Gram-Schmidt orthonormalized with respect to the "vi"
!
! <vi|wk_schmidt> = 0 
!
! |wk_schmidt> = |wk> - \sum_{i=1}^ni (<vi|wk>/<vi|wi>) |wi> 
!
! according to Eq. (5), (6) of Computers Structures, Vol 56, No. 4, pp 605-613, 1995
!
! https://doi.org/10.1016/0045-7949(94)00565-K
 END_DOC
 integer, intent(in) :: n,ni
 double precision, intent(in) :: wi(n,ni),vi(n,ni),wk(n)
 double precision, intent(out):: wk_schmidt(n)
 double precision :: vi_wk,u_dot_v,tmp,u_dot_u
 double precision, allocatable :: sii(:)
 integer :: i,j
 allocate( sii(ni) )
 wk_schmidt = wk
 do i = 1, ni
  sii(i) = u_dot_v(vi(1,i),wi(1,i),n)
 enddo
! do i = 1, n
!  print*,i,'wk',wk(i)
! enddo
! print*,''
! print*,''
 do i = 1, ni
!  print*,'i',i
  ! Gram-Schmidt 
  vi_wk = u_dot_v(vi(1,i),wk,n)
  vi_wk = vi_wk / sii(i)
!  print*,''
  do j = 1, n
!   print*,j,vi_wk,wi(j,i)
   wk_schmidt(j) -= vi_wk * wi(j,i)
  enddo
 enddo
 tmp = u_dot_u(wk_schmidt,n)
 tmp = 1.d0/dsqrt(tmp)
 wk_schmidt = tmp * wk_schmidt
! do j = 1, n
!  print*,j,'wk_scc',wk_schmidt(j)
! enddo
! pause
end
