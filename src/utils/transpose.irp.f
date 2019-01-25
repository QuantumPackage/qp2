!DIR$ attributes forceinline :: transpose
recursive subroutine transpose(A,LDA,B,LDB,d1,d2)
 implicit none
 BEGIN_DOC
! Transpose input matrix A into output matrix B
 END_DOC
 integer, intent(in) :: d1, d2, LDA, LDB
 real, intent(in) :: A(LDA,d2)
 real, intent(out) :: B(LDB,d1)

 integer :: i,j,k
 if ( d2 < 32 ) then
   do j=1,d1
     !DIR$ LOOP COUNT (16)
     do i=1,d2
      B(i,j  ) = A(j  ,i)
     enddo
    enddo
   return
 else if (d1 > d2) then
 !DIR$ forceinline
   k=d1/2
 !DIR$ forceinline recursive
   call transpose(A(1,1),LDA,B(1,1),LDB,k,d2)
 !DIR$ forceinline recursive
   call transpose(A(k+1,1),LDA,B(1,k+1),LDB,d1-k,d2)
   return
 else
 !DIR$ forceinline
   k=d2/2
 !DIR$ forceinline recursive
   call transpose(A(1,k+1),LDA,B(k+1,1),LDB,d1,d2-k)
 !DIR$ forceinline recursive
   call transpose(A(1,1),LDA,B(1,1),LDB,d1,k)
   return
 endif

end

!DIR$ attributes forceinline :: dtranspose
recursive subroutine dtranspose(A,LDA,B,LDB,d1,d2)
 implicit none
 BEGIN_DOC
! Transpose input matrix A into output matrix B
 END_DOC
 integer, intent(in) :: d1, d2, LDA, LDB
 double precision, intent(in) :: A(LDA,d2)
 double precision, intent(out) :: B(LDB,d1)


! do j=1,d1
!   do i=1,d2
!    B(i,j  ) = A(j  ,i)
!   enddo
! enddo
! return

 integer :: i,j,k
 if ( d2 < 32 ) then
   do j=1,d1
     !DIR$ LOOP COUNT (16)
     do i=1,d2
      B(i,j  ) = A(j  ,i)
     enddo
    enddo
   return
 else if (d1 > d2) then
 !DIR$ forceinline
   k=d1/2
 !DIR$ forceinline recursive
   call dtranspose(A(1,1),LDA,B(1,1),LDB,k,d2)
 !DIR$ forceinline recursive
   call dtranspose(A(k+1,1),LDA,B(1,k+1),LDB,d1-k,d2)
   return
 else
 !DIR$ forceinline
   k=d2/2
 !DIR$ forceinline recursive
   call dtranspose(A(1,k+1),LDA,B(k+1,1),LDB,d1,d2-k)
 !DIR$ forceinline recursive
   call dtranspose(A(1,1),LDA,B(1,1),LDB,d1,k)
   return
 endif

end

