module cfunctions
        use, intrinsic :: ISO_C_BINDING
      interface
         subroutine printcfglist(nint, ncfgs, cfglist) bind(C, name='printCFGList')
           import C_INT32_T, C_INT64_T
         integer(kind=C_INT32_T) :: nint
         integer(kind=C_INT32_T) :: ncfgs
         integer(kind=C_INT64_T) :: cfglist(nint,2,ncfgs)
       end subroutine printcfglist
      end interface
      interface
         subroutine getApqIJMatrixDims(Isomo, Jsomo, MS, rowsout, colsout) &
              bind(C, name='getApqIJMatrixDims')
           import C_INT32_T, C_INT64_T
           integer(kind=C_INT64_T),value,intent(in) :: Isomo ! CSFI
           integer(kind=C_INT64_T),value,intent(in) :: Jsomo ! CSFJ
           integer(kind=C_INT64_T),value,intent(in) :: MS    ! Ms = 2*Spin
           integer(kind=C_INT32_T),intent(out):: rowsout
           integer(kind=C_INT32_T),intent(out):: colsout
         end subroutine getApqIJMatrixDims
      end interface
      interface
         subroutine getApqIJMatrixDriver(Isomo, Jsomo, orbp, orbq,  &
              MS, NMO, CSFICSFJApqIJ, rowsmax, colsmax) bind(C, name='getApqIJMatrixDriverArrayInp')
           import C_INT32_T, C_INT64_T, C_DOUBLE
           integer(kind=C_INT64_T),value,intent(in) :: Isomo
           integer(kind=C_INT64_T),value,intent(in) :: Jsomo
           integer(kind=C_INT32_T),value,intent(in) :: orbp
           integer(kind=C_INT32_T),value,intent(in) :: orbq
           integer(kind=C_INT64_T),value,intent(in) :: MS
           integer(kind=C_INT64_T),value,intent(in) :: NMO
           integer(kind=C_INT32_T),intent(in) :: rowsmax
           integer(kind=C_INT32_T),intent(in) :: colsmax
           real   (kind=C_DOUBLE ),intent(out) :: CSFICSFJApqIJ(rowsmax,colsmax)
           !integer(kind=C_INT32_T),dimension(rowApqIJ,colApqIJ) :: ApqIJ
         end subroutine getApqIJMatrixDriver
      end interface
       interface
         subroutine getCSFtoDETTransformationMatrix(Isomo,&
              MS, rowsmax, colsmax, csftodetmatrix) bind(C, name='convertCSFtoDetBasis')
           import C_INT32_T, C_INT64_T, C_DOUBLE
           integer(kind=C_INT64_T),value,intent(in) :: Isomo
           integer(kind=C_INT64_T),value,intent(in) :: MS
           integer(kind=C_INT32_T),intent(in) :: rowsmax
           integer(kind=C_INT32_T),intent(in) :: colsmax
           real   (kind=C_DOUBLE ),intent(out) :: csftodetmatrix(rowsmax,colsmax)
         end subroutine getCSFtoDETTransformationMatrix
      end interface
       interface
         subroutine gramSchmidt(A, m, n, B) bind(C, name='gramSchmidt')
           import C_INT32_T, C_INT64_T, C_DOUBLE
           integer(kind=C_INT32_T),value,intent(in) :: m
           integer(kind=C_INT32_T),value,intent(in) :: n
           real (kind=C_DOUBLE ),intent(in) :: A(m,n)
           real (kind=C_DOUBLE ),intent(out) :: B(m,n)
         end subroutine gramSchmidt
      end interface
       interface
         subroutine gramSchmidt_qp(A, m, n, B) bind(C, name='gramSchmidt_qp')
           import C_INT32_T, C_INT64_T, C_DOUBLE
           integer(kind=C_INT32_T),value,intent(in) :: m
           integer(kind=C_INT32_T),value,intent(in) :: n
           real (kind=C_DOUBLE ),intent(in) :: A(m,n)
           real (kind=C_DOUBLE ),intent(out) :: B(m,n)
         end subroutine gramSchmidt_qp
      end interface
end module cfunctions

subroutine f_dgemm(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC) &
          bind(C, name='f_dgemm')
  use iso_c_binding
  implicit none
  character, intent(in), value :: TRANSA, TRANSB
  integer, intent(in), value   :: M,N,K,LDA,LDB,LDC
  double precision, intent(in), value :: ALPHA, BETA
  double precision, intent(in) :: A(LDA,*), B(LDB,*)
  double precision, intent(out) :: C(LDC,*)
  call dgemm(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
end subroutine


