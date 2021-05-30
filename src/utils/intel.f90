module intel
  use, intrinsic :: iso_c_binding
  interface
    subroutine ippsSortRadixIndexGetBufferSize(len, dataType, pBufSize) bind(C, name='ippsSortRadixIndexGetBufferSize')
     use  iso_c_binding
     integer, intent(in), value :: len
     integer, intent(in), value :: dataType
     integer, intent(out) :: pBufSize
    end
  end interface
  interface
    subroutine ippsSortAscend_32s_I(pSrc, len) bind(C, name='ippsSortAscend_32s_I')
     use  iso_c_binding
     integer, intent(in), value :: len
     integer, intent(inout) :: pSrc(len)
    end
  end interface
  interface
    subroutine ippsSortRadixAscend_32s_I(pSrc, len, pTmp) bind(C, name='ippsSortRadixAscend_32s_I')
     use  iso_c_binding
     integer, intent(in), value :: len
     integer, intent(inout) :: pSrc(len)
     integer, intent(inout) :: pTmp(len)
    end
  end interface
  interface
    subroutine ippsSortRadixIndexAscend_32s(pSrc, srcStrideBytes, pDstIndx, len, pTmpIndx) bind(C, name='ippsSortRadixIndexAscend_32s')
     use  iso_c_binding
     integer, intent(in), value :: len
     integer, intent(inout) :: pSrc(len)
     integer, intent(in), value :: srcStrideBytes
     integer, intent(inout) :: pDstIndx(len)
     integer, intent(inout) :: pTmpIndx(len)
    end
  end interface
  interface
    subroutine ippsSortRadixIndexAscend_32f(pSrc, srcStrideBytes, pDstIndx, len, pTmpIndx) bind(C,name='ippsSortRadixIndexAscend_32f')
     use  iso_c_binding
     integer, intent(in), value :: len
     real   , intent(inout) :: pSrc(len)
     integer, intent(in), value :: srcStrideBytes
     integer, intent(inout) :: pDstIndx(len)
     integer, intent(inout) :: pTmpIndx(len)
    end
  end interface
  interface
    subroutine ippsSortIndexAscend_32f_I(pSrcDst, pDstIndx, len) bind(C,name='ippsSortIndexAscend_32f_I')
     use  iso_c_binding
     real(4), intent(in) :: pSrcDst(*)
     integer(4), intent(inout) :: pDstIndx(*)
     integer(4), intent(in), value :: len
    end
  end interface
  interface
    subroutine ippsSortIndexAscend_32s_I(pSrcDst, pDstIndx, len) bind(C,name='ippsSortIndexAscend_32s_I')
     use  iso_c_binding
     integer(4), intent(in) :: pSrcDst(*)
     integer(4), intent(inout) :: pDstIndx(*)
     integer(4), intent(in), value :: len
    end
  end interface
  interface
    subroutine ippsSortIndexAscend_64f_I(pSrcDst, pDstIndx, len) bind(C,name='ippsSortIndexAscend_64f_I')
     use  iso_c_binding
     real(8), intent(in) :: pSrcDst(*)
     integer(4), intent(inout) :: pDstIndx(*)
     integer(4), intent(in), value :: len
    end
  end interface
end module
