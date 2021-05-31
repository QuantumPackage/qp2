module intel
  use, intrinsic :: iso_c_binding
  interface
    subroutine ippsSortAscend_32s_I(pSrc, len) bind(C, name='ippsSortAscend_32s_I')
     use  iso_c_binding
     integer, intent(in), value :: len
     integer, intent(inout) :: pSrc(len)
    end
  end interface
  interface
    subroutine ippsSortAscend_32f_I(pSrc, len) bind(C, name='ippsSortAscend_32f_I')
     use  iso_c_binding
     integer, intent(in), value :: len
     real, intent(inout) :: pSrc(len)
    end
  end interface
  interface
    subroutine ippsSortAscend_64s_I(pSrc, len) bind(C, name='ippsSortAscend_64s_I')
     use  iso_c_binding
     integer*8, intent(in), value :: len
     integer, intent(inout) :: pSrc(len)
    end
  end interface
  interface
    subroutine ippsSortAscend_64f_I(pSrc, len) bind(C, name='ippsSortAscend_64f_I')
     use  iso_c_binding
     double precision, intent(in), value :: len
     real, intent(inout) :: pSrc(len)
    end
  end interface

  interface
    subroutine ippsSortRadixIndexGetBufferSize(len, dataType, pBufSize) bind(C, name='ippsSortRadixIndexGetBufferSize')
     use  iso_c_binding
     integer, intent(in), value :: len
     integer, intent(in), value :: dataType
     integer, intent(out) :: pBufSize
    end
  end interface

  interface
    subroutine ippsSortRadixAscend_16s_I(pSrc, len, pTmp) bind(C, name='ippsSortRadixAscend_16s_I')
     use  iso_c_binding
     integer, intent(in), value :: len
     integer*2, intent(inout) :: pSrc(len)
     character, intent(inout) :: pTmp(len)
    end
  end interface
  interface
    subroutine ippsSortRadixAscend_32s_I(pSrc, len, pTmp) bind(C, name='ippsSortRadixAscend_32s_I')
     use  iso_c_binding
     integer, intent(in), value :: len
     integer, intent(inout) :: pSrc(len)
     character, intent(inout) :: pTmp(len)
    end
  end interface
  interface
    subroutine ippsSortRadixAscend_32f_I(pSrc, len, pTmp) bind(C, name='ippsSortRadixAscend_32f_I')
     use  iso_c_binding
     integer, intent(in), value :: len
     real, intent(inout) :: pSrc(len)
     character, intent(inout) :: pTmp(len)
    end
  end interface
  interface
    subroutine ippsSortRadixAscend_64s_I(pSrc, len, pTmp) bind(C, name='ippsSortRadixAscend_64s_I')
     use  iso_c_binding
     integer, intent(in), value :: len
     integer*8, intent(inout) :: pSrc(len)
     character, intent(inout) :: pTmp(len)
    end
  end interface
  interface
    subroutine ippsSortRadixAscend_64f_I(pSrc, len, pTmp) bind(C, name='ippsSortRadixAscend_64f_I')
     use  iso_c_binding
     integer, intent(in), value :: len
     double precision, intent(inout) :: pSrc(len)
     character, intent(inout) :: pTmp(len)
    end
  end interface

  interface
    subroutine ippsSortRadixIndexAscend_16s(pSrc, srcStrideBytes, pDstIndx, len, pTmpIndx) bind(C, name='ippsSortRadixIndexAscend_16s')
     use  iso_c_binding
     integer, intent(in), value :: len
     integer*2, intent(inout) :: pSrc(len)
     integer, intent(in), value :: srcStrideBytes
     integer, intent(inout) :: pDstIndx(len)
     character, intent(inout) :: pTmpIndx(len)
    end
  end interface
  interface
    subroutine ippsSortRadixIndexAscend_32s(pSrc, srcStrideBytes, pDstIndx, len, pTmpIndx) bind(C, name='ippsSortRadixIndexAscend_32s')
     use  iso_c_binding
     integer, intent(in), value :: len
     integer, intent(inout) :: pSrc(len)
     integer, intent(in), value :: srcStrideBytes
     integer, intent(inout) :: pDstIndx(len)
     character, intent(inout) :: pTmpIndx(len)
    end
  end interface
  interface
    subroutine ippsSortRadixIndexAscend_32f(pSrc, srcStrideBytes, pDstIndx, len, pTmpIndx) bind(C,name='ippsSortRadixIndexAscend_32f')
     use  iso_c_binding
     integer, intent(in), value :: len
     real   , intent(inout) :: pSrc(len)
     integer, intent(in), value :: srcStrideBytes
     integer, intent(inout) :: pDstIndx(len)
     character, intent(inout) :: pTmpIndx(len)
    end
  end interface
  interface
    subroutine ippsSortRadixIndexAscend_64s(pSrc, srcStrideBytes, pDstIndx, len, pTmpIndx) bind(C, name='ippsSortRadixIndexAscend_64s')
     use  iso_c_binding
     integer, intent(in), value :: len
     integer*8, intent(inout) :: pSrc(len)
     integer, intent(in), value :: srcStrideBytes
     integer, intent(inout) :: pDstIndx(len)
     character, intent(inout) :: pTmpIndx(len)
    end
  end interface
  interface
    subroutine ippsSortRadixIndexAscend_64f(pSrc, srcStrideBytes, pDstIndx, len, pTmpIndx) bind(C,name='ippsSortRadixIndexAscend_64f')
     use  iso_c_binding
     integer, intent(in), value :: len
     real*8 , intent(inout) :: pSrc(len)
     integer, intent(in), value :: srcStrideBytes
     integer, intent(inout) :: pDstIndx(len)
     character, intent(inout) :: pTmpIndx(len)
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
  interface
    subroutine ippsSortIndexAscend_64s_I(pSrcDst, pDstIndx, len) bind(C,name='ippsSortIndexAscend_64s_I')
     use  iso_c_binding
     integer(8), intent(in) :: pSrcDst(*)
     integer(4), intent(inout) :: pDstIndx(*)
     integer(4), intent(in), value :: len
    end
  end interface
  interface
    subroutine ippsSortIndexAscend_16s_I(pSrcDst, pDstIndx, len) bind(C,name='ippsSortIndexAscend_16s_I')
     use  iso_c_binding
     integer(2), intent(in) :: pSrcDst(*)
     integer(4), intent(inout) :: pDstIndx(*)
     integer(4), intent(in), value :: len
    end
  end interface
end module
