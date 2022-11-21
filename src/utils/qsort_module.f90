module qsort_module
  use iso_c_binding
  
  interface
     
     subroutine i2sort_c(A, iorder, isize) bind(C, name="qsort_int16_t")
       use iso_c_binding
       integer(c_int32_t), value :: isize
       integer(c_int32_t)        :: iorder(isize)
       integer   (c_int16_t)         :: A(isize)
     end subroutine i2sort_c
     
     subroutine i2sort_noidx_c(A, isize) bind(C, name="qsort_int16_t_noidx")
       use iso_c_binding
       integer(c_int32_t), value :: isize
       integer   (c_int16_t)         :: A(isize)
     end subroutine i2sort_noidx_c
     
     
     
     subroutine i2sort_big_c(A, iorder, isize) bind(C, name="qsort_int16_t_big")
       use iso_c_binding
       integer(c_int64_t), value :: isize
       integer(c_int64_t)        :: iorder(isize)
       integer   (c_int16_t)         :: A(isize)
     end subroutine i2sort_big_c
     
     subroutine i2sort_noidx_big_c(A, isize) bind(C, name="qsort_int16_t_noidx_big")
       use iso_c_binding
       integer(c_int64_t), value :: isize
       integer   (c_int16_t)         :: A(isize)
     end subroutine i2sort_noidx_big_c
     
     
     
     subroutine isort_c(A, iorder, isize) bind(C, name="qsort_int32_t")
       use iso_c_binding
       integer(c_int32_t), value :: isize
       integer(c_int32_t)        :: iorder(isize)
       integer   (c_int32_t)         :: A(isize)
     end subroutine isort_c
     
     subroutine isort_noidx_c(A, isize) bind(C, name="qsort_int32_t_noidx")
       use iso_c_binding
       integer(c_int32_t), value :: isize
       integer   (c_int32_t)         :: A(isize)
     end subroutine isort_noidx_c
     
     
     
     subroutine isort_big_c(A, iorder, isize) bind(C, name="qsort_int32_t_big")
       use iso_c_binding
       integer(c_int64_t), value :: isize
       integer(c_int64_t)        :: iorder(isize)
       integer   (c_int32_t)         :: A(isize)
     end subroutine isort_big_c
     
     subroutine isort_noidx_big_c(A, isize) bind(C, name="qsort_int32_t_noidx_big")
       use iso_c_binding
       integer(c_int64_t), value :: isize
       integer   (c_int32_t)         :: A(isize)
     end subroutine isort_noidx_big_c
     
     
     
     subroutine i8sort_c(A, iorder, isize) bind(C, name="qsort_int64_t")
       use iso_c_binding
       integer(c_int32_t), value :: isize
       integer(c_int32_t)        :: iorder(isize)
       integer   (c_int64_t)         :: A(isize)
     end subroutine i8sort_c
     
     subroutine i8sort_noidx_c(A, isize) bind(C, name="qsort_int64_t_noidx")
       use iso_c_binding
       integer(c_int32_t), value :: isize
       integer   (c_int64_t)         :: A(isize)
     end subroutine i8sort_noidx_c
     
     
     
     subroutine i8sort_big_c(A, iorder, isize) bind(C, name="qsort_int64_t_big")
       use iso_c_binding
       integer(c_int64_t), value :: isize
       integer(c_int64_t)        :: iorder(isize)
       integer   (c_int64_t)         :: A(isize)
     end subroutine i8sort_big_c
     
     subroutine i8sort_noidx_big_c(A, isize) bind(C, name="qsort_int64_t_noidx_big")
       use iso_c_binding
       integer(c_int64_t), value :: isize
       integer   (c_int64_t)         :: A(isize)
     end subroutine i8sort_noidx_big_c
     
     
     
     subroutine dsort_c(A, iorder, isize) bind(C, name="qsort_double")
       use iso_c_binding
       integer(c_int32_t), value :: isize
       integer(c_int32_t)        :: iorder(isize)
       real   (c_double)         :: A(isize)
     end subroutine dsort_c
     
     subroutine dsort_noidx_c(A, isize) bind(C, name="qsort_double_noidx")
       use iso_c_binding
       integer(c_int32_t), value :: isize
       real   (c_double)         :: A(isize)
     end subroutine dsort_noidx_c
     
     
     
     subroutine dsort_big_c(A, iorder, isize) bind(C, name="qsort_double_big")
       use iso_c_binding
       integer(c_int64_t), value :: isize
       integer(c_int64_t)        :: iorder(isize)
       real   (c_double)         :: A(isize)
     end subroutine dsort_big_c
     
     subroutine dsort_noidx_big_c(A, isize) bind(C, name="qsort_double_noidx_big")
       use iso_c_binding
       integer(c_int64_t), value :: isize
       real   (c_double)         :: A(isize)
     end subroutine dsort_noidx_big_c
     
     
     
     subroutine sort_c(A, iorder, isize) bind(C, name="qsort_float")
       use iso_c_binding
       integer(c_int32_t), value :: isize
       integer(c_int32_t)        :: iorder(isize)
       real   (c_float)         :: A(isize)
     end subroutine sort_c
     
     subroutine sort_noidx_c(A, isize) bind(C, name="qsort_float_noidx")
       use iso_c_binding
       integer(c_int32_t), value :: isize
       real   (c_float)         :: A(isize)
     end subroutine sort_noidx_c
     
     
     
     subroutine sort_big_c(A, iorder, isize) bind(C, name="qsort_float_big")
       use iso_c_binding
       integer(c_int64_t), value :: isize
       integer(c_int64_t)        :: iorder(isize)
       real   (c_float)         :: A(isize)
     end subroutine sort_big_c
     
     subroutine sort_noidx_big_c(A, isize) bind(C, name="qsort_float_noidx_big")
       use iso_c_binding
       integer(c_int64_t), value :: isize
       real   (c_float)         :: A(isize)
     end subroutine sort_noidx_big_c
     
     
     
  end interface

end module qsort_module


subroutine i2sort(A, iorder, isize) 
  use qsort_module
  use iso_c_binding
  integer(c_int32_t)        :: isize
  integer(c_int32_t)        :: iorder(isize)
  integer   (c_int16_t)         :: A(isize)
  call i2sort_c(A, iorder, isize)
end subroutine i2sort

subroutine i2sort_noidx(A, isize)
  use iso_c_binding
  use qsort_module
  integer(c_int32_t) :: isize
  integer   (c_int16_t)    :: A(isize)
  call i2sort_noidx_c(A, isize)
end subroutine i2sort_noidx



subroutine i2sort_big(A, iorder, isize) 
  use qsort_module
  use iso_c_binding
  integer(c_int64_t)        :: isize
  integer(c_int64_t)        :: iorder(isize)
  integer   (c_int16_t)         :: A(isize)
  call i2sort_big_c(A, iorder, isize)
end subroutine i2sort_big

subroutine i2sort_noidx_big(A, isize)
  use iso_c_binding
  use qsort_module
  integer(c_int64_t) :: isize
  integer   (c_int16_t)    :: A(isize)
  call i2sort_noidx_big_c(A, isize)
end subroutine i2sort_noidx_big



subroutine isort(A, iorder, isize) 
  use qsort_module
  use iso_c_binding
  integer(c_int32_t)        :: isize
  integer(c_int32_t)        :: iorder(isize)
  integer   (c_int32_t)         :: A(isize)
  call isort_c(A, iorder, isize)
end subroutine isort

subroutine isort_noidx(A, isize)
  use iso_c_binding
  use qsort_module
  integer(c_int32_t) :: isize
  integer   (c_int32_t)    :: A(isize)
  call isort_noidx_c(A, isize)
end subroutine isort_noidx



subroutine isort_big(A, iorder, isize) 
  use qsort_module
  use iso_c_binding
  integer(c_int64_t)        :: isize
  integer(c_int64_t)        :: iorder(isize)
  integer   (c_int32_t)         :: A(isize)
  call isort_big_c(A, iorder, isize)
end subroutine isort_big

subroutine isort_noidx_big(A, isize)
  use iso_c_binding
  use qsort_module
  integer(c_int64_t) :: isize
  integer   (c_int32_t)    :: A(isize)
  call isort_noidx_big_c(A, isize)
end subroutine isort_noidx_big



subroutine i8sort(A, iorder, isize) 
  use qsort_module
  use iso_c_binding
  integer(c_int32_t)        :: isize
  integer(c_int32_t)        :: iorder(isize)
  integer   (c_int64_t)         :: A(isize)
  call i8sort_c(A, iorder, isize)
end subroutine i8sort

subroutine i8sort_noidx(A, isize)
  use iso_c_binding
  use qsort_module
  integer(c_int32_t) :: isize
  integer   (c_int64_t)    :: A(isize)
  call i8sort_noidx_c(A, isize)
end subroutine i8sort_noidx



subroutine i8sort_big(A, iorder, isize) 
  use qsort_module
  use iso_c_binding
  integer(c_int64_t)        :: isize
  integer(c_int64_t)        :: iorder(isize)
  integer   (c_int64_t)         :: A(isize)
  call i8sort_big_c(A, iorder, isize)
end subroutine i8sort_big

subroutine i8sort_noidx_big(A, isize)
  use iso_c_binding
  use qsort_module
  integer(c_int64_t) :: isize
  integer   (c_int64_t)    :: A(isize)
  call i8sort_noidx_big_c(A, isize)
end subroutine i8sort_noidx_big



subroutine dsort(A, iorder, isize) 
  use qsort_module
  use iso_c_binding
  integer(c_int32_t)        :: isize
  integer(c_int32_t)        :: iorder(isize)
  real   (c_double)         :: A(isize)
  call dsort_c(A, iorder, isize)
end subroutine dsort

subroutine dsort_noidx(A, isize)
  use iso_c_binding
  use qsort_module
  integer(c_int32_t) :: isize
  real   (c_double)    :: A(isize)
  call dsort_noidx_c(A, isize)
end subroutine dsort_noidx



subroutine dsort_big(A, iorder, isize) 
  use qsort_module
  use iso_c_binding
  integer(c_int64_t)        :: isize
  integer(c_int64_t)        :: iorder(isize)
  real   (c_double)         :: A(isize)
  call dsort_big_c(A, iorder, isize)
end subroutine dsort_big

subroutine dsort_noidx_big(A, isize)
  use iso_c_binding
  use qsort_module
  integer(c_int64_t) :: isize
  real   (c_double)    :: A(isize)
  call dsort_noidx_big_c(A, isize)
end subroutine dsort_noidx_big



subroutine sort(A, iorder, isize) 
  use qsort_module
  use iso_c_binding
  integer(c_int32_t)        :: isize
  integer(c_int32_t)        :: iorder(isize)
  real   (c_float)         :: A(isize)
  call sort_c(A, iorder, isize)
end subroutine sort

subroutine sort_noidx(A, isize)
  use iso_c_binding
  use qsort_module
  integer(c_int32_t) :: isize
  real   (c_float)    :: A(isize)
  call sort_noidx_c(A, isize)
end subroutine sort_noidx



subroutine sort_big(A, iorder, isize) 
  use qsort_module
  use iso_c_binding
  integer(c_int64_t)        :: isize
  integer(c_int64_t)        :: iorder(isize)
  real   (c_float)         :: A(isize)
  call sort_big_c(A, iorder, isize)
end subroutine sort_big

subroutine sort_noidx_big(A, isize)
  use iso_c_binding
  use qsort_module
  integer(c_int64_t) :: isize
  real   (c_float)    :: A(isize)
  call sort_noidx_big_c(A, isize)
end subroutine sort_noidx_big
