#include <simple_gpu.F90>

module gpu_affinity_mod
  use iso_c_binding
  use omp_lib
  implicit none

  interface
    subroutine build_gpu_affinity_map() &
        bind(C, name="build_gpu_affinity_map_")
    end subroutine

    integer(c_int) function get_closest_gpu() &
        bind(C, name="get_closest_gpu_")
      use iso_c_binding
    end function
  end interface

end module

