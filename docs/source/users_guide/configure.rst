.. _configure:

=========
configure
=========

.. program:: configure

The :command:`configure` script prepares the compilation of |qp|. It
initializes the git submodules, checks that the external dependencies
are installed, and generates the |Ninja| build file
(:file:`config/build.ninja`) from a compilation configuration file.

Usage
=====

.. code:: bash

   ./configure
   ./configure -c <file>
   ./configure -c <file> -g <gpu>
   ./configure -i <package>
   ./configure -h

Options
-------

.. option:: -c <file>

   Define the COMPILATION configuration file, in
   :file:`${QP_ROOT}/config/`. This generates the
   :file:`config/build.ninja` file that drives the compilation. See
   "Choosing the compilation configuration file" below.

.. option:: -g [nvidia|intel|none]

   Choose the GPU acceleration. This creates the
   :file:`src/gpu_arch` symbolic link to the matching plugin in
   :file:`plugins/local/`. Without this option, the link points to the
   :file:`gpu_x86` plugin, the CPU implementation of the GPU routines.
   See the GPU acceleration section below.

.. option:: -i <package>

   Install <package>. Use at your own risk: no support will be provided
   for the installation of external dependencies. The supported packages
   and the manual fallbacks are described on the installation page.

.. option:: -h

   Print the help message.

.. note::

   Running :command:`configure` needs a network connection: it
   initializes the git submodules and updates the
   :file:`external/qp2-dependencies` repository. If the
   :command:`irpf90` executable is not found in the :envvar:`PATH`,
   configure builds it from the :file:`external/irpf90` submodule.

GPU acceleration
================

GPU acceleration works well with |CCSD|, but for |CI| workflows it is
usually not worth it. It becomes worthwhile on machines where the CPU
is not very powerful, for example a laptop with an Nvidia GPU: using
the Cholesky vectors in single precision can then show a performance
boost. Do not run with the GPU path on supercomputers; it is a waste
of resources. Prefer the highly efficient multicore CPUs.

Choosing the compilation configuration file
===========================================

Each file in the :file:`${QP_ROOT}/config/` directory describes how to
compile |qp| on a given machine: the Fortran and C compilers
(:samp:`FC` and :samp:`CC`), the BLAS/LAPACK libraries
(:samp:`LAPACK_LIB`), and the |IRPF90| options
(:samp:`IRPF90` and :samp:`IRPF90_FLAGS`) in the ``[COMMON]`` section.
The :samp:`MODE` option of the ``[OPTION]`` section selects one of the
``[OPT]``, ``[PROFILE]`` or ``[DEBUG]`` sections, and the ``[OPENMP]``
section is appended when :samp:`OPENMP` is 1. The flags of the selected
sections are concatenated in this order.

Pick the file that matches your machine:

Your machine or situation                                  Configuration file
=========================================================== ===================================
x86-64 Linux, GNU compilers, system BLAS/LAPACK (default)   gfortran.cfg
x86-64 Linux, AVX enabled                                   gfortran_avx.cfg
x86-64 Linux, AVX2 enabled                                  gfortran_avx2.cfg
x86-64 Linux, OpenBLAS instead of the system BLAS           gfortran_openblas.cfg
x86-64 Linux, AVX2 and OpenBLAS                             gfortran_avx2_openblas.cfg
Intel MKL installed, recommended (set MKLROOT)              gfortran_mkl.cfg
Intel MKL installed, debug build (set MKLROOT)              gfortran_debug_mkl.cfg
MPI-parallel build (mpif90)                                 gfortran_mpi.cfg
Debug build (MODE: DEBUG)                                   gfortran_debug.cfg
Pinned to gfortran-10, debug build                          gfortran10.cfg
macOS (Homebrew gfortran), alpha                            gfortran_macos.cfg
ARM Linux with the ARMPL BLAS                               gfortran_armpl.cfg
ThreadSanitizer data-race debugging                         gfortran_avx_thread_sanitizer.cfg
LLVM flang                                                  flang_avx.cfg
NVIDIA HPC SDK (not recommended)                            nvfortran.cfg
Intel Fortran (ifort) clusters (not recommended)            ifort_2019_*.cfg, ifort_2021_*.cfg
Cray, GNU compilers                                         cray_gfortran.cfg
Cray, Intel compilers                                       cray_intel.cfg
Bull clusters                                               bull.cfg
Continuous integration only (not for users)                 travis.cfg

The file names follow the pattern
``<compiler>[_<instruction-set>][_<blas>][_mpi][_debug]``:
:file:`gfortran_avx2_openblas.cfg` uses gfortran with AVX2 and OpenBLAS,
and the ifort files carry the compiler year and the node instruction set
(:file:`ifort_2021_rome.cfg` is ifort 2021 tuned for AMD Rome nodes).
The :samp:`_mpi` variants use the MPI wrappers, the :samp:`_debug`
variants set :samp:`MODE : DEBUG`, and the :samp:`_notz` variants define
``WITHOUT_TRAILZ`` and ``WITHOUT_SHIFTRL``, which disable the code paths
based on those instructions. All the ifort configurations
link Intel MKL.

macOS is not yet supported: the :file:`gfortran_macos.cfg` configuration
is there because Mac support is in development, but it is still in alpha
state.

The :samp:`popcnt` instruction accelerates the programs a lot, and the
performance gain from enabling the SSE4.2, AVX or AVX2 instruction sets
is large. Check which instruction sets your CPU supports before
choosing a configuration file:

.. code:: bash

   lscpu | grep -o 'avx2\|avx' | sort -u

The default :file:`gfortran.cfg` configuration works on any machine,
but the performance of the resulting binaries is poor. The recommended
setup for x86-64 machines is gfortran with Intel MKL
(:file:`gfortran_mkl.cfg`). The Intel Fortran (:command:`ifx`) and
NVIDIA Fortran (:command:`nvfortran`) configurations are not
recommended: :command:`ifx` gives poor performance, and
:command:`nvfortran` gave a bad experience in 2021 and has not been
tested since.

MPI is not encouraged by default. The :samp:`_mpi` configurations only
accelerate very large CIPSI runs: the memory is not distributed across the
nodes, and single-node runs are fully OpenMP. MPI runs with |qp| are
exceptional, and reserved for a few particular use cases.

Do not modify the shipped configuration files
=============================================

The configuration files shipped in :file:`config/` are tracked by git.
Modifying one of them in place means that the next
:command:`git pull` will reject your changes or produce merge
conflicts.

Copy the file that matches your machine, and modify the copy:

.. code:: bash

   cp config/gfortran.cfg config/gfortran_mybox.cfg
   # edit config/gfortran_mybox.cfg
   ./configure -c config/gfortran_mybox.cfg

New configuration files are ignored by git (:file:`config/.gitignore`
contains :samp:`*.cfg`), so your copy will not dirty the
:command:`git status` output and cannot conflict with the updates you
pull. The shipped files remain tracked and keep working as they are.

Example
=======

.. code:: bash

   ./configure
   ./configure -c config/gfortran.cfg
   ./configure -c config/gfortran_mybox.cfg -g nvidia
