.. _module_basis: 
 
.. program:: basis 
 
.. default-role:: option 
 
======
basis
======

Description of the |AO| basis set.

The structure of the basis set (number of shells and primitives,
angular momenta, exponents, contraction coefficients) is read from the
|EZFIO| database through the EZFIO parameters of this module: `basis`,
`shell_num`, `shell_ang_mom`, `shell_prim_num`, `prim_num`,
`prim_expo` and `prim_coef`. The provider
:c:data:`prim_normalization_factor` is the normalization factor of
each primitive Gaussian. When :option:`basis primitives_normalized`
is |false|, it is set to 1; when it is |true| and the factors are not
stored in the |EZFIO| file, they are computed numerically.
 
 
 
EZFIO parameters 
---------------- 
 
.. option:: basis
 
    Name of the |AO| basis set
 
 
.. option:: basis_type
 
    Type of basis set. (Gaussian|Slater|Numerical)
 
 
.. option:: slater_r_power
 
    Power of r in Slater functions
 
 
.. option:: shell_num
 
    Number of shells
 
 
.. option:: nucleus_shell_num
 
    Number of shells per nucleus
 
 
.. option:: shell_normalization_factor
 
    Normalization factor applied to the whole shell, ex $1/\sqrt{ <d_{z^2}|d_{z^2}>}$
 
 
.. option:: shell_ang_mom
 
    Angular momentum of each shell
 
 
.. option:: shell_prim_num
 
    Number of primitives in a shell
 
 
.. option:: shell_index
 
    Index of the shell for each primitive
 
 
.. option:: basis_nucleus_index
 
    Nucleus on which the shell is centered
 
 
.. option:: prim_normalization_factor
 
    Normalization factor applied to each primitive
 
 
.. option:: prim_num
 
    Total number of primitives
 
 
.. option:: prim_coef
 
    Primitive coefficients
 
 
.. option:: prim_expo
 
    Exponents in the shell
 
 
.. option:: primitives_normalized
 
    If true, assume primitive basis functions are normalized
 
    Default: true
 
.. option:: ao_normalized
 
    If true, normalize the basis functions
 
    Default: false
 
 
Providers 
--------- 
 
.. c:var:: double prim_normalization_factor


    File : :file:`basis/basis.irp.f`

    .. code:: fortran

        double precision, allocatable	:: prim_normalization_factor	(prim_num)


    Number of primitives per |AO|

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`mpi_master`
       * :c:data:`output_wall_time_0`
       * :c:data:`prim_expo`
       * :c:data:`prim_num`
       * :c:data:`primitives_normalized`
       * :c:data:`shell_ang_mom`
       * :c:data:`shell_index`
       * :c:data:`shell_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized`

