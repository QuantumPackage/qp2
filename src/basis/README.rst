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
