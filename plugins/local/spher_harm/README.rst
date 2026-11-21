==========
spher_harm
==========

Evaluation of the spherical harmonics :math:`Y_{lm}` in real space.

The main routine is `spher_harm_func_r3(r, l, m, re_ylm, im_ylm)`,
which evaluates :math:`Y_{lm}` at a point :math:`r` of the real
space. The evaluation is done in spherical coordinates, with the
associated Legendre polynomials `plgndr` (`assoc_gaus_pol.irp.f`).
The test routine `test_spher_harm` (`routines_test.irp.f`)
documents and exercises the different variants of the computation.
