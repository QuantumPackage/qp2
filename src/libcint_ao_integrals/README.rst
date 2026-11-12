====================
libcint_ao_integrals
====================

|AO| integrals computed by `libcint <https://github.com/sunqm/libcint>`_.

This module is an alternative engine for the overlap, kinetic energy,
nucleus-electron attraction and two-electron repulsion integrals. It is
selected at run time by :c:data:`use_libcint`, which defaults to
true; when it is false, or when the basis cannot be represented in
libcint, every integral is computed by the built-in |QP| code and the
results are unchanged.

The choice is a property of the |EZFIO| database, not of the binary, so
the same :ref:`qp_run` executable serves both engines and a calculation
can be replayed with the other engine by flipping one logical.

What libcint cannot represent
=============================

The following cases fall back to the built-in code, with a message on the
master rank:

* complex |AOs| (:c:data:`use_cgtos` is true);
* periodic systems (:c:data:`is_periodic` is true);
* angular momenta above 12;
* a basis in which the |AOs| of a shell do not share the same contraction
  (|QP| normalizes each |AO| separately, libcint has one coefficient set
  per shell). This rejects a TREXIO file whose primitive normalization
  depends on the AO powers (per-component normalization);
* a basis in which an |AO| is not on the shell's center, or whose
  Cartesian powers are not in libcint's order (x, y, z, xx, xy, ...).
  Both are stored per |AO|, so libcint would silently return permuted
  blocks;
* non-integer or negative nuclear charges, for the nucleus-electron
  operator only: the two-electron integrals stay on libcint. libcint
  replaces the charge by :math:`-|\text{charge}|`, so a negative charge
  would change sign.

Coefficient convention
======================

libcint multiplies every integral by :math:`f(l_i) f(l_j)` (one electron)
or :math:`f(l_i) f(l_j) f(l_k) f(l_l)` (two electrons), where
:math:`f(0) = 0.282094791773878143`, :math:`f(1) = 0.488602511902919921`
and :math:`f(l \ge 2) = 1` (`CINTcommon_fac_sp`). The contraction
coefficients stored in ``env`` are divided by that factor, so the
integrals come out in the |QP| convention. The remaining per-|AO|
normalization, which libcint cannot express, is applied to each block as
:c:data:`libcint_scale`.

Building
========

libcint must be built with::

  ./configure -i libcint

which pins v6.2.1, builds a static `libcint.a` with ``I8=off`` and
installs it in ``$QP_ROOT/lib``. ``I8=on`` changes the width of every
``FINT`` argument and silently corrupts every integral: do not link
against a system libcint.

The design of the mapping is written out in :file:`docs/libcint_design.md`.

Tests
=====

:ref:`test_libcint` prints every |AO| integral tagged by operator, so
that a run with :c:data:`use_libcint` true and a run with false can be
compared: the two must agree to :math:`10^{-10}` absolute for the
overlaps (QP's built-in overlap has a hard-coded :math:`10^{-10}` floor
when a primitive pair underflows) and :math:`10^{-11}` for the other
operators.
