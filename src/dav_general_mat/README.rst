===============
dav_general_mat
===============

Davidson diagonalization of a general symmetric matrix.

Unlike the Davidson diagonalization of the Hamiltonian (see the
:ref:`module_davidson` module), this module diagonalizes any symmetric
matrix, with two ways to provide the matrix-vector product:

* :c:func:`davidson_general` takes the matrix directly
* :c:func:`davidson_general_ext_rout` takes a subroutine that
  computes :math:`v = H u`

Both have dressed variants, :c:func:`davidson_general_ext_rout_dressed`
and :c:func:`davidson_general_ext_rout_diag_dressed`, that add the
dressing operator of the :ref:`module_dressing` module, as used in
the |CIPSI| iteration, and non symmetric variants that do not assume
the symmetry of the matrix.

The program :ref:`test_dav` is a small example of the use of these
routines (see :file:`test_dav.irp.f`).
