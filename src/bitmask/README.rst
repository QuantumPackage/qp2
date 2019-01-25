==============
bitmask module
==============

The central part of this module is the :file:`bitmasks_module.f90` file. It contains
the constants that will be used to define on which kind of integer the bitmasks
will be defined.

In the program, to represent a determinant as a pair of bitstrings,
the determinant should be defined as

.. code-block:: fortran

  use bitmasks
  integer(bit_kind)  :: determinant(N_int,2)


:file:`bitmasks_routines.irp.f` contains helper routines to manipulate bitmask, like
transforming a bit string to a list of integers for example.


`bit_kind_shift`, `bit_kind_size` and `bit_kind` are supposed to be consistent::

   2**bit_kind_shift = bit_kind_size
   bit_kind = bit_kind_size / 8


For an example of how to use the bitmaks, see the file :file:`example.irp.f`.
