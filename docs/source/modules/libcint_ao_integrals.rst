.. _module_libcint_ao_integrals: 
 
.. program:: libcint_ao_integrals 
 
.. default-role:: option 
 
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
  per shell). This rejects files written with the legacy GAMESS
  convention, `ezfio_convention < 20250211`;
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
 
 
 
Providers 
--------- 
 
.. c:var:: int libcint_atm


    File : :file:`libcint_ao_integrals/libcint_data.irp.f`

    .. code:: fortran

        integer, allocatable	:: libcint_atm	(ATM_SLOTS*nucl_num)
        integer, allocatable	:: libcint_bas	(BAS_SLOTS*shell_num)
        double precision, allocatable	:: libcint_env	(libcint_nenv)


    libcint description of the |AO| basis.
    
    - `atm` has one point-charge entry per nucleus;
    - `bas` has one shell per |QP| shell, with `NCTR_OF = 1`;
    - `env` holds the coordinates, the exponents and the contraction
      coefficients, the latter divided by `CINTcommon_fac_sp(l)` so that
      libcint returns the raw integrals in the |QP| convention. Only the
      coefficients of the first |AO| of each shell are stored; the
      per-|AO| normalization is carried by :f:func:`libcint_scale`.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized`
       * :c:data:`ao_coef_normalized_ordered`
       * :c:data:`ao_coef_normalized_ordered`
       * :c:data:`ao_first_of_shell`
       * :c:data:`basis_nucleus_index`
       * :c:data:`libcint_nenv`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`
       * :c:data:`shell_ang_mom`
       * :c:data:`shell_num`
       * :c:data:`shell_prim_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`libcint_available`
       * :c:data:`libcint_env_lr`

 
.. c:var:: bool libcint_available


    File : :file:`libcint_ao_integrals/libcint_data.irp.f`

    .. code:: fortran

        logical	:: libcint_available	


    True if the |AO| basis and the operators it is used with can be
    represented in libcint. If false, every integral falls back to the
    built-in |QP| code.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized`
       * :c:data:`ao_coef_normalized`
       * :c:data:`ao_coef_normalized_ordered`
       * :c:data:`ao_coef_normalized_ordered`
       * :c:data:`ao_first_of_shell`
       * :c:data:`ao_l`
       * :c:data:`ao_l`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`ao_prim_num_max`
       * :c:data:`basis_nucleus_index`
       * :c:data:`is_periodic`
       * :c:data:`libcint_atm`
       * :c:data:`libcint_env_lr`
       * :c:data:`libcint_scale`
       * :c:data:`mpi_master`
       * :c:data:`nucl_num`
       * :c:data:`prim_num`
       * :c:data:`shell_ang_mom`
       * :c:data:`shell_num`
       * :c:data:`shell_prim_num`
       * :c:data:`use_cgtos`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`libcint_nuc_available`
       * :c:data:`use_libcint_effective`

 
.. c:var:: int libcint_bas


    File : :file:`libcint_ao_integrals/libcint_data.irp.f`

    .. code:: fortran

        integer, allocatable	:: libcint_atm	(ATM_SLOTS*nucl_num)
        integer, allocatable	:: libcint_bas	(BAS_SLOTS*shell_num)
        double precision, allocatable	:: libcint_env	(libcint_nenv)


    libcint description of the |AO| basis.
    
    - `atm` has one point-charge entry per nucleus;
    - `bas` has one shell per |QP| shell, with `NCTR_OF = 1`;
    - `env` holds the coordinates, the exponents and the contraction
      coefficients, the latter divided by `CINTcommon_fac_sp(l)` so that
      libcint returns the raw integrals in the |QP| convention. Only the
      coefficients of the first |AO| of each shell are stored; the
      per-|AO| normalization is carried by :f:func:`libcint_scale`.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized`
       * :c:data:`ao_coef_normalized_ordered`
       * :c:data:`ao_coef_normalized_ordered`
       * :c:data:`ao_first_of_shell`
       * :c:data:`basis_nucleus_index`
       * :c:data:`libcint_nenv`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`
       * :c:data:`shell_ang_mom`
       * :c:data:`shell_num`
       * :c:data:`shell_prim_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`libcint_available`
       * :c:data:`libcint_env_lr`

 
.. c:var:: double libcint_env


    File : :file:`libcint_ao_integrals/libcint_data.irp.f`

    .. code:: fortran

        integer, allocatable	:: libcint_atm	(ATM_SLOTS*nucl_num)
        integer, allocatable	:: libcint_bas	(BAS_SLOTS*shell_num)
        double precision, allocatable	:: libcint_env	(libcint_nenv)


    libcint description of the |AO| basis.
    
    - `atm` has one point-charge entry per nucleus;
    - `bas` has one shell per |QP| shell, with `NCTR_OF = 1`;
    - `env` holds the coordinates, the exponents and the contraction
      coefficients, the latter divided by `CINTcommon_fac_sp(l)` so that
      libcint returns the raw integrals in the |QP| convention. Only the
      coefficients of the first |AO| of each shell are stored; the
      per-|AO| normalization is carried by :f:func:`libcint_scale`.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized`
       * :c:data:`ao_coef_normalized_ordered`
       * :c:data:`ao_coef_normalized_ordered`
       * :c:data:`ao_first_of_shell`
       * :c:data:`basis_nucleus_index`
       * :c:data:`libcint_nenv`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`
       * :c:data:`shell_ang_mom`
       * :c:data:`shell_num`
       * :c:data:`shell_prim_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`libcint_available`
       * :c:data:`libcint_env_lr`

 
.. c:var:: double libcint_env_lr


    File : :file:`libcint_ao_integrals/libcint_data.irp.f`

    .. code:: fortran

        double precision, allocatable	:: libcint_env_lr	(libcint_nenv)


    Environment for the :math:`\mathrm{erf}(\mu r_{12})/r_{12}` operator.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`libcint_atm`
       * :c:data:`libcint_atm`
       * :c:data:`libcint_nenv`
       * :c:data:`mu_erf`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`libcint_available`

 
.. c:var:: int libcint_natm


    File : :file:`libcint_ao_integrals/libcint_data.irp.f`

    .. code:: fortran

        integer	:: libcint_natm	


    Number of entries in the :c:data:`libcint_atm` table, i.e. the number
    of nuclei.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_num`


 
.. c:var:: int libcint_nbas


    File : :file:`libcint_ao_integrals/libcint_data.irp.f`

    .. code:: fortran

        integer	:: libcint_nbas	


    Number of entries in the :c:data:`libcint_bas` table, i.e. the number
    of |QP| shells.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`shell_num`


 
.. c:var:: int libcint_nenv


    File : :file:`libcint_ao_integrals/libcint_data.irp.f`

    .. code:: fortran

        integer	:: libcint_nenv	


    Size of the :c:data:`env` array passed to libcint: 20 header slots,
    3 coordinates per nucleus, and the exponents and coefficients of every
    shell.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_num`
       * :c:data:`prim_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`libcint_atm`
       * :c:data:`libcint_env_lr`

 
.. c:var:: bool libcint_nuc_available


    File : :file:`libcint_ao_integrals/libcint_data.irp.f`

    .. code:: fortran

        logical	:: libcint_nuc_available	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`libcint_available`
       * :c:data:`mpi_master`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_num`
       * :c:data:`use_libcint`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_n_e`

 
.. c:var:: double libcint_scale


    File : :file:`libcint_ao_integrals/libcint_data.irp.f`

    .. code:: fortran

        double precision, allocatable	:: libcint_scale	(ao_num)


    Ratio of the |AO| normalization factor to the one of the first |AO| of
    the shell. Equal to 1 for `s` and `p` shells.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized`
       * :c:data:`ao_coef_normalized`
       * :c:data:`ao_first_of_shell`
       * :c:data:`ao_num`
       * :c:data:`ao_shell`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`libcint_available`

 
.. c:var:: bool use_libcint_effective


    File : :file:`libcint_ao_integrals/libcint_data.irp.f`

    .. code:: fortran

        logical	:: use_libcint_effective	


    :c:data:`use_libcint` restricted to the cases libcint can handle.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`libcint_available`
       * :c:data:`use_libcint`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_kinetic_integrals`
       * :c:data:`ao_overlap`
       * :c:data:`ao_two_e_integral_erf_schwartz`
       * :c:data:`ao_two_e_integral_schwartz`
       * :c:data:`ao_two_e_integrals_erf_in_map`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`cholesky_ao_num`

 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: void libcint_ao_kinetic(void)


    File : :file:`libcint_ao_integrals/libcint_one_e.irp.f`

    .. code:: fortran

        subroutine libcint_ao_kinetic(mat)


    Kinetic energy matrix in the |AO| basis, from libcint.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_kinetic_integrals`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`libcint_ao_one_e`

 
.. c:function:: void libcint_ao_nuc(void)


    File : :file:`libcint_ao_integrals/libcint_one_e.irp.f`

    .. code:: fortran

        subroutine libcint_ao_nuc(mat)


    Nucleus-electron attraction matrix in the |AO| basis, from libcint.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_n_e`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`libcint_ao_one_e`

 
.. c:function:: void libcint_ao_one_e(void)


    File : :file:`libcint_ao_integrals/libcint_one_e.irp.f`

    .. code:: fortran

        subroutine libcint_ao_one_e(op, mat)


    One-electron |AO| integrals with libcint: overlap (`LIB_OP_OVLP`),
    kinetic energy (`LIB_OP_KIN`) or nucleus-electron attraction
    (`LIB_OP_NUC`).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_first_of_shell`
       * :c:data:`ao_l`
       * :c:data:`ao_l`
       * :c:data:`ao_num`
       * :c:data:`libcint_atm`
       * :c:data:`libcint_atm`
       * :c:data:`libcint_atm`
       * :c:data:`libcint_natm`
       * :c:data:`libcint_nbas`
       * :c:data:`libcint_scale`
       * :c:data:`shell_ang_mom`
       * :c:data:`shell_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`libcint_ao_kinetic`
       * :c:func:`libcint_ao_nuc`
       * :c:func:`libcint_ao_overlap`

 
.. c:function:: void libcint_ao_overlap(void)


    File : :file:`libcint_ao_integrals/libcint_one_e.irp.f`

    .. code:: fortran

        subroutine libcint_ao_overlap(mat)


    Overlap matrix in the |AO| basis, from libcint.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_overlap`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`libcint_ao_one_e`

 
.. c:function:: void libcint_ao_two_e_integral(void)


    File : :file:`libcint_ao_integrals/libcint_two_e.irp.f`

    .. code:: fortran

        double precision function libcint_ao_two_e_integral(i, j, k, l, lr)


    :math:`\langle ij | kl \rangle` with libcint, plain or long-range
    Coulomb operator.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_first_of_shell`
       * :c:data:`ao_shell`
       * :c:data:`libcint_atm`
       * :c:data:`libcint_atm`
       * :c:data:`libcint_atm`
       * :c:data:`libcint_env_lr`
       * :c:data:`libcint_natm`
       * :c:data:`libcint_nbas`
       * :c:data:`libcint_scale`
       * :c:data:`mu_erf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`libcint_two_e_block`

