========
Symmetry
========

Point-group detection and irrep labels for QP2: molecular point group, MO
irreps, determinant irreps, and state symmetry weights ``w_a = <Psi|P^a|Psi>``.

Needed Modules
==============
* ``devel/symmetry``     — ``nuclei``, ``mo_basis``, ``ao_one_e_ints``
* ``devel/symmetry_det`` — ``symmetry``, ``determinants``

Everything is specified, with the equations and the measured tests, in
``SPEC.md`` next to this file. ``NOTES-core-integration.md`` is the plan for
moving this into the core of QP2.

What the modules do
===================

``devel/symmetry``
------------------

* ``frame.irp.f`` — canonical frame (SPEC 4): centre of mass, principal axes
  recomputed with ``dsyev``, plane test and normal, ``pg_rot`` from the input
  frame to the canonical frame. Convention: ``z`` is the principal axis, ``x``
  is the plane normal for a planar molecule, and ``C2'`` is the ``C2`` class
  containing ``x`` — the Mulliken/Cotton labelling that literature uses.
* ``group.irp.f`` — the group is *enumerated*, not looked up: generators are
  taken from the geometry (atoms, midpoints, inertia axes, cross products),
  closed under composition, and identified among the 204 named finite point
  groups. An element is identified by its atom permutation together with the
  sign of its 3x3 matrix (``pg_perm``, ``pg_elem_det``), never by comparing
  floating-point matrices.
* ``chartab.irp.f`` — a pure parser for the 211 character tables in
  ``Symmetry/``. The group name comes from the file name, since the first line
  of some files is wrong.
* ``columns.irp.f`` — maps the group elements onto the *columns* of the table.
  A column may be the union of several geometric classes and their inverses
  (``T``: eight ``C3`` in two classes, one column), so every sum in the module
  runs over columns, with ``pg_elem_col`` giving the column of each element.
* ``ao_transform.irp.f`` — exact AO and MO representation matrices: for each
  element, the permutation of the contracted shells and the ``L`` matrix of
  the rotated spherical harmonics, then ``X_g = C^T S M_g C``.
* ``mo_irrep.irp.f`` — MO labels from the universal projector
  ``P^a = (d_a/(mult_a*|G|)) sum_c chi_a(c) U_c``. MOs are grouped into
  *blocks* (connected components of ``|X_g(j,i)|^2 > pg_thresh``, tested on the
  weight, not the amplitude); a block is one irrep space and its members share
  one label. Also ``pg_write_mo_symmetry``.
* ``symmetry_driver.irp.f`` — ``program symmetry``, binary
  ``symmetry_driver``: prints the frame, the group, the classes, the table,
  every MO label and gates G1-G6.

``devel/symmetry_det``
----------------------

* ``det_irrep.irp.f`` — determinant labels. ``<D'|U(g)|D> =
  prod_b det(X_g[A'_b,A_b]) prod_b det(X_g[B'_b,B_b])`` and the projector is
  the *same* ``P^a`` as for the MOs. Determinants with the same per-block spin
  counts form a *family*; ``T^a`` is a projector on the family and
  ``pg_det_purity = max_a T^a_DD``. A determinant is pure iff its family has
  one member.
* ``state_weights.irp.f`` — ``w_a(s) = sum_{D,D'} c_D(s) T^a_{DD'} c_D'(s)``,
  exact off-diagonal terms included, plus ``pg_state_irrep`` and the mixed
  determinant fraction.
* ``symmetry_det_driver.irp.f`` — ``symmetry_det_driver``: prints the family
  table, the determinant table (hidden above 200 determinants), the per-irrep
  weights of every state, and gates D1-D6.

Keywords
========

======================  ========  ==========================================
``pg_tol``              1e-4      Distance (bohr) under which two atoms of
                                  the same element are considered the same
                                  point.
``pg_thresh``           1e-6      Weight threshold: an off-diagonal
                                  ``|X_g(j,i)|^2`` above it connects two MOs
                                  into the same block, and a purity above
                                  ``1 - pg_thresh`` makes a label definite.
``pg_data_dir``         empty     Directory of the character tables. Empty or
                                  ``None`` means ``$QP_ROOT/data/Symmetry``.
======================  ========  ==========================================

The defaults live in ``settings.irp.f``, not in ``EZFIO.cfg``: under IRPF90 the
``default:`` line of a keyword is ignored and the provider stops when the key
is absent (SPEC 15). A virgin EZFIO therefore works.

Install and run
===============

The 211 character tables are installed outside the source tree::

    make -C ~/qp2/src/symmetry install     # -> $QP_ROOT/data/Symmetry
    make -C ~/qp2/src/symmetry uninstall

``QP_ROOT`` must be set; the target stops with a message when it is not.
Nothing in QP2 installs plugin data.

Activate the modules with symlinks (**absolute** paths: the relative form
resolves one level too high) and rebuild::

    ln -s ~/qp2/plugins/qp_plugins_scemama/devel/symmetry     ~/qp2/src/symmetry
    ln -s ~/qp2/plugins/qp_plugins_scemama/devel/symmetry_det ~/qp2/src/symmetry_det
    cd ~/qp2 && ./scripts/compilation/qp_create_ninja update
    cd ~/qp2 && bash -c 'source ./quantum_package.rc; \
        export CPATH=$QP_ROOT/include LIBRARY_PATH=$QP_ROOT/lib; ninja'

``qp_create_ninja update`` is required the first time (new ``.irp.f`` files);
after that plain ``ninja`` is enough. Sourcing ``quantum_package.rc`` is not
optional even for the Fortran build: the generator calls ``ei_handler.py``.

Then, outside the source tree, use the normal ``qp_run`` path::

    qp_create_ezfio -a -b cc-pvdz -o water.ezfio water.xyz   # -a: bohr
    qp_run scf water.ezfio
    qp_run symmetry_driver water.ezfio
    qp_run symmetry_det_driver water.ezfio

``qp_run`` runs ``qp_edit -c`` before the program, so the two drivers are only
trustworthy if checking the input does not modify it. That used to be false:
``qp_edit -c`` echoed every section back into EZFIO, and the coordinates went
through ``%16.8f`` Angstrom, i.e. a ``1.9e-8`` bohr grid. A geometry created in
bohr, where an exact symmetry relation holds to ``1e-16``, came back with
residuals of ``1e-8`` and the detection failed on a *good* input. The fix (see
SPEC 15) is four edits in QP2, not in this plugin:

* ``scripts/ezfio_interface/qp_edit_template``: the section writer takes
  ``?apply``; ``-c`` still parses each section (so a malformed input is still
  reported) but no longer writes it.
* ``ocaml/Point3d.ml``: Angstrom to bohr is a division by the unit, not a
  multiplication by its reciprocal, and coordinates print as ``%.17g``.
* ``ocaml/Charge.ml``: nuclear charges print as ``%.17g``, so ghost and
  fractional charges survive.
* ``ocaml/Input_nuclei_by_hand.ml``: the xyz block prints the charge with
  ``Charge.to_string`` instead of truncating it with ``Charge.to_int``.

What ``qp_edit -c`` still does is *materialize missing defaults*: every
``read_*`` that finds its keyword absent writes the default. On a fresh
directory that creates ``determinants/{bit_kind,n_int,mo_label,expected_s2,
psi_coef,psi_det,state_average_weight}`` and the ``ao_one_e_ints`` keywords; on
a second run it changes nothing. That behaviour is kept on purpose: it is what
``qp_reset --dets`` relies on to leave a usable directory.

To store the labels in the EZFIO (``mo_basis/mo_symmetry``, "same integer =
same irrep"), ask for it explicitly::

    PG_WRITE_MO_SYMMETRY=1 qp_run symmetry_driver water.ezfio

Nothing is written when a gate failed or when an MO has no definite irrep, so
the field is always a checked partition.

Reading the output
==================

* ``i  a  irrep  block  dim  purity`` — ``a`` is the irrep index (1 is the
  totally symmetric row by construction, 0 means "no definite irrep"), and the
  members of a degenerate block share one label: inside a block ``P^a`` is the
  identity, so no partner of a degenerate irrep has an individual label.
* ``sum_D c_D^2`` per irrep is the weight of the *determinants* labelled with
  that irrep; ``w_a`` from the state table is the weight of the *state* and is
  what an excited-state assignment should use. For a state built inside a
  degenerate shell the first is zero while the second is 1.
* A deficit ``sum_a w_a < sum_D c_D^2`` means ``psi_coef`` is not normalized in
  the determinant basis; it never means that the group is broken (SPEC 10.2).
* Gates: G1-G6 for the MO labels, D1-D6 for the determinants. Each gate prints
  its measured number and the bound; a gate failure exits 1 only after every
  number of the run has been printed.
* A closed-shell occupation of a partially filled degenerate shell has no
  one-determinant irrep. That is physics, not a bug: the run prints the family
  size, the purity, the split weights, and skips D6 for that wave function.

Tests
=====

``chartab_test`` parses the tables without any EZFIO::

    ~/qp2/src/symmetry/chartab_test ~/qp2/data/Symmetry/*     # PARSED 211/211

It is the one program called directly: it takes a list of table files, not an
EZFIO directory, so ``qp_run`` (which always passes the EZFIO path) cannot
drive it. It touches no geometry, so there is nothing to round-trip.
