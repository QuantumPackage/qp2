=====================
Symmetry determinant
=====================

Determinant and state symmetry: part two of ``devel/symmetry``.

The design, the equations, the conventions and the measured tests are in
``../symmetry/SPEC.md`` (sections 9 to 12); the keywords, the character-table
install and the general build recipe are in ``../symmetry/README.rst``. This
file only says what this module adds.

Needed Modules
==============
* ``symmetry``     -- point group, character tables, MO irreps
* ``determinants`` -- ``psi_det``, ``psi_coef``, the sorted-bit index

What it provides
================
``pg_det_irrep``, ``pg_det_purity``, ``pg_det_block``
    Irrep index, purity ``max_a T^a_DD`` and family size of every determinant.
    Irrep 1 is the totally symmetric row; 0 means no definite irrep. A
    determinant is pure if and only if its family -- the determinants with the
    same per-block spin counts -- has one member.

``pg_det_family``, ``pg_n_family``, ``pg_family_first``
    Family id (order of first appearance), number of families, and the
    ``psi_det`` index of each family representative, i.e. the member the
    projectors are built for.

``pg_state_irrep_weight``, ``pg_state_irrep``, ``pg_state_purity``
    ``w_a(s) = <Psi_s|P^a|Psi_s>`` with the off-diagonal ``T^a_{D'D}`` terms
    included, the dominant irrep (0 when none has weight ``> 1-pg_thresh``)
    and its weight. The same ``P^a`` as the MO labels, never a second formula.

``pg_det_mixed_frac``, ``pg_state_norm``
    Weight of state ``s`` carried by determinants with no definite irrep, and
    ``sum_D c_D(s)^2``, which equals ``sum_a w_a(s)``. A deficit between the
    two means ``psi_coef`` is not normalized; it never means that the CI space
    is not group invariant (SPEC 10.2).

``symmetry_det_driver``
    Prints the family table, the determinant table (truncated above 200
    determinants), the weights of every state, and the gates D1-D6. Each gate
    prints its measured value; the run exits 1 only after every number has
    been printed.

Physics note
============
An open degenerate shell has no one-determinant irrep, and neither has a
closed-shell occupation of a partially filled shell: with one electron in a
two-dimensional ``e`` block the determinant is a 50/50 mixture of ``A1`` and
``E``. That is the SCF solution breaking the point group, not a labelling
error. The determinant is reported with purity 0.5 and no label; the state
weights carry the symmetry instead, and the driver says so explicitly. On an
Abelian group every determinant is pure, whatever its open shells.

Build and run
=============
::

    ln -s ~/qp2/plugins/qp_plugins_scemama/devel/symmetry_det ~/qp2/src/symmetry_det
    cd ~/qp2 && ./scripts/compilation/qp_create_ninja update
    cd ~/qp2 && bash -c 'source ./quantum_package.rc; \
        export CPATH=$QP_ROOT/include LIBRARY_PATH=$QP_ROOT/lib; \
        ninja -f config/build.ninja module_symmetry_det'
    qp_run symmetry_det_driver water.ezfio

The symlink must be absolute: the relative form resolves one level too high.
``qp_create_ninja update`` is needed only when a new ``.irp.f`` file appears.

``qp_run`` is the normal path and runs ``qp_edit -c`` first. That is safe only
since the ``qp_edit`` fix described in ``../symmetry/README.rst`` (SPEC 15):
before it, checking the input rewrote every coordinate at ``%16.8f`` Angstrom.
The driver sets ``read_wf`` itself, so an EZFIO with no stored wave function is
read as the Hartree-Fock determinant (``N_det = 1``).

The MO labels must be clean first: the driver stops with the list of offending
MOs when a MO has no definite irrep (SPEC 9.5). Run
``src/symmetry/symmetry_driver`` on the same EZFIO and read its G1-G6 gates
before interpreting a determinant label.

Tests
=====
The two drivers are exercised by ``run_all.sh`` and ``run_all_det.sh`` in
``~/qp2_sym_tests`` (14 EZFIO fixtures each). Both scripts copy each fixture to
a scratch directory and drive it through ``qp_run``, so the fixtures are never
modified. One fixture is expected to exit 1 and is not a regression:
``water-6-31g_star`` has two mixed MOs, which is the SPEC 9.5 precondition.
``virgin`` (geometry, basis and electrons only) used to stop upstream on a
missing ``lin_dep_cutoff``; now that ``qp_edit -c`` materializes that default,
it runs to completion and labels the MOs that QP2 builds when ``mo_coef`` is
absent (``ao_ortho_canonical_coef``). All its gates pass (G1 ``6.6e-15``, G2
``6.8e-15``), so the driver also works before any SCF has been run.
The multi-determinant cases are ``detwork/water-cisd`` (679 determinants, all
``A1``) and ``detwork/methane-4st`` (5005 determinants, four states: ``A1``,
``T2``, ``T2``, ``T2``, largest family 81).
