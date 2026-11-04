.. _module_symmetry_mo: 
 
.. program:: symmetry_mo 
 
.. default-role:: option 
 
===========
symmetry_mo
===========

Molecular-orbital irreps. This module turns the point group found by
:ref:`module_symmetry` into labels for the molecular orbitals, writes them
to the EZFIO at the end of SCF/DFT, and provides the labels used by the
determinant and state layer (:ref:`module_symmetry_det`).

Needed Modules
==============
* :ref:`module_symmetry` -- point group, character tables
* :ref:`module_mo_basis`, :ref:`module_ao_basis`,
  :ref:`module_ao_one_e_ints` -- the MO coefficients and the AO overlap
* :ref:`module_utils`

What it adds
============

``pg_xg`` and the AO transformation
-----------------------------------

For a group element :math:`g`, the module builds its matrix
:math:`M_g` in the AO basis. Contracted shells are permuted by the atom
permutation of the element (``pg_shell_target``); inside a shell, the
block is the symmetric power :math:`\mathrm{Sym}^l(R_g)` of the rotation,
rescaled by the ratio of the radial factors of the primitives. Nothing is
hardcoded about the angular convention: the monomial representation is
expanded from :math:`R_g` directly.

The AO basis must support the group: it must be invariant under it, and the
coefficients of a shell must factorize as :math:`A_j B_p` (a radial factor
times an angular factor). ``pg_basis_ierr`` is 0 when this holds and 1, 2 or
3 otherwise, with ``pg_basis_msg`` naming the atom, shell, angular momentum
and element. The diagnostic is reported, never thrown from a provider.

The matrix of the element in the MO basis is

.. math::

   X_g = C^T S M_g C ,

where :math:`C` is the MO coefficient matrix and :math:`S` the AO overlap.
It is built once per element by ``pg_get_mo_rep``; the determinant layer
keeps the set in ``pg_xg``.

Symmetry blocks and MO labels
-----------------------------

A *symmetry block* is a connected component of the graph on the MOs with an
edge :math:`(i,j)` when :math:`|X_g(j,i)|^2 >` :option:`symmetry pg_thresh` for some
element :math:`g`. The test is on a weight, the same number and the same
units as the label cut, because the assembly defect is a weight defect: the
spurious amplitude between two MOs of different irreps is the square root of
the defect. One block is one irrep space; two copies of the same irrep at
different energies stay separate blocks.

The projector on irrep :math:`a` is the same formula everywhere in the
module:

.. math::

   P^a = \frac{d_a}{m_a |G|} \sum_c \chi_a(c) U_c ,

with :math:`d_a` the dimension (``pg_irrep_dim``), :math:`m_a` 2 when the
table row is the sum of a complex-conjugate pair and 1 otherwise
(``pg_irrep_mult``), and :math:`\chi_a(c)` the characters of the file column
:math:`c`. The weight of irrep :math:`a` in MO :math:`i` is
``pg_mo_weight(i,a) = p_a(i) = <i|P^a|i>``. Then

* ``pg_mo_irrep(i)`` is the argmax of :math:`p_a(i)`;
* it is 0 when :math:`\max_a p_a(i) < 1 -` :option:`symmetry pg_thresh`: a mixed MO is
  not given a best guess;
* ``pg_mo_purity(i)`` is :math:`\max_a p_a(i)`, equal to 1 for a symmetry
  adapted MO and also 1 for every partner of a degenerate block, where
  :math:`P^a` is the identity;
* ``pg_irrep_weight(a)`` is the row sum of the weights, i.e. the number of
  MOs carrying irrep :math:`a`. It is never recomputed from a second trace
  formula: the gates test the very numbers the report prints.

The report at the end of SCF/DFT
--------------------------------

``pg_report_mo_symmetry`` is called at the end of :ref:`module_hartree_fock`
and :ref:`module_kohn_sham`. It prints the point group and the number of MOs
per irrep, and writes ``mo_symmetry`` when every MO has a definite irrep.
It never stops and prints nothing when :option:`symmetry auto_symmetry` is false.

``pg_write_mo_symmetry`` writes ``mo_basis/mo_symmetry``: the integer irrep
index of each MO, with the convention "same integer = same irrep". A mixed
MO has index 0 and nothing is written in that case, because writing 0 would
declare all mixed MOs to belong to one irrep. The in-memory provider is
invalidated with ``soft_touch`` after the write, so the current process and
the EZFIO agree.

Guards
------

``pg_symmetry_ready()`` is true when the group was detected, its table is
usable and the AO basis supports it. ``pg_ao_ok`` is the same statement for
the basis alone. ``pg_mo_labels_ok()`` is true when every MO has a definite
irrep. These functions never stop: they are what keeps a failed detection
from aborting a CI run. A provider of this module is only entered when its
guard has already passed (IRPF90 hoists a provider build to the entry of the
procedure that names it, so the guard is always written in the caller, never
in the worker).

Report program
==============

``qp_run symmetry_mo DIR`` (source: ``symmetry_mo.irp.f``, ``program
symmetry``) prints:

* the canonical frame and the measured geometry residual;
* the group, its order, its classes and the character table;
* the AO-basis gates and ``delta_MO``;
* every MO label and the gates G1-G6.

With ``PG_WRITE_MO_SYMMETRY=1`` in the environment it also writes
``mo_basis/mo_symmetry``. Nothing is written when a gate failed or when an
MO has no definite irrep, so a written field is always a checked partition.

Gates
=====

Every gate prints its measured value and its bound; a failing gate is
reported and the run continues, so one run shows every number, and the
program exits 1 at the very end.

.. list-table::
   :header-rows: 1

   * - Gate
     - Test
     - Bound
   * - G-iso
     - :math:`M_g^T S M_g = S` for every element, with the raw AO overlap
     - ``1e-7``
   * - G1
     - :math:`C^T S C = 1` (orthonormality of the MO set)
     - ``1e-7``
   * - G2
     - :math:`X_g X_h = X_{gh}`, products taken from the permutation table
     - ``1e-7``
   * - G3
     - every :math:`\mathrm{tr}\,P^a` is an integer multiple of :math:`d_a`
     - ``1e-6``
   * - G4
     - :math:`\sum_a \mathrm{tr}\,P^a =` ``mo_num``
     - ``1e-6``
   * - G5
     - :math:`0 \le p_a(i) \le 1`, :math:`\sum_a p_a(i) = 1`,
       :math:`\max_a p_a(i) = 1`, with the window ``max(1e-8, 10*delta_MO)``
     - window
   * - G6
     - :math:`P^a P^a v = P^a v` and :math:`\sum_a P^a v = v` on a fixed
       pseudo-random vector
     - ``1e-6``

``delta_MO = max_g |X_g^T X_g - 1|`` is printed but not gated: it is G-iso
transported into the MO metric, and it is the accuracy floor of every purity
of this module. G3 is the authoritative test of the column assignment of the
character table: a swap of ``sv`` and ``sd`` in an Abelian group keeps every
orthogonality sum intact but changes the number of MOs per irrep.

The ``1e-6`` bounds are not machine precision and must not be lowered:
the representation matrices are built from *representatives* of the group
elements, which are only as accurate as the geometry tolerance. Measured on
NH3, ``X_g X_gh`` is off by ``2.3e-8`` at ``pg_tol = 1e-4`` bohr, and
:math:`P^a` is a linear combination of the :math:`X_g`, so
:math:`||P^2 - P||` grows with the same floor. A wrong column assignment, a
wrong :math:`m_a` or a missing irrep is O(1) in these gates.

Mixed MOs
=========

A molecular orbital that is not symmetry adapted gets index 0 and purity
below 1. This is not a bug in the labels: it means the SCF converged to a
solution that breaks the point group, that the geometry is not symmetric
within :option:`symmetry pg_tol`, or that the EZFIO database holds no
``mo_basis/mo_coef`` and QP2 falls back to the orthonormalized atomic
orbitals. The report names the mixed MOs with the two largest weights of
each, the number of blocks and ``delta_MO``, so that the cause can be
identified.

The determinant and state layers keep working on such a basis. ``X_g`` is
block diagonal on the blocks of this module by construction, so the
per-block occupations of a determinant never change under the group: the
family projector stays exact, and the irrep of a determinant or of a state
is an exact probability for the stored determinant space.
:ref:`module_symmetry_det` prints those probabilities. Only the label is
lost, never the symmetry information.

What a mixed basis does disable is everything that needs a label instead of
a distribution: state following in Davidson and the :option:`irrep_states`
filter are off, because they call ``pg_labels_usable``, which tests
``pg_mo_labels_ok``; and ``mo_basis/mo_symmetry`` is not written. No program
stops.

``qp_run symmetry_mo`` is a diagnostic: on such a basis it prints the whole
report and then exits 1, because gate G5 is precisely "the MO set is
symmetry adapted". That exit code is the verdict on the MO set, and is
expected on a fixture that is meant to exercise the mixed path. The CI
programs and ``qp_run symmetry_det`` are unaffected; ``symmetry_det`` exits
on its own gates.

The cost of the mixed path is one diagonal projector per printed
determinant, and the usual family projectors for the states; both are
estimated in :ref:`module_symmetry_det`.
 
 
 
Programs 
-------- 
 
 * :ref:`symmetry_mo` 
 
Providers 
--------- 
 
.. c:var:: bool pg_ao_ok


    File : :file:`auto_report.irp.f`

    .. code:: fortran

        logical	:: pg_ao_ok	


    True when the AO basis is invariant under the detected point group
    and every shell's coefficients factorize as A_j*B_p, i.e. when the
    MO representation matrices of ao_transform.irp.f are exact
    (:c:data:`pg_basis_ierr` = 0).
    
    :c:data:`pg_basis_ierr` is the full test: :c:func:`pg_check_ao_factor` alone would
    miss the shell-to-shell matching failures (ierr 1 and 2), which only
    :c:func:`pg_ao_rep_ex` would turn into a zero matrix, later and with a second
    message.
    
    Being a provider, the diagnostic is printed once per process: an early
    return still sets the built flag, so the message is cached with the
    value and never re-printed until a TOUCH.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`pg_shell_target`


 
.. c:var:: int pg_basis_ierr


    File : :file:`ao_transform.irp.f`

    .. code:: fortran

        integer, allocatable	:: pg_shell_target	(shell_num,pg_n_elem)
        integer	:: pg_basis_ierr	
        character*(512)	:: pg_basis_msg	


    ``pg_shell_target(s,g)`` is the unique shell on the atom
    ``pg_perm(basis_nucleus_index(s),g)`` whose contraction equals the one of
    shell ``s``. :c:data:`pg_basis_ierr` is 0 when the AO basis is invariant under the
    point group and the AO coefficients factorize as A_j*B_p inside every
    shell; otherwise it is 1 (no unique image shell), 2 (ambiguous image
    shell) or 3 (the primitive factor depends on the AO powers) and
    :c:data:`pg_basis_msg` says which atom, shell, angular momentum and element.
    
    Reporting is left to the caller (:c:data:`pg_ao_ok` and the gate report): a
    provider that stops cannot be reported on.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized`
       * :c:data:`ao_first_of_shell`
       * :c:data:`ao_prim_num`
       * :c:data:`ao_prim_num_max`
       * :c:data:`basis_nucleus_index`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_rot_mat`
       * :c:data:`shell_ang_mom`
       * :c:data:`shell_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_ao_ok`

 
.. c:var:: char pg_basis_msg


    File : :file:`ao_transform.irp.f`

    .. code:: fortran

        integer, allocatable	:: pg_shell_target	(shell_num,pg_n_elem)
        integer	:: pg_basis_ierr	
        character*(512)	:: pg_basis_msg	


    ``pg_shell_target(s,g)`` is the unique shell on the atom
    ``pg_perm(basis_nucleus_index(s),g)`` whose contraction equals the one of
    shell ``s``. :c:data:`pg_basis_ierr` is 0 when the AO basis is invariant under the
    point group and the AO coefficients factorize as A_j*B_p inside every
    shell; otherwise it is 1 (no unique image shell), 2 (ambiguous image
    shell) or 3 (the primitive factor depends on the AO powers) and
    :c:data:`pg_basis_msg` says which atom, shell, angular momentum and element.
    
    Reporting is left to the caller (:c:data:`pg_ao_ok` and the gate report): a
    provider that stops cannot be reported on.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized`
       * :c:data:`ao_first_of_shell`
       * :c:data:`ao_prim_num`
       * :c:data:`ao_prim_num_max`
       * :c:data:`basis_nucleus_index`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_rot_mat`
       * :c:data:`shell_ang_mom`
       * :c:data:`shell_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_ao_ok`

 
.. c:var:: double pg_col_diag


    File : :file:`mo_irrep.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pg_col_diag	(mo_num,pg_n_col)
        integer, allocatable	:: pg_mo_block	(mo_num)
        integer	:: pg_n_mo_block	
        integer, allocatable	:: pg_mo_block_dim	(mo_num)


    ``pg_col_diag(i,c)`` = sum over the ELEMENTS of column ``c`` of X_g(i,i), and
    the symmetry blocks of the MO set.
    
    A block is a connected component of the graph on the MOs with an edge
    (i,j), i /= j, when :math:`|X_g(j,i)|^2` > pg_thresh for SOME element g. The edge
    test is on a WEIGHT, the same number and the same units as the label cut:
    the assembly defect delta_MO is a weight defect, so the spurious amplitude
    between two MOs of different irreps is sqrt(delta_MO). A
    block is an irrep space: two different copies of the same irrep stay
    separate because X_g is block diagonal in the symmetry-adapted basis.
    
    ``pg_mo_block(i)`` is the component id, renumbered 1..:c:data:`pg_n_mo_block` in
    order of first appearance so that the report is deterministic;
    ``pg_mo_block_dim(i)`` is the size of the block of i.
    
    Nothing here is stored across the group: one X_g at a time, one pass.
    The guard on mo_num == 0 keeps a nuclei-only EZFIO (no ao_basis, hence no
    MOs) from reaching any of this.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`pg_n_col`
       * :c:data:`pg_thresh`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_mo_weight`

 
.. c:var:: int pg_irrep_weight


    File : :file:`mo_irrep.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pg_mo_weight	(mo_num,pg_n_irrep)
        integer, allocatable	:: pg_mo_irrep	(mo_num)
        double precision, allocatable	:: pg_mo_purity	(mo_num)
        integer, allocatable	:: pg_irrep_weight	(pg_n_irrep)


    ``pg_mo_weight(i,a)`` = p_a(i) = (d_a/(mult_a*:math:`|G|`)) * sum_c chi_a(c)
    pg_col_diag(i,c), the expectation value of the projector P^a on MO i.
    
    ``pg_mo_irrep(i)`` = argmax_a p_a(i), or 0 when max_a p_a(i) <
    1 - pg_thresh: a mixed MO must not get a best-guess label.
    ``pg_mo_purity(i)`` = max_a p_a(i), 1 for a symmetry-adapted MO and 1 for
    every partner of a degenerate block (P^a is the identity there).
    
    ``pg_irrep_weight(a)`` = nint(sum_i p_a(i)) = tr P^a, the number of times
    irrep a appears. It is the SUM OF THE ROWS OF :c:data:`pg_mo_weight`, never a
    second trace formula built from the characters: two mathematically equal
    expressions would only ever disagree through a bug in one of them, and
    gates G3/G4 must test the very numbers the report prints.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_thresh`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_n_mo_mixed`

 
.. c:var:: int pg_mo_block


    File : :file:`mo_irrep.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pg_col_diag	(mo_num,pg_n_col)
        integer, allocatable	:: pg_mo_block	(mo_num)
        integer	:: pg_n_mo_block	
        integer, allocatable	:: pg_mo_block_dim	(mo_num)


    ``pg_col_diag(i,c)`` = sum over the ELEMENTS of column ``c`` of X_g(i,i), and
    the symmetry blocks of the MO set.
    
    A block is a connected component of the graph on the MOs with an edge
    (i,j), i /= j, when :math:`|X_g(j,i)|^2` > pg_thresh for SOME element g. The edge
    test is on a WEIGHT, the same number and the same units as the label cut:
    the assembly defect delta_MO is a weight defect, so the spurious amplitude
    between two MOs of different irreps is sqrt(delta_MO). A
    block is an irrep space: two different copies of the same irrep stay
    separate because X_g is block diagonal in the symmetry-adapted basis.
    
    ``pg_mo_block(i)`` is the component id, renumbered 1..:c:data:`pg_n_mo_block` in
    order of first appearance so that the report is deterministic;
    ``pg_mo_block_dim(i)`` is the size of the block of i.
    
    Nothing here is stored across the group: one X_g at a time, one pass.
    The guard on mo_num == 0 keeps a nuclei-only EZFIO (no ao_basis, hence no
    MOs) from reaching any of this.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`pg_n_col`
       * :c:data:`pg_thresh`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_mo_weight`

 
.. c:var:: int pg_mo_block_dim


    File : :file:`mo_irrep.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pg_col_diag	(mo_num,pg_n_col)
        integer, allocatable	:: pg_mo_block	(mo_num)
        integer	:: pg_n_mo_block	
        integer, allocatable	:: pg_mo_block_dim	(mo_num)


    ``pg_col_diag(i,c)`` = sum over the ELEMENTS of column ``c`` of X_g(i,i), and
    the symmetry blocks of the MO set.
    
    A block is a connected component of the graph on the MOs with an edge
    (i,j), i /= j, when :math:`|X_g(j,i)|^2` > pg_thresh for SOME element g. The edge
    test is on a WEIGHT, the same number and the same units as the label cut:
    the assembly defect delta_MO is a weight defect, so the spurious amplitude
    between two MOs of different irreps is sqrt(delta_MO). A
    block is an irrep space: two different copies of the same irrep stay
    separate because X_g is block diagonal in the symmetry-adapted basis.
    
    ``pg_mo_block(i)`` is the component id, renumbered 1..:c:data:`pg_n_mo_block` in
    order of first appearance so that the report is deterministic;
    ``pg_mo_block_dim(i)`` is the size of the block of i.
    
    Nothing here is stored across the group: one X_g at a time, one pass.
    The guard on mo_num == 0 keeps a nuclei-only EZFIO (no ao_basis, hence no
    MOs) from reaching any of this.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`pg_n_col`
       * :c:data:`pg_thresh`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_mo_weight`

 
.. c:var:: int pg_mo_irrep


    File : :file:`mo_irrep.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pg_mo_weight	(mo_num,pg_n_irrep)
        integer, allocatable	:: pg_mo_irrep	(mo_num)
        double precision, allocatable	:: pg_mo_purity	(mo_num)
        integer, allocatable	:: pg_irrep_weight	(pg_n_irrep)


    ``pg_mo_weight(i,a)`` = p_a(i) = (d_a/(mult_a*:math:`|G|`)) * sum_c chi_a(c)
    pg_col_diag(i,c), the expectation value of the projector P^a on MO i.
    
    ``pg_mo_irrep(i)`` = argmax_a p_a(i), or 0 when max_a p_a(i) <
    1 - pg_thresh: a mixed MO must not get a best-guess label.
    ``pg_mo_purity(i)`` = max_a p_a(i), 1 for a symmetry-adapted MO and 1 for
    every partner of a degenerate block (P^a is the identity there).
    
    ``pg_irrep_weight(a)`` = nint(sum_i p_a(i)) = tr P^a, the number of times
    irrep a appears. It is the SUM OF THE ROWS OF :c:data:`pg_mo_weight`, never a
    second trace formula built from the characters: two mathematically equal
    expressions would only ever disagree through a bug in one of them, and
    gates G3/G4 must test the very numbers the report prints.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_thresh`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_n_mo_mixed`

 
.. c:var:: double pg_mo_purity


    File : :file:`mo_irrep.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pg_mo_weight	(mo_num,pg_n_irrep)
        integer, allocatable	:: pg_mo_irrep	(mo_num)
        double precision, allocatable	:: pg_mo_purity	(mo_num)
        integer, allocatable	:: pg_irrep_weight	(pg_n_irrep)


    ``pg_mo_weight(i,a)`` = p_a(i) = (d_a/(mult_a*:math:`|G|`)) * sum_c chi_a(c)
    pg_col_diag(i,c), the expectation value of the projector P^a on MO i.
    
    ``pg_mo_irrep(i)`` = argmax_a p_a(i), or 0 when max_a p_a(i) <
    1 - pg_thresh: a mixed MO must not get a best-guess label.
    ``pg_mo_purity(i)`` = max_a p_a(i), 1 for a symmetry-adapted MO and 1 for
    every partner of a degenerate block (P^a is the identity there).
    
    ``pg_irrep_weight(a)`` = nint(sum_i p_a(i)) = tr P^a, the number of times
    irrep a appears. It is the SUM OF THE ROWS OF :c:data:`pg_mo_weight`, never a
    second trace formula built from the characters: two mathematically equal
    expressions would only ever disagree through a bug in one of them, and
    gates G3/G4 must test the very numbers the report prints.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_thresh`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_n_mo_mixed`

 
.. c:var:: double pg_mo_weight


    File : :file:`mo_irrep.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pg_mo_weight	(mo_num,pg_n_irrep)
        integer, allocatable	:: pg_mo_irrep	(mo_num)
        double precision, allocatable	:: pg_mo_purity	(mo_num)
        integer, allocatable	:: pg_irrep_weight	(pg_n_irrep)


    ``pg_mo_weight(i,a)`` = p_a(i) = (d_a/(mult_a*:math:`|G|`)) * sum_c chi_a(c)
    pg_col_diag(i,c), the expectation value of the projector P^a on MO i.
    
    ``pg_mo_irrep(i)`` = argmax_a p_a(i), or 0 when max_a p_a(i) <
    1 - pg_thresh: a mixed MO must not get a best-guess label.
    ``pg_mo_purity(i)`` = max_a p_a(i), 1 for a symmetry-adapted MO and 1 for
    every partner of a degenerate block (P^a is the identity there).
    
    ``pg_irrep_weight(a)`` = nint(sum_i p_a(i)) = tr P^a, the number of times
    irrep a appears. It is the SUM OF THE ROWS OF :c:data:`pg_mo_weight`, never a
    second trace formula built from the characters: two mathematically equal
    expressions would only ever disagree through a bug in one of them, and
    gates G3/G4 must test the very numbers the report prints.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_thresh`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_n_mo_mixed`

 
.. c:var:: int pg_n_mo_block


    File : :file:`mo_irrep.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pg_col_diag	(mo_num,pg_n_col)
        integer, allocatable	:: pg_mo_block	(mo_num)
        integer	:: pg_n_mo_block	
        integer, allocatable	:: pg_mo_block_dim	(mo_num)


    ``pg_col_diag(i,c)`` = sum over the ELEMENTS of column ``c`` of X_g(i,i), and
    the symmetry blocks of the MO set.
    
    A block is a connected component of the graph on the MOs with an edge
    (i,j), i /= j, when :math:`|X_g(j,i)|^2` > pg_thresh for SOME element g. The edge
    test is on a WEIGHT, the same number and the same units as the label cut:
    the assembly defect delta_MO is a weight defect, so the spurious amplitude
    between two MOs of different irreps is sqrt(delta_MO). A
    block is an irrep space: two different copies of the same irrep stay
    separate because X_g is block diagonal in the symmetry-adapted basis.
    
    ``pg_mo_block(i)`` is the component id, renumbered 1..:c:data:`pg_n_mo_block` in
    order of first appearance so that the report is deterministic;
    ``pg_mo_block_dim(i)`` is the size of the block of i.
    
    Nothing here is stored across the group: one X_g at a time, one pass.
    The guard on mo_num == 0 keeps a nuclei-only EZFIO (no ao_basis, hence no
    MOs) from reaching any of this.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`pg_n_col`
       * :c:data:`pg_thresh`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_mo_weight`

 
.. c:var:: int pg_n_mo_mixed


    File : :file:`mo_irrep.irp.f`

    .. code:: fortran

        integer	:: pg_n_mo_mixed	


    Number of MOs with :c:data:`pg_mo_irrep` = 0, i.e. ``max_a p_a < 1 - pg_thresh``.
    This is the single definition every report and the warning of the
    determinant labelling use; it used to be recomputed in five places.
    
    It counts MOs, never blocks: the partners of a degenerate block of a
    symmetry-adapted set all carry a label, so they contribute 0 here. A
    nonzero count means that the basis is not symmetry adapted, and then
    the determinant and state irreps are probabilities, not labels.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`pg_mo_weight`


 
.. c:var:: int pg_shell_target


    File : :file:`ao_transform.irp.f`

    .. code:: fortran

        integer, allocatable	:: pg_shell_target	(shell_num,pg_n_elem)
        integer	:: pg_basis_ierr	
        character*(512)	:: pg_basis_msg	


    ``pg_shell_target(s,g)`` is the unique shell on the atom
    ``pg_perm(basis_nucleus_index(s),g)`` whose contraction equals the one of
    shell ``s``. :c:data:`pg_basis_ierr` is 0 when the AO basis is invariant under the
    point group and the AO coefficients factorize as A_j*B_p inside every
    shell; otherwise it is 1 (no unique image shell), 2 (ambiguous image
    shell) or 3 (the primitive factor depends on the AO powers) and
    :c:data:`pg_basis_msg` says which atom, shell, angular momentum and element.
    
    Reporting is left to the caller (:c:data:`pg_ao_ok` and the gate report): a
    provider that stops cannot be reported on.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized`
       * :c:data:`ao_first_of_shell`
       * :c:data:`ao_prim_num`
       * :c:data:`ao_prim_num_max`
       * :c:data:`basis_nucleus_index`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_rot_mat`
       * :c:data:`shell_ang_mom`
       * :c:data:`shell_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_ao_ok`

 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: void pg_ao_ratio(void)


    File : :file:`ao_transform.irp.f`

    .. code:: fortran

        double precision function pg_ao_ratio(p, q)


    ``B_p/B_q`` from ``ao_coef_normalized(p,j) = A_j B_p``: the ``j`` dependence
    cancels, so any primitive gives the same number. ``jref`` is the primitive
    carrying the largest ``ao_coef_normalized(q,j)`` so that the division is by
    the largest coefficient available (:c:func:`pg_check_ao_factor` has already
    proved that the ratio does not depend on ``j``).
    
    This is exact in both EZFIO conventions. In the pre-2025 one
    :c:data:`ao_coef_normalization_factor` is NOT B_p: the primitive factor there
    carries a power-dependent part, so using it would make a d shell not
    rotation covariant by an O(1) factor (found by G-iso on ammonia/6-31G*).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized`
       * :c:data:`ao_prim_num`

 
.. c:function:: void pg_ao_rep(void)


    File : :file:`ao_transform.irp.f`

    .. code:: fortran

        subroutine pg_ao_rep(ig, m)


    ``m`` = M_g in the AO basis, a dense (ao_num x ao_num) matrix, block
    diagonal by shell: shell ``s`` maps onto ``pg_shell_target(s,g)`` and its
    block is Sym^l(R_g) rescaled by B_p/B_q.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_get_mo_rep`
       * :c:func:`pg_report_ao`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pg_ao_rep_ex`

 
.. c:function:: void pg_ao_rep_ex(void)


    File : :file:`ao_transform.irp.f`

    .. code:: fortran

        subroutine pg_ao_rep_ex(ig, m)


    Block-diagonal builder: the block of shell ``s`` is Sym^l(R_g) rescaled by
    B_p/B_q and written at the rows of shell ``pg_shell_target(s,g)``.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_first_of_shell`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`pg_detection_failed`
       * :c:data:`pg_rot_mat`
       * :c:data:`pg_shell_target`
       * :c:data:`shell_ang_mom`
       * :c:data:`shell_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_ao_rep`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pg_monomial_rep`

 
.. c:function:: void pg_check_ao_factor(void)


    File : :file:`ao_transform.irp.f`

    .. code:: fortran

        subroutine pg_check_ao_factor(ierr, msg)


    Verifies, for every shell, that the ratios
    ``ao_coef_normalized(p,j)/ao_coef_normalized(p0,j)`` are independent of the
    primitive ``j`` (``p0`` = first AO of the shell). Cross-multiplied so that a
    zero coefficient needs no special case. This is the j-independence test
    of the representation: no double factorial, no hardcoded convention.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized`
       * :c:data:`ao_first_of_shell`
       * :c:data:`ao_prim_num`
       * :c:data:`basis_nucleus_index`
       * :c:data:`shell_ang_mom`
       * :c:data:`shell_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_shell_target`

 
.. c:function:: void pg_g2_pair(void)


    File : :file:`symmetry_mo.irp.f`

    .. code:: fortran

        subroutine pg_g2_pair(g, h, dev)


    ``dev`` = max ``X_g X_h - X_{gh}`` for one pair. The product element comes
    from :c:data:`pg_elem_mul`, which is built from the permutations: no matrix is
    ever matched.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`pg_elem_mul`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_report_ao`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pg_get_mo_rep`

 
.. c:function:: void pg_get_mo_rep(void)


    File : :file:`ao_transform.irp.f`

    .. code:: fortran

        subroutine pg_get_mo_rep(ig, x)


    ``x`` = C^T S M_g C, the matrix of group element ``ig`` in the MO basis.
    Nothing is kept between calls: the caller loops on
    the elements and holds one ``x`` at a time.
    
    ``matmul``, not ``dgemm``: the matrices are at most a few hundred wide in
    practice and correctness matters more than the constant. Three dgemms
    replace it verbatim if a production run ever makes this the hot spot.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_col_diag`
       * :c:func:`pg_g2_pair`
       * :c:func:`pg_report_ao`
       * :c:func:`pg_report_mo`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pg_ao_rep`

 
.. c:function:: void pg_mo_labels_ok(void)


    File : :file:`auto_report.irp.f`

    .. code:: fortran

        logical function pg_mo_labels_ok()


    True when every MO carries a definite irrep, i.e. when
    ``pg_mo_irrep(i)`` is positive for all i. Called only after
    :c:func:`pg_symmetry_ready` returned true, so the representation matrices are
    already built.
    
    This is the guard that keeps a labelling of a basis that is not
    symmetry adapted away from the Davidson fast path and from the
    ``mo_symmetry`` write. The reports themselves accept such a basis and
    print irrep probabilities instead.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`pg_mo_weight`

 
.. c:function:: void pg_mo_top2(void)


    File : :file:`mo_irrep.irp.f`

    .. code:: fortran

        subroutine pg_mo_top2(i, a1, a2)


    ``a1`` and ``a2`` are the indices of the two largest entries of
    :c:data:`pg_mo_weight` on MO ``i``; ``a2 = 0`` when the group has a single
    irrep. The reports use them to show the composition of a mixed MO
    (``A1 0.9856 B2 0.0144``) in one line; nothing else depends on them.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`pg_mo_weight`
       * :c:data:`pg_n_col`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_report_mo`

 
.. c:function:: void pg_monomial_rep(void)


    File : :file:`ao_transform.irp.f`

    .. code:: fortran

        subroutine pg_monomial_rep(l, r, np, powp, nq, powq, t)


    Builds ``t`` such that ``U_g v_p = sum_q t(q,p) v_q`` with
    ``v_p(r) = r_x^a r_y^b r_z^c`` and ``U_g f(r) = f(R^-1 r)``.
    ``R(j,i)`` appears with the *row* index running over x,y,z: expanding
    ``(sum_j R(j,i) r_j)^p_i``, the coefficient of ``r_j`` is ``R(j,i)``.
    
    The product is accumulated as p_i successive multiplications by the
    linear form ``sum_j R(j,i) r_j``, which needs no factorials and no
    hardcoded degree: it is the same expansion as the multinomial formula,
    only cheaper to get right.
    
    ``powp``, ``powq`` and ``t`` are declared with their *full* extent (nmax = 55
    monomials for lmax = 9), not with (3,np) and (nq,np): Fortran associates
    an explicit-shape dummy with the actual argument by sequence, so a dummy
    (3,np) would stride by np while the caller's array strides by nmax, and
    every column above the first would land at the wrong offset. The extents
    here must stay equal to those of :c:func:`pg_ao_rep`.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_ao_rep_ex`

 
.. c:function:: void pg_perm_string(void)


    File : :file:`symmetry_mo.irp.f`

    .. code:: fortran

        subroutine pg_perm_string(perm, nat, s)


    Comma-free list of ``perm(1:nat)`` in ``s``; truncated with "..." when ``s``
    is too short.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_report`

 
.. c:function:: void pg_ready_columns(void)


    File : :file:`auto_report.irp.f`

    .. code:: fortran

        subroutine pg_ready_columns(ok)


    ``ok`` is true when the AO basis supports the group and the character
    table of :c:data:`pg_group_name` could be parsed and validated. Never stops.
    
    :c:data:`pg_n_irrep` prints the missing-table diagnostic itself, and a table
    that fails a validation check sets :c:data:`pg_detection_failed` and leaves
    the count at 0.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_ao_ok`
       * :c:data:`pg_detection_failed`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_col`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_symmetry_ready`

 
.. c:function:: void pg_ready_detection(void)


    File : :file:`auto_report.irp.f`

    .. code:: fortran

        subroutine pg_ready_detection(ok)


    ``ok`` is true when a group name could be measured. Never stops.
    
    :c:data:`pg_inertia_eval` comes first: it rejects a geometry with no
    nucleus, or all nuclei at the centre, before any enumeration. The
    group name is a second step because a geometry that only looks
    symmetric gives up later, in the invariants, and an empty name must
    not be reported as ready.
    
    C1 deliberately passes: it is a legitimate group with a one-row
    table, so the report has one line per irrep and no extra branch.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_detection_failed`
       * :c:data:`pg_group_name`
       * :c:data:`pg_inertia_eval`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_symmetry_ready`

 
.. c:function:: void pg_report(void)


    File : :file:`symmetry_mo.irp.f`

    .. code:: fortran

        subroutine pg_report(write_mo)


    Prints the canonical frame, the group, its classes and the table.
    ``write_mo`` forwards the "-w" request of the driver to :c:func:`pg_report_mo`.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`
       * :c:data:`pg_center`
       * :c:data:`pg_class_size`
       * :c:data:`pg_geom_residual`
       * :c:data:`pg_group_name`
       * :c:data:`pg_has_i`
       * :c:data:`pg_is_planar`
       * :c:data:`pg_n_class`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_n_max`
       * :c:data:`pg_plane_normal`
       * :c:data:`pg_rot`
       * :c:data:`pg_rot_mat`
       * :c:data:`pg_thresh`
       * :c:data:`pg_tol`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`symmetry`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pg_character_table_path`
       * :c:func:`pg_decode_kind`
       * :c:func:`pg_perm_string`
       * :c:func:`pg_report_ao`
       * :c:func:`pg_report_mo`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`mo_symmetry`

 
.. c:function:: void pg_report_ao(void)


    File : :file:`symmetry_mo.irp.f`

    .. code:: fortran

        subroutine pg_report_ao(delta_mo)


    Prints the AO basis summary, then the gates:
    
    ``G-iso``
        M_g^T S M_g = S. Authoritative test of the AO transform, with S the
        RAW ao_overlap of the EZFIO, tolerance 1e-7; needs no MOs, so it
        runs on a nuclei+AO EZFIO.
    
    ``G1``
        C^T S C = 1. Orthonormality of the MO set in the AO metric; no group
        element enters it, and it runs whether or not the MO space is
        complete.
    
    ``G2``
        X_g X_h = X_gh. Only meaningful when mo_num == ao_num.
    
    delta_MO = max_g ``X_g^T X_g - 1`` is printed and NOT gated: it is G-iso
    transported into the MO metric (see the identity at the measurement).
    
    Each gate prints the measured deviation and stops on failure. The 1e-7
    bound of G-iso is not machine precision: :c:data:`ao_overlap` itself is only
    group symmetric to ~3e-9 (see the comment at the measurement below).
    
    ``delta_mo`` is returned so that the MO report can size its G5 window:
    it is the accuracy floor of every MO purity.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`ao_shell`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`pg_elem_mul`
       * :c:data:`pg_geom_residual`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_shell_target`
       * :c:data:`pg_thresh`
       * :c:data:`shell_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_report`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_has_ao_basis_ao_num`
       * :c:func:`pg_ao_rep`
       * :c:func:`pg_g2_pair`
       * :c:func:`pg_get_mo_rep`

 
.. c:function:: void pg_report_mo(void)


    File : :file:`symmetry_mo.irp.f`

    .. code:: fortran

        subroutine pg_report_mo(delta_mo, write_mo)


    Prints the label of every MO and the last four MO purity gates:
    
    ``G3``
        every tr P^a is an integer (1e-6) and a multiple of d_a, i.e. n_a
        copies of irrep a. Authoritative test of the column matching:
        a sv/sd swap in an Abelian group changes n_a and is invisible to
        every earlier test.
    
    ``G4``
        sum_a tr P^a = mo_num (1e-6).
    
    ``G5``
        0 <= p_a(i) <= 1 (w5), sum_a p_a(i) = 1 (1e-6) and
        max_a p_a(i) = 1 (w5), with w5 = max(1e-8, 10*delta_MO). The last
        one is the MO-set-is-symmetry-adapted test.
    
    ``G6``
        P^a P^a v = P^a v and sum_a P^a v = v (1e-6) for a fixed
        pseudo-random v: one pass over the group to build P^a v, then a
        second pass for P^a(P^a v). One X_g per element, never stored.
        For an irrep with w_a = 0 the ratio is pure representation noise,
        so the deviation of that projector is tested absolutely against
        the same 1e-6 (stricter than the ratio).
    
    Every gate prints its measured value; a failing gate prints FAILED and
    the subroutine keeps going, so that one run shows every number. It stops
    with status 1 only at the very end.
    
    The 1e-6 bound is not machine precision and must not be lowered: every
    X_g is built from the stored *representative* of an element, which is
    only as good as the geometry tolerance (measured: 2.3e-8 on X_g X_gh for
    the C3 axis of NH3 at pg_tol = 1e-4 bohr), and P^a is a linear
    combination of X_g, so :math:`||P^2 - P||` grows as (d_a/mult_a)^2 times that
    floor. A wrong column matching, a wrong mult_a or a
    missing irrep is O(1) in these gates.
    
    ``delta_mo`` is the representation floor measured by the AO report (the
    diagnostic max_g ``X_g^T X_g - 1``, implied by G-iso). It sizes the G5
    window, so that a well-conditioned but wide basis cannot fail G5 on
    integral noise alone.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_mo_weight`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_n_mo_mixed`
       * :c:data:`pg_shell_target`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_report`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_has_ao_basis_ao_num`
       * :c:func:`pg_get_mo_rep`
       * :c:func:`pg_mo_top2`
       * :c:func:`pg_write_mo_symmetry`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`mo_symmetry`

 
.. c:function:: void pg_report_mo_symmetry(void)


    File : :file:`auto_report.irp.f`

    .. code:: fortran

        subroutine pg_report_mo_symmetry()


    Called at the end of SCF/DFT: prints the point group and the number of
    MOs per irrep, and writes ``mo_basis/mo_symmetry`` when every MO has a
    definite irrep. Never stops, prints nothing when :c:data:`auto_symmetry` is
    false.
    
    This body names :c:func:`pg_symmetry_ready` only, so the guard runs before any
    provider of the chain is entered.

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pg_report_mo_symmetry_work`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`mo_symmetry`

 
.. c:function:: void pg_report_mo_symmetry_work(void)


    File : :file:`auto_report.irp.f`

    .. code:: fortran

        subroutine pg_report_mo_symmetry_work()


    Worker of :c:func:`pg_report_mo_symmetry`: names the chain, so it may only be
    called once :c:func:`pg_symmetry_ready` returned true.
    
    The count printed per irrep is :c:data:`pg_irrep_weight`, the row sum of the
    projector weights: it is the very number the gates test, and for a
    symmetry-adapted set it equals the number of MOs carrying the label.
    
    A mixed MO has :c:data:`pg_mo_irrep` = 0 and nothing is written: :c:data:`mo_symmetry`
    only carries "same integer = same irrep", so writing 0 would declare
    all the mixed MOs to be one single irrep.
    
    The mixed count runs on every rank; only the printing and the EZFIO
    write are under :c:data:`mpi_master`, as in save_wavefunction_general.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`mpi_master`
       * :c:data:`pg_geom_residual`
       * :c:data:`pg_group_name`
       * :c:data:`pg_mo_weight`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_n_mo_mixed`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_report_mo_symmetry`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pg_write_mo_symmetry`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`mo_symmetry`

 
.. c:function:: void pg_shell_same(void)


    File : :file:`ao_transform.irp.f`

    .. code:: fortran

        logical function pg_shell_same(s1, s2)


    True when the two shells have the same angular momentum and the same
    contraction: same number of primitives, same *sorted* exponents and the
    primitive coefficients that go with them, to 1e-10.
    ``ao_coef`` (raw), never :c:data:`ao_coef_normalized`: the latter legitimately
    differs between AOs of the same shell when their powers differ.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_first_of_shell`
       * :c:data:`ao_prim_num`
       * :c:data:`ao_prim_num_max`
       * :c:data:`shell_ang_mom`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pg_sorted_contraction`

 
.. c:function:: void pg_sorted_contraction(void)


    File : :file:`ao_transform.irp.f`

    .. code:: fortran

        subroutine pg_sorted_contraction(s, ex, cf, n)


    Exponents of shell ``s`` in increasing order, each with its own raw
    coefficient. The pair must be sorted together: a basis that lists its
    primitives in another order is still the same contraction.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef`
       * :c:data:`ao_expo`
       * :c:data:`ao_first_of_shell`
       * :c:data:`ao_prim_num`
       * :c:data:`ao_prim_num_max`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_shell_same`

 
.. c:function:: void pg_symmetry_ready(void)


    File : :file:`auto_report.irp.f`

    .. code:: fortran

        logical function pg_symmetry_ready()


    True when the point group was detected, its character table is
    usable and the AO basis supports the group: the MOs can be labelled
    and ``mo_basis/mo_symmetry`` written.
    
    Never stops. :c:data:`auto_symmetry` = F is silent; a failed detection has
    already printed its reason.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`auto_symmetry`
       * :c:data:`pg_detection_failed`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pg_ready_columns`
       * :c:func:`pg_ready_detection`

 
.. c:function:: void pg_write_mo_symmetry(void)


    File : :file:`mo_irrep.irp.f`

    .. code:: fortran

        subroutine pg_write_mo_symmetry()


    Writes :c:data:`pg_mo_irrep` into ``mo_basis/mo_symmetry`` of the EZFIO.
    
    :c:data:`mo_symmetry` is an existing QP2 provider (src/mo_basis/mo_class.irp.f)
    and IRPF90 forbids filling it from another module, so it cannot be a
    provider of this one. The user's program must call
    ``call pg_write_mo_symmetry()`` explicitly; no provider writes EZFIO.
    
    A mixed MO has ``pg_mo_irrep = 0``, and :c:data:`mo_symmetry` only carries "same
    integer = same irrep": writing 0 would silently declare all mixed MOs to
    be one irrep. Nothing is written in that case, and the count is printed.
    
    :c:data:`mo_symmetry` is usually built (all-1 default) before this point, so the
    EZFIO write alone would leave the in-memory provider stale: the
    ``soft_touch`` invalidates it without building it when it is not built.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_symmetry`
       * :c:data:`pg_mo_weight`
       * :c:data:`pg_n_mo_mixed`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_report_mo`
       * :c:func:`pg_report_mo_symmetry_work`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_set_mo_basis_mo_symmetry`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`mo_symmetry`

