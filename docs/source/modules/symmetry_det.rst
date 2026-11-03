.. _module_symmetry_det: 
 
.. program:: symmetry_det 
 
.. default-role:: option 
 
============
symmetry_det
============

Determinant and CI-state symmetry: the third layer of the symmetry support,
after the point group (:ref:`module_symmetry`) and the MO labels
(:ref:`module_symmetry_mo`). It labels the determinants of the wave function,
projects the CI states on the irreps, prints the state table, and implements
the :option:`symmetry irrep_states` filter and the irrep-aware state following used by
Davidson.

Needed Modules
==============
* :ref:`module_symmetry_mo` -- point group, character tables, MO irreps
* :ref:`module_determinants` -- ``psi_det``, ``psi_coef``, the sorted-bit
  index
* :ref:`module_utils`

The module deliberately does **not** need :ref:`module_davidson`. The
dependency goes the other way: ``davidson`` needs ``symmetry_det``. This is
what lets the CI programs use the labels without a cycle.

Determinants: families and projectors
=====================================

Two determinants belong to the same *family* when they have the same number
of spin-up and spin-down electrons in every symmetry block. A family is the
set of determinants that a group element can mix, so the projector is built
once per family, not once per determinant.

For an element :math:`g` and two members :math:`D, D'` of a family,

.. math::

   \langle D' | U(g) | D \rangle =
   \prod_b \det\left( X_g[A'_b, A_b] \right)
   \prod_b \det\left( X_g[B'_b, B_b] \right) ,

where :math:`A_b` and :math:`B_b` are the alpha and beta occupied MOs of the
determinant in block :math:`b`, and :math:`X_g` is the MO representation
matrix of :ref:`module_symmetry_mo`. Each block contributes a square
submatrix because the family has the same count in every block.

The projector on irrep :math:`a` is the same formula as everywhere else:

.. math::

   T^a = \frac{d_a}{m_a |G|} \sum_c \chi_a(c) U_c ,

restricted to the family. ``pg_det_irrep(D)`` is the argmax of
:math:`T^a_{DD}` over :math:`a`; it is 0 when
:math:`\max_a T^a_{DD} < 1 -` :option:`symmetry pg_thresh`. ``pg_det_purity(D)`` is
that maximum, 1 for a pure determinant and :math:`1/d_a` for one member of
an open degenerate shell.

A determinant is pure **if and only if its family has one member**: a family
of size greater than one exists only because some block holds between 1 and
:math:`d_b - 1` electrons, and every member then shares that partial
occupation. ``pg_det_block`` is the family size, ``pg_det_family`` the
family id, and ``pg_family_first`` the ``psi_det`` index of the
representative the module enumerates.

A truncated CI space is handled exactly: a family member that is not in
``psi_det`` has coefficient 0, so the missing part of the projector
contributes nothing. The label of a determinant never depends on which
other members happen to be selected.

CI states
=========

The weight of irrep :math:`a` in state :math:`s` is the exact diagonal
matrix element of the projector in the state,

.. math::

   w_a(s) = \sum_{D,D'} c_D(s)\, T^a_{DD'}\, c_{D'}(s) ,

off-diagonal family terms included. ``pg_state_irrep(s)`` is the argmax, or
0 when no weight reaches ``1 - pg_thresh``;
``pg_state_purity(s)`` is the maximum. ``pg_det_mixed_frac(s)`` is the
weight of the state carried by determinants with no definite irrep, and
``pg_state_norm(s) = sum_D c_D(s)^2``.

``sum_a w_a(s)`` equals ``pg_state_norm(s)`` for every state. A deficit is
only a bookkeeping statement about ``psi_coef`` (it is not normalized in the
determinant basis); it never means that the CI space is not group
invariant. The sum rule is reported, not gated: the exact mathematical
content is already in the projector construction.

The state table
---------------

``pg_print_state_irreps(energy, s2)`` is called by :ref:`cis`, :ref:`cisd`
and :ref:`fci` after diagonalization. It prints one line per state::

    symmetry: state irreps (N_det = 679)
    symmetry:    state      energy (Ha)       S^2  irrep   purity
    symmetry:         1    -76.1140581951    0.0000  A1   1.000

and nothing when :option:`symmetry state_symmetry` is false, when the detection
failed, or when the MO basis has mixed orbitals (one line says so). The
:math:`S^2` column is the caller's value, so a symmetry label and a spin
label are always read together.

The table uses the exact weights above, with the projectors for the families
of size greater than one. For a pure determinant the diagonal is a Kronecker
delta and the contribution is taken from the label, which is why the cost is
one pass over the determinants per state plus one projector per multi-member
family.

Fast labels for Davidson
------------------------

Davidson cannot afford the exact family projectors at every iteration. It
uses ``pg_det_state_labels``, which sums ``coef(i,k)**2`` over the
determinants of each irrep and over the unlabelled determinants, and returns
a relative purity and a label. Everything is a ratio because Davidson hands
in vectors that are not normalized. A state whose weight sits on
multi-member families gets label 0 and is not constrained; the exact
``pg_state_irrep_weight`` remains the value the report prints.

The labels are used to follow states: the label of the current vector is
compared with the label of the previous wave function, so state following
does not swap two states of different irreps.

The ``irrep_states`` filter
===========================

:option:`symmetry irrep_states` is a list of irrep names. When it is non-empty and
the symmetry machinery is usable, the diagonalization keeps only the states
whose known irrep is in the list. The semantics are deliberately narrow:

* it is applied **after** the :math:`S^2` selection and, in Davidson, after
  the state-following attribution. It never widens a selection; it is a
  filter on the states that already survived;
* a state with no definite irrep (label 0) is always kept: a filter must not
  drop a state it cannot judge;
* the kept states are moved to the front of the CI vector set, preserving
  their energy order;
* when fewer states than requested are allowed, the run completes and says
  so: ``symmetry: irrep filter: only N state(s) in irrep_states for n_states
  = M``. This message is printed once per diagonalization call, so a CI run
  with several calls shows it several times;
* every state whose known irrep is rejected prints
  ``symmetry: irrep filter: state <k> is <name>, not in irrep_states``;
* :option:`symmetry state_symmetry` does not disable the filter: a printing switch
  must not silently change a computation;
* CSF-based Davidson (``csf_based = true``) cannot be filtered. It prints
  ``symmetry: irrep_states is ignored when csf_based is true`` and runs
  unchanged;
* an unknown irrep name prints the available names and stops: malformed
  input is the one place where stopping is allowed.

On the Lapack path (small determinant spaces, ``N_det`` at most
``n_det_max_full``), the filter is a *stable permutation* of the list that
the :math:`S^2` selection already produced. That list is "states of the
expected :math:`S^2`, then fill states". The filter promotes the allowed
fill states too, so asking for ``A1`` with ``n_states = 3`` can return two
singlets and a triplet: it is the third ``A1`` state of the computed slots,
not necessarily the third singlet. The energies are unchanged; only the
order and the count of kept states are affected.

State following and the filter compose as "keep the states that follow,
then move the allowed ones forward". The filter therefore acts on the states
that passed the :math:`S^2` test, not on the whole subspace; when the
requested irreps are not among them, raise
:option:`davidson_keywords n_states_diag` (the multiplier of
:option:`determinants n_states` that sets the number of computed roots) or
:option:`determinants n_states`, and re-run.

Report program
==============

``qp_run symmetry_det DIR`` prints the families, the determinant table (the
individual lines are hidden above 200 determinants), the per-irrep weights
of every state and the gates D1-D6. It exits 1 if a gate fails, after every
number has been printed, so it can be chained in a test script.

The program reads the wave function stored in the EZFIO. It is meaningful
after a CI calculation; after an SCF only, it labels the Hartree-Fock
determinant.

Gates
=====

.. list-table::
   :header-rows: 1

   * - Gate
     - Test
   * - D1
     - :math:`\sum_a T^a = 1` on every family (validates the characters,
       :math:`m_a` and the family structure at once)
   * - D2
     - :math:`(T^a)^2 = T^a` and :math:`T^a` Hermitian, every family,
       every irrep
   * - D3
     - :math:`-10^{-8} \le` purity :math:`\le 1 + 10^{-8}` for every
       determinant
   * - D4
     - :math:`|\sum_a w_a(s) - \sum_D c_D(s)^2| \le 10^{-8}` (the
       bookkeeping sum rule)
   * - D5
     - every labelled state has :math:`w_{\max} \ge 1 - 10^{-6}`; the mixed
       determinant fraction is reported
   * - D6
     - the dominant determinant of state 1 has irrep 1 for a closed-shell
       Hartree-Fock wave function, when that determinant is pure

Families of size one are skipped by D1 and D2: there :math:`T^a` is a
Kronecker delta by construction, and the label already comes from it.

Physics note
============

An open degenerate shell has no one-determinant irrep, and neither has a
closed-shell occupation of a partially filled shell: with one electron in a
two-dimensional ``e`` block, the determinant is a 50/50 mixture of two
irreps. That is the SCF solution breaking the point group, not a labelling
error. The determinant is reported with a purity below 1 and no label; the
state weights carry the symmetry instead. On an Abelian group every
determinant is pure, whatever its open shells.

Limits
======

* a determinant family is enumerated only up to ``pg_family_cap`` (1000)
  members;
* a symmetry block may hold at most 12 occupied MOs in the determinant
  projector;
* when a limit is hit, ``pg_label_ierr`` and ``pg_label_msg`` say which
  determinant was refused, and the CI programs run without symmetry labels.
  Nothing hangs and nothing is silently truncated.
 
 
 
Programs 
-------- 
 
 * :ref:`symmetry_det` 
 
Providers 
--------- 
 
.. c:var:: int pg_block_dim


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        integer, allocatable	:: pg_block_start	(pg_n_mo_block)
        integer, allocatable	:: pg_block_mo	(mo_num)
        integer, allocatable	:: pg_block_dim	(pg_n_mo_block)


    ``pg_block_mo(pg_block_start(b) : pg_block_start(b)+pg_block_dim(b)-1)``
    holds the MOs of symmetry block ``b`` in ascending MO order.
    
    The MOs of one block are NOT contiguous in the MO list (the same
    irrep can appear at several energies: the A1g MOs of benzene are
    spread over the whole spectrum), which is why the family enumeration
    of section 9.2 goes through this indirection and never through
    ``pg_mo_block(i)`` arithmetic.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`pg_col_diag`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_det_irrep`
       * :c:data:`pg_state_irrep_weight`

 
.. c:var:: int pg_block_mo


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        integer, allocatable	:: pg_block_start	(pg_n_mo_block)
        integer, allocatable	:: pg_block_mo	(mo_num)
        integer, allocatable	:: pg_block_dim	(pg_n_mo_block)


    ``pg_block_mo(pg_block_start(b) : pg_block_start(b)+pg_block_dim(b)-1)``
    holds the MOs of symmetry block ``b`` in ascending MO order.
    
    The MOs of one block are NOT contiguous in the MO list (the same
    irrep can appear at several energies: the A1g MOs of benzene are
    spread over the whole spectrum), which is why the family enumeration
    of section 9.2 goes through this indirection and never through
    ``pg_mo_block(i)`` arithmetic.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`pg_col_diag`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_det_irrep`
       * :c:data:`pg_state_irrep_weight`

 
.. c:var:: int pg_block_start


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        integer, allocatable	:: pg_block_start	(pg_n_mo_block)
        integer, allocatable	:: pg_block_mo	(mo_num)
        integer, allocatable	:: pg_block_dim	(pg_n_mo_block)


    ``pg_block_mo(pg_block_start(b) : pg_block_start(b)+pg_block_dim(b)-1)``
    holds the MOs of symmetry block ``b`` in ascending MO order.
    
    The MOs of one block are NOT contiguous in the MO list (the same
    irrep can appear at several energies: the A1g MOs of benzene are
    spread over the whole spectrum), which is why the family enumeration
    of section 9.2 goes through this indirection and never through
    ``pg_mo_block(i)`` arithmetic.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`pg_col_diag`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_det_irrep`
       * :c:data:`pg_state_irrep_weight`

 
.. c:var:: int pg_det_block


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        integer, allocatable	:: pg_det_irrep	(N_det)
        double precision, allocatable	:: pg_det_purity	(N_det)
        integer, allocatable	:: pg_det_block	(N_det)
        integer, allocatable	:: pg_det_family	(N_det)
        integer	:: pg_n_family	
        integer, allocatable	:: pg_family_first	(N_det)


    ``pg_det_irrep(D)`` = argmax_a T^a_{DD}, or 0 when max_a T^a_{DD} <
    1 - pg_thresh: a determinant with one electron in a degenerate block
    (:math:`|e_x>`, :math:`|e_y>`) is not an eigenfunction of the projectors and must not
    get a best-guess label (SPEC.md 9.3).
    
    ``pg_det_purity(D)`` = max_a T^a_{DD} = sum_{D'} :math:`(T^a_{D'D})^2` at the
    argmax, 1 for a pure determinant and 1/d_a for one member of an open
    degenerate shell.
    
    ``pg_det_block(D)`` = :math:`|F|`, the dimension of the family of D;
    ``pg_det_family(D)`` = its family id, in order of first appearance;
    ``pg_family_first(f)`` = the psi_det index of the first member of
    family f, i.e. the representative the other modules enumerate.
    
    A determinant is pure if and only if its family has size 1: a family
    of size > 1 exists only because a block has 1 <= n <= d_b-1, and then
    every member shares that partial occupation. The state weights of
    section 10 use this to skip the projectors for pure determinants.
    
    :c:data:`pg_n_family` <= N_det. A determinant whose family is handled by an
    earlier member costs one enumeration and one binary search per member
    and nothing else.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`pg_block_start`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_elem_col`
       * :c:data:`pg_family_cap`
       * :c:data:`pg_label_ierr`
       * :c:data:`pg_mo_weight`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_thresh`
       * :c:data:`psi_det`
       * :c:data:`psi_det_sorted_bit`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_state_irrep_weight`

 
.. c:var:: int pg_det_family


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        integer, allocatable	:: pg_det_irrep	(N_det)
        double precision, allocatable	:: pg_det_purity	(N_det)
        integer, allocatable	:: pg_det_block	(N_det)
        integer, allocatable	:: pg_det_family	(N_det)
        integer	:: pg_n_family	
        integer, allocatable	:: pg_family_first	(N_det)


    ``pg_det_irrep(D)`` = argmax_a T^a_{DD}, or 0 when max_a T^a_{DD} <
    1 - pg_thresh: a determinant with one electron in a degenerate block
    (:math:`|e_x>`, :math:`|e_y>`) is not an eigenfunction of the projectors and must not
    get a best-guess label (SPEC.md 9.3).
    
    ``pg_det_purity(D)`` = max_a T^a_{DD} = sum_{D'} :math:`(T^a_{D'D})^2` at the
    argmax, 1 for a pure determinant and 1/d_a for one member of an open
    degenerate shell.
    
    ``pg_det_block(D)`` = :math:`|F|`, the dimension of the family of D;
    ``pg_det_family(D)`` = its family id, in order of first appearance;
    ``pg_family_first(f)`` = the psi_det index of the first member of
    family f, i.e. the representative the other modules enumerate.
    
    A determinant is pure if and only if its family has size 1: a family
    of size > 1 exists only because a block has 1 <= n <= d_b-1, and then
    every member shares that partial occupation. The state weights of
    section 10 use this to skip the projectors for pure determinants.
    
    :c:data:`pg_n_family` <= N_det. A determinant whose family is handled by an
    earlier member costs one enumeration and one binary search per member
    and nothing else.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`pg_block_start`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_elem_col`
       * :c:data:`pg_family_cap`
       * :c:data:`pg_label_ierr`
       * :c:data:`pg_mo_weight`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_thresh`
       * :c:data:`psi_det`
       * :c:data:`psi_det_sorted_bit`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_state_irrep_weight`

 
.. c:var:: int pg_det_irrep


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        integer, allocatable	:: pg_det_irrep	(N_det)
        double precision, allocatable	:: pg_det_purity	(N_det)
        integer, allocatable	:: pg_det_block	(N_det)
        integer, allocatable	:: pg_det_family	(N_det)
        integer	:: pg_n_family	
        integer, allocatable	:: pg_family_first	(N_det)


    ``pg_det_irrep(D)`` = argmax_a T^a_{DD}, or 0 when max_a T^a_{DD} <
    1 - pg_thresh: a determinant with one electron in a degenerate block
    (:math:`|e_x>`, :math:`|e_y>`) is not an eigenfunction of the projectors and must not
    get a best-guess label (SPEC.md 9.3).
    
    ``pg_det_purity(D)`` = max_a T^a_{DD} = sum_{D'} :math:`(T^a_{D'D})^2` at the
    argmax, 1 for a pure determinant and 1/d_a for one member of an open
    degenerate shell.
    
    ``pg_det_block(D)`` = :math:`|F|`, the dimension of the family of D;
    ``pg_det_family(D)`` = its family id, in order of first appearance;
    ``pg_family_first(f)`` = the psi_det index of the first member of
    family f, i.e. the representative the other modules enumerate.
    
    A determinant is pure if and only if its family has size 1: a family
    of size > 1 exists only because a block has 1 <= n <= d_b-1, and then
    every member shares that partial occupation. The state weights of
    section 10 use this to skip the projectors for pure determinants.
    
    :c:data:`pg_n_family` <= N_det. A determinant whose family is handled by an
    earlier member costs one enumeration and one binary search per member
    and nothing else.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`pg_block_start`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_elem_col`
       * :c:data:`pg_family_cap`
       * :c:data:`pg_label_ierr`
       * :c:data:`pg_mo_weight`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_thresh`
       * :c:data:`psi_det`
       * :c:data:`psi_det_sorted_bit`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_state_irrep_weight`

 
.. c:var:: double pg_det_mixed_frac


    File : :file:`state_weights.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pg_state_irrep_weight	(pg_n_irrep,N_states)
        integer, allocatable	:: pg_state_irrep	(N_states)
        double precision, allocatable	:: pg_state_purity	(N_states)
        double precision, allocatable	:: pg_det_mixed_frac	(N_states)
        double precision, allocatable	:: pg_state_norm	(N_states)


    ``pg_state_irrep_weight(a,s)`` = w_a(s), the weight of irrep ``a`` in
    state ``s``: the diagonal matrix element of the projector of section 8
    in the state, not an occupation-number guess.
    
    ``pg_state_irrep(s)`` = argmax_a w_a(s), or 0 when the state has no
    definite irrep (spin contamination, an incomplete family, a broken
    MO basis): a best guess would be a lie.
    ``pg_state_purity(s)`` = max_a w_a(s).
    
    ``pg_det_mixed_frac(s)`` = sum of c_D(s)^2 over the determinants whose
    :c:data:`pg_det_irrep` is 0: the weight of the state carried by determinants
    that belong to no symmetry group on their own.
    
    ``pg_state_norm(s)`` = sum_D c_D(s)^2, the exact value of
    sum_a w_a(s) (SPEC.md 10.2). A deficit means psi_coef is not
    normalized in the determinant basis; it never means that the CI space
    is not group invariant.
    
    A determinant with :c:data:`pg_det_block` = 1 is pure and its family is
    itself, so T^a_{DD} = delta_{a,pg_det_irrep(D)} exactly: those
    determinants are accumulated from the label, and the projectors are
    built only for the families of size > 1 (which are exactly the
    determinants with no definite irrep).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`pg_block_start`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_elem_col`
       * :c:data:`pg_label_ierr`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_thresh`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_sorted_bit`


 
.. c:var:: double pg_det_purity


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        integer, allocatable	:: pg_det_irrep	(N_det)
        double precision, allocatable	:: pg_det_purity	(N_det)
        integer, allocatable	:: pg_det_block	(N_det)
        integer, allocatable	:: pg_det_family	(N_det)
        integer	:: pg_n_family	
        integer, allocatable	:: pg_family_first	(N_det)


    ``pg_det_irrep(D)`` = argmax_a T^a_{DD}, or 0 when max_a T^a_{DD} <
    1 - pg_thresh: a determinant with one electron in a degenerate block
    (:math:`|e_x>`, :math:`|e_y>`) is not an eigenfunction of the projectors and must not
    get a best-guess label (SPEC.md 9.3).
    
    ``pg_det_purity(D)`` = max_a T^a_{DD} = sum_{D'} :math:`(T^a_{D'D})^2` at the
    argmax, 1 for a pure determinant and 1/d_a for one member of an open
    degenerate shell.
    
    ``pg_det_block(D)`` = :math:`|F|`, the dimension of the family of D;
    ``pg_det_family(D)`` = its family id, in order of first appearance;
    ``pg_family_first(f)`` = the psi_det index of the first member of
    family f, i.e. the representative the other modules enumerate.
    
    A determinant is pure if and only if its family has size 1: a family
    of size > 1 exists only because a block has 1 <= n <= d_b-1, and then
    every member shares that partial occupation. The state weights of
    section 10 use this to skip the projectors for pure determinants.
    
    :c:data:`pg_n_family` <= N_det. A determinant whose family is handled by an
    earlier member costs one enumeration and one binary search per member
    and nothing else.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`pg_block_start`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_elem_col`
       * :c:data:`pg_family_cap`
       * :c:data:`pg_label_ierr`
       * :c:data:`pg_mo_weight`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_thresh`
       * :c:data:`psi_det`
       * :c:data:`psi_det_sorted_bit`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_state_irrep_weight`

 
.. c:var:: int pg_elem_col


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        integer, allocatable	:: pg_elem_col	(pg_n_elem)


    Inverse of :c:data:`pg_col_member`: ``pg_elem_col(g)`` is the file column that
    element ``g`` belongs to. :c:data:`pg_char` is indexed by column and both
    projectors sum over elements, so this is the map between the two.
    
    A column can stand for several geometric classes (T, Td, Oh, Ih have
    classes of 4, 6, 8, ... elements); the character is constant on the
    column, which is why section 8 sums over columns and not over classes.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_n_col`
       * :c:data:`pg_n_elem`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_det_irrep`
       * :c:data:`pg_state_irrep_weight`

 
.. c:var:: int pg_family_cap


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        integer	:: pg_family_cap	


    Largest family the module agrees to enumerate. A family has
    ``prod_b C(d_b,n_b)`` members per spin, so a determinant with several
    partially occupied degenerate shells grows combinatorially. Above
    this size the module stops and says which determinant it was,
    instead of hanging or silently truncating the projector.

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_det_irrep`

 
.. c:var:: int pg_family_first


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        integer, allocatable	:: pg_det_irrep	(N_det)
        double precision, allocatable	:: pg_det_purity	(N_det)
        integer, allocatable	:: pg_det_block	(N_det)
        integer, allocatable	:: pg_det_family	(N_det)
        integer	:: pg_n_family	
        integer, allocatable	:: pg_family_first	(N_det)


    ``pg_det_irrep(D)`` = argmax_a T^a_{DD}, or 0 when max_a T^a_{DD} <
    1 - pg_thresh: a determinant with one electron in a degenerate block
    (:math:`|e_x>`, :math:`|e_y>`) is not an eigenfunction of the projectors and must not
    get a best-guess label (SPEC.md 9.3).
    
    ``pg_det_purity(D)`` = max_a T^a_{DD} = sum_{D'} :math:`(T^a_{D'D})^2` at the
    argmax, 1 for a pure determinant and 1/d_a for one member of an open
    degenerate shell.
    
    ``pg_det_block(D)`` = :math:`|F|`, the dimension of the family of D;
    ``pg_det_family(D)`` = its family id, in order of first appearance;
    ``pg_family_first(f)`` = the psi_det index of the first member of
    family f, i.e. the representative the other modules enumerate.
    
    A determinant is pure if and only if its family has size 1: a family
    of size > 1 exists only because a block has 1 <= n <= d_b-1, and then
    every member shares that partial occupation. The state weights of
    section 10 use this to skip the projectors for pure determinants.
    
    :c:data:`pg_n_family` <= N_det. A determinant whose family is handled by an
    earlier member costs one enumeration and one binary search per member
    and nothing else.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`pg_block_start`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_elem_col`
       * :c:data:`pg_family_cap`
       * :c:data:`pg_label_ierr`
       * :c:data:`pg_mo_weight`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_thresh`
       * :c:data:`psi_det`
       * :c:data:`psi_det_sorted_bit`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_state_irrep_weight`

 
.. c:var:: int pg_label_ierr


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        integer	:: pg_label_ierr	
        character*(512)	:: pg_label_msg	


    0 when the last labelling attempt succeeded; otherwise
    
    ``1``
        there is no determinant to label
    
    ``2``
        there is no MO or no symmetry block
    
    ``3``
        the MO basis is not block diagonal
    
    ``4``
        a determinant family is larger than pg_family_cap
    
    ``5``
        the family enumeration is inconsistent
    
    ``6``
        one symmetry block holds more than 12 occupied MOs
    
    :c:data:`pg_label_msg` is the one-sentence reason. Reporting is left to the
    caller: a provider that stops cannot be reported on.

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_det_irrep`
       * :c:data:`pg_state_irrep_weight`

 
.. c:var:: char pg_label_msg


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        integer	:: pg_label_ierr	
        character*(512)	:: pg_label_msg	


    0 when the last labelling attempt succeeded; otherwise
    
    ``1``
        there is no determinant to label
    
    ``2``
        there is no MO or no symmetry block
    
    ``3``
        the MO basis is not block diagonal
    
    ``4``
        a determinant family is larger than pg_family_cap
    
    ``5``
        the family enumeration is inconsistent
    
    ``6``
        one symmetry block holds more than 12 occupied MOs
    
    :c:data:`pg_label_msg` is the one-sentence reason. Reporting is left to the
    caller: a provider that stops cannot be reported on.

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_det_irrep`
       * :c:data:`pg_state_irrep_weight`

 
.. c:var:: int pg_n_family


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        integer, allocatable	:: pg_det_irrep	(N_det)
        double precision, allocatable	:: pg_det_purity	(N_det)
        integer, allocatable	:: pg_det_block	(N_det)
        integer, allocatable	:: pg_det_family	(N_det)
        integer	:: pg_n_family	
        integer, allocatable	:: pg_family_first	(N_det)


    ``pg_det_irrep(D)`` = argmax_a T^a_{DD}, or 0 when max_a T^a_{DD} <
    1 - pg_thresh: a determinant with one electron in a degenerate block
    (:math:`|e_x>`, :math:`|e_y>`) is not an eigenfunction of the projectors and must not
    get a best-guess label (SPEC.md 9.3).
    
    ``pg_det_purity(D)`` = max_a T^a_{DD} = sum_{D'} :math:`(T^a_{D'D})^2` at the
    argmax, 1 for a pure determinant and 1/d_a for one member of an open
    degenerate shell.
    
    ``pg_det_block(D)`` = :math:`|F|`, the dimension of the family of D;
    ``pg_det_family(D)`` = its family id, in order of first appearance;
    ``pg_family_first(f)`` = the psi_det index of the first member of
    family f, i.e. the representative the other modules enumerate.
    
    A determinant is pure if and only if its family has size 1: a family
    of size > 1 exists only because a block has 1 <= n <= d_b-1, and then
    every member shares that partial occupation. The state weights of
    section 10 use this to skip the projectors for pure determinants.
    
    :c:data:`pg_n_family` <= N_det. A determinant whose family is handled by an
    earlier member costs one enumeration and one binary search per member
    and nothing else.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`pg_block_start`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_elem_col`
       * :c:data:`pg_family_cap`
       * :c:data:`pg_label_ierr`
       * :c:data:`pg_mo_weight`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_thresh`
       * :c:data:`psi_det`
       * :c:data:`psi_det_sorted_bit`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_state_irrep_weight`

 
.. c:var:: int pg_state_irrep


    File : :file:`state_weights.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pg_state_irrep_weight	(pg_n_irrep,N_states)
        integer, allocatable	:: pg_state_irrep	(N_states)
        double precision, allocatable	:: pg_state_purity	(N_states)
        double precision, allocatable	:: pg_det_mixed_frac	(N_states)
        double precision, allocatable	:: pg_state_norm	(N_states)


    ``pg_state_irrep_weight(a,s)`` = w_a(s), the weight of irrep ``a`` in
    state ``s``: the diagonal matrix element of the projector of section 8
    in the state, not an occupation-number guess.
    
    ``pg_state_irrep(s)`` = argmax_a w_a(s), or 0 when the state has no
    definite irrep (spin contamination, an incomplete family, a broken
    MO basis): a best guess would be a lie.
    ``pg_state_purity(s)`` = max_a w_a(s).
    
    ``pg_det_mixed_frac(s)`` = sum of c_D(s)^2 over the determinants whose
    :c:data:`pg_det_irrep` is 0: the weight of the state carried by determinants
    that belong to no symmetry group on their own.
    
    ``pg_state_norm(s)`` = sum_D c_D(s)^2, the exact value of
    sum_a w_a(s) (SPEC.md 10.2). A deficit means psi_coef is not
    normalized in the determinant basis; it never means that the CI space
    is not group invariant.
    
    A determinant with :c:data:`pg_det_block` = 1 is pure and its family is
    itself, so T^a_{DD} = delta_{a,pg_det_irrep(D)} exactly: those
    determinants are accumulated from the label, and the projectors are
    built only for the families of size > 1 (which are exactly the
    determinants with no definite irrep).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`pg_block_start`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_elem_col`
       * :c:data:`pg_label_ierr`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_thresh`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_sorted_bit`


 
.. c:var:: double pg_state_irrep_weight


    File : :file:`state_weights.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pg_state_irrep_weight	(pg_n_irrep,N_states)
        integer, allocatable	:: pg_state_irrep	(N_states)
        double precision, allocatable	:: pg_state_purity	(N_states)
        double precision, allocatable	:: pg_det_mixed_frac	(N_states)
        double precision, allocatable	:: pg_state_norm	(N_states)


    ``pg_state_irrep_weight(a,s)`` = w_a(s), the weight of irrep ``a`` in
    state ``s``: the diagonal matrix element of the projector of section 8
    in the state, not an occupation-number guess.
    
    ``pg_state_irrep(s)`` = argmax_a w_a(s), or 0 when the state has no
    definite irrep (spin contamination, an incomplete family, a broken
    MO basis): a best guess would be a lie.
    ``pg_state_purity(s)`` = max_a w_a(s).
    
    ``pg_det_mixed_frac(s)`` = sum of c_D(s)^2 over the determinants whose
    :c:data:`pg_det_irrep` is 0: the weight of the state carried by determinants
    that belong to no symmetry group on their own.
    
    ``pg_state_norm(s)`` = sum_D c_D(s)^2, the exact value of
    sum_a w_a(s) (SPEC.md 10.2). A deficit means psi_coef is not
    normalized in the determinant basis; it never means that the CI space
    is not group invariant.
    
    A determinant with :c:data:`pg_det_block` = 1 is pure and its family is
    itself, so T^a_{DD} = delta_{a,pg_det_irrep(D)} exactly: those
    determinants are accumulated from the label, and the projectors are
    built only for the families of size > 1 (which are exactly the
    determinants with no definite irrep).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`pg_block_start`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_elem_col`
       * :c:data:`pg_label_ierr`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_thresh`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_sorted_bit`


 
.. c:var:: double pg_state_norm


    File : :file:`state_weights.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pg_state_irrep_weight	(pg_n_irrep,N_states)
        integer, allocatable	:: pg_state_irrep	(N_states)
        double precision, allocatable	:: pg_state_purity	(N_states)
        double precision, allocatable	:: pg_det_mixed_frac	(N_states)
        double precision, allocatable	:: pg_state_norm	(N_states)


    ``pg_state_irrep_weight(a,s)`` = w_a(s), the weight of irrep ``a`` in
    state ``s``: the diagonal matrix element of the projector of section 8
    in the state, not an occupation-number guess.
    
    ``pg_state_irrep(s)`` = argmax_a w_a(s), or 0 when the state has no
    definite irrep (spin contamination, an incomplete family, a broken
    MO basis): a best guess would be a lie.
    ``pg_state_purity(s)`` = max_a w_a(s).
    
    ``pg_det_mixed_frac(s)`` = sum of c_D(s)^2 over the determinants whose
    :c:data:`pg_det_irrep` is 0: the weight of the state carried by determinants
    that belong to no symmetry group on their own.
    
    ``pg_state_norm(s)`` = sum_D c_D(s)^2, the exact value of
    sum_a w_a(s) (SPEC.md 10.2). A deficit means psi_coef is not
    normalized in the determinant basis; it never means that the CI space
    is not group invariant.
    
    A determinant with :c:data:`pg_det_block` = 1 is pure and its family is
    itself, so T^a_{DD} = delta_{a,pg_det_irrep(D)} exactly: those
    determinants are accumulated from the label, and the projectors are
    built only for the families of size > 1 (which are exactly the
    determinants with no definite irrep).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`pg_block_start`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_elem_col`
       * :c:data:`pg_label_ierr`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_thresh`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_sorted_bit`


 
.. c:var:: double pg_state_purity


    File : :file:`state_weights.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pg_state_irrep_weight	(pg_n_irrep,N_states)
        integer, allocatable	:: pg_state_irrep	(N_states)
        double precision, allocatable	:: pg_state_purity	(N_states)
        double precision, allocatable	:: pg_det_mixed_frac	(N_states)
        double precision, allocatable	:: pg_state_norm	(N_states)


    ``pg_state_irrep_weight(a,s)`` = w_a(s), the weight of irrep ``a`` in
    state ``s``: the diagonal matrix element of the projector of section 8
    in the state, not an occupation-number guess.
    
    ``pg_state_irrep(s)`` = argmax_a w_a(s), or 0 when the state has no
    definite irrep (spin contamination, an incomplete family, a broken
    MO basis): a best guess would be a lie.
    ``pg_state_purity(s)`` = max_a w_a(s).
    
    ``pg_det_mixed_frac(s)`` = sum of c_D(s)^2 over the determinants whose
    :c:data:`pg_det_irrep` is 0: the weight of the state carried by determinants
    that belong to no symmetry group on their own.
    
    ``pg_state_norm(s)`` = sum_D c_D(s)^2, the exact value of
    sum_a w_a(s) (SPEC.md 10.2). A deficit means psi_coef is not
    normalized in the determinant basis; it never means that the CI space
    is not group invariant.
    
    A determinant with :c:data:`pg_det_block` = 1 is pure and its family is
    itself, so T^a_{DD} = delta_{a,pg_det_irrep(D)} exactly: those
    determinants are accumulated from the label, and the projectors are
    built only for the families of size > 1 (which are exactly the
    determinants with no definite irrep).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`pg_block_start`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_elem_col`
       * :c:data:`pg_label_ierr`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_thresh`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_sorted_bit`


 
.. c:var:: double pg_xg


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pg_xg	(mo_num,mo_num,pg_n_elem)


    ``pg_xg(:,:,g)`` = C^T S M_g C, the matrix of group element ``g`` in the
    MO basis (SPEC.md section 7), built once.
    
    SPEC.md 9.4 asked for one scratch matrix at a time, but the
    determinant projectors need every element for every family, and a
    rebuild is a full AO transform (three matmuls over ao_num). :math:`|G|` is at
    most 120 and the MO basis of a plugin-sized system a few hundred
    wide, so the set costs :math:`8|G|N_{MO}^2` bytes: 0.8 MB for benzene in
    6-31G. The report prints the number. A production version can keep
    only the rows and columns of the partially occupied blocks; nothing
    below cares, because :c:func:`pg_elem_me` is the only reader.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`pg_n_elem`


 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: void pg_binom(void)


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        integer function pg_binom(dim, n)


    C(dim,n), 0 outside 0<=n<=dim. The product is built so that every
    division is exact (r*(dim-n+k)/k is an integer at each step), which
    is what keeps C(5,2) = 10 and C(12,6) = 924 worth of exactness in a
    four-byte integer.

 
.. c:function:: void pg_det_counts(void)


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        subroutine pg_det_counts(idet, nsel)


    ``nsel(s,b)`` = number of spin-``s`` electrons of determinant ``idet`` in
    symmetry block ``b``. Two determinants are in the same family if and
    only if their ``nsel`` agree (SPEC.md 9.2).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`pg_col_diag`
       * :c:data:`psi_det`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_det_gates`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_state_irrep_weight`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`

 
.. c:function:: void pg_det_gates(void)


    File : :file:`symmetry_det.irp.f`

    .. code:: fortran

        logical function pg_det_gates()


    ``D1``
        sum_a T^a = 1 on every family (validates chi, mult_a and the
        family structure at once).
    
    ``D2``
        (T^a)^2 = T^a and T^a Hermitian, every family, every a.
    
    ``D3``
        -1e-8 <= purity <= 1+1e-8 for every determinant.
    
    ``D4``
        ``sum_a w_a(s) - sum_D c_D(s)^2`` <= 1e-8: the accumulation loop
        reproduces the sum rule of section 10.2.
    
    ``D5``
        every labelled state has w_max >= 1-1e-6, and the mixed
        determinant fraction is reported.
    
    ``D6``
        the dominant determinant of state 1 has irrep 1 for a closed-shell
        Hartree-Fock wave function (the cross-check that section 9 and
        section 8 live in the same labelling). Gated only when the
        determinant's family has one member: a block that is partially
        filled makes even a closed-shell occupation a mixture, so such a
        configuration has no one-determinant irrep to compare with. That
        case is printed, and the state weights carry the answer.
    
    Families of size 1 are skipped by D1/D2: there T^a = delta is exact
    by construction, and the label already comes from that diagonal.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_n_col`
       * :c:data:`pg_state_irrep_weight`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pg_det_counts`
       * :c:func:`pg_family_fill`
       * :c:func:`pg_family_projector`

 
.. c:function:: void pg_det_report(void)


    File : :file:`symmetry_det.irp.f`

    .. code:: fortran

        subroutine pg_det_report


    Prints the labels and the gates, and stops with status 1 if any gate
    fails. Everything measured is printed before the verdict: one run
    shows the whole story.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_detection_failed`
       * :c:data:`pg_family_cap`
       * :c:data:`pg_group_name`
       * :c:data:`pg_label_ierr`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_state_irrep_weight`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`read_wf`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`symmetry_det`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`

 
.. c:function:: void pg_det_state_labels(void)


    File : :file:`state_labels.irp.f`

    .. code:: fortran

        subroutine pg_det_state_labels(ndet, coef, lda, nst, label, purity, mixed)


    Fast irrep labels of ``nst`` CI vectors stored in ``coef``.
    
    For each column k: ``purity(k)`` = wmax/tot, ``mixed(k)`` = wmix/tot and
    ``label(k)`` = argmax when wmax >= (1 - pg_thresh)*tot, else 0, where
    w(a) sums coef(i,k)**2 over the determinants of irrep a, and wmix sums
    the determinants whose :c:data:`pg_det_irrep` is 0 (mixed families).
    
    Everything is a ratio: the callers hand in vectors that are not
    normalized (Davidson's in-loop ``u_in`` is W*y), so an absolute test
    would label nothing. A zero column gives label 0, purity 0, mixed 0.
    
    The output arrays are intent(inout), not intent(out): they are
    assumed-size dummies, and an intent(out) dummy would leave the whole
    actual argument undefined whatever bounds the caller used. The per-k
    zeroing below defines them.
    
    Fast path, not the exact formula: a state whose weight sits on
    multi-member determinant families gets label 0 and is not constrained.
    :c:data:`pg_state_irrep_weight` (stage 3) stays the exact value, off-diagonal
    family terms included, and is what the report prints.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_det_irrep`
       * :c:data:`pg_label_ierr`
       * :c:data:`pg_n_col`
       * :c:data:`pg_thresh`

 
.. c:function:: void pg_elem_me(void)


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        subroutine pg_elem_me(ig, np, mask, ip, jp, val)


    ``val`` = :math:`<D_ip|U(g_ig)|D_jp>` = prod_b det( X_g[A'_b,A_b] ) * prod_b det( X_g[B'_b,B_b] ).
    
    Every block contributes a square submatrix because the family has
    the same count in every block. The counts are computed first and the
    submatrix filled afterwards: the two must not share a running index,
    or a block with a zero count inherits the previous block's value.
    
    An over-full or inconsistent pair returns exactly 0 rather than a
    number, and an exact zero is propagated exactly, so the caller may
    skip the accumulation of that pair.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_block_start`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_label_ierr`
       * :c:data:`pg_xg`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_family_projector`

 
.. c:function:: void pg_family_fill(void)


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        subroutine pg_family_fill(nsel, np, mask)


    ``mask(s,b,ip)`` is the bit mask of the MOs of block ``b`` occupied with
    spin ``s`` by family member ``ip``, bit ``k`` standing for the
    MO ``pg_block_mo(pg_block_start(b)+k)``.
    
    An odometer over the per-(spin,block) combination ranks: ``scalar``
    index q in [0,np-1] decomposed in mixed radix, each digit unranked to
    a subset. Cost np * sum_b d_b, paid once per family.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_block_start`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_label_ierr`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_det_gates`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_state_irrep_weight`

 
.. c:function:: void pg_family_index(void)


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        subroutine pg_family_index(np, mask, idx, imin)


    ``idx(ip)`` is the index in :c:data:`psi_det` of family member ``ip``,
    or 0 when that member is not in the wave function (a truncated CI
    space: the missing member simply has coefficient 0, SPEC.md 10.2).
    ``imin`` is the smallest nonzero index.
    
    The lookup is QP2's own binary search on :c:data:`psi_det_sorted_bit`, which
    is why the members are built as bit strings through :c:func:`list_to_bitstring`
    and not compared determinantal by determinantal.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`pg_block_start`
       * :c:data:`pg_col_diag`
       * :c:data:`psi_det_sorted_bit`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_det_irrep`
       * :c:data:`pg_state_irrep_weight`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`list_to_bitstring`

 
.. c:function:: void pg_family_projector(void)


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        subroutine pg_family_projector(np, mask, t)


    ``t(ip,jp,a)`` = T^a_{ip,jp} = (d_a/(mult_a*:math:`|G|`)) sum_g chi_a(col(g)) :math:`<D_ip|U(g)|D_jp>`.
    
    This is SPEC.md 9.1 verbatim: the same coefficients, the same
    characters, the same ``mult_a`` as :c:data:`pg_mo_weight`. T^a is a Hermitian
    projector on the family (gates D1/D2), and its diagonal gives the
    determinant labels.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_col_diag`
       * :c:data:`pg_elem_col`
       * :c:data:`pg_label_ierr`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_elem`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_det_gates`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_state_irrep_weight`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pg_elem_me`

 
.. c:function:: void pg_family_size(void)


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        subroutine pg_family_size(nsel, np)


    ``np`` = prod_{s,b} C(d_b, nsel(s,b)), the number of determinants with
    the same per-block spin counts. A block that is empty or full
    contributes one choice; the product is cut as soon as it exceeds
    :c:data:`pg_family_cap` so that the caller reports instead of overflowing.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_block_start`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_family_cap`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_det_irrep`

 
.. c:function:: void pg_labels_usable(void)


    File : :file:`state_labels.irp.f`

    .. code:: fortran

        logical function pg_labels_usable()


    Cheap guard for the stage-4 symmetry work: true when there is a
    determinant, the point group was detected and the MOs carry definite
    irreps. Never stops and never builds the determinant labelling: it
    names ``N_det`` and the two stage-2 functions only, so a caller can test
    it before paying for anything.
    
    It does not test :c:data:`state_symmetry`: that keyword means "print the
    table", and a printing switch must not disable the irrep filter
    (SPEC.md 9.4f A1). :c:data:`auto_symmetry` stays the master switch, through
    :c:func:`pg_symmetry_ready`.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`

 
.. c:function:: void pg_print_state_irreps(void)


    File : :file:`state_report.irp.f`

    .. code:: fortran

        subroutine pg_print_state_irreps(energy, s2)


    Prints the irrep of each CI state. Never stops. :c:data:`state_symmetry` = F
    is silent; a failed detection has already printed its reason; a mixed
    MO basis is reported in one line. ``energy`` and ``s2`` are the caller's
    arrays, at least ``N_states`` long.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`mpi_master`
       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`state_symmetry`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pg_state_irreps_label`

 
.. c:function:: void pg_print_state_irreps_work(void)


    File : :file:`state_report.irp.f`

    .. code:: fortran

        subroutine pg_print_state_irreps_work(energy, s2)


    The table. :c:data:`pg_state_irrep_weight` delivers the weights, the labels,
    the purities and the norm in one provider block.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`pg_n_col`
       * :c:data:`pg_state_irrep_weight`
       * :c:data:`pg_thresh`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_state_irreps_label`

 
.. c:function:: void pg_small_det(void)


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        double precision function pg_small_det(a, n)


    det(a(1:n,1:n)), ``a`` untouched, exact zero returns exactly 0.

 
.. c:function:: void pg_state_apply(void)


    File : :file:`state_labels.irp.f`

    .. code:: fortran

        subroutine pg_state_apply(label, n, ok)


    Narrow ``ok(1:n)`` to the states whose known irrep is in :c:data:`irrep_states`.
    
    ``label(k)`` = 0 means "unknown": such a state is kept, because a filter
    must not drop a state it may not be subject to. It never widens, so it
    composes with any earlier selection (the S^2 test and the state
    following attribution mask in Davidson).
    
    :c:data:`pg_irrep_allowed` is provided here, not in the callers: IRPF90 hoists
    the explicit PROVIDE to the entry of this routine, so the mask and the
    detection chain are built at the first filtered call and by no
    caller's entry block.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_irrep_allowed`

 
.. c:function:: void pg_state_irreps_label(void)


    File : :file:`state_report.irp.f`

    .. code:: fortran

        subroutine pg_state_irreps_label(energy, s2)


    Forces the determinant labelling and gates the report on
    :c:data:`pg_label_ierr`. Separate from the dispatcher because IRPF90 hoists
    both the implicit needs and an explicit PROVIDE to the procedure
    entry: a body holding the PROVIDE would label before the guards.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_det_irrep`
       * :c:data:`pg_label_ierr`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_print_state_irreps`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pg_print_state_irreps_work`

 
.. c:function:: void pg_state_warn_filtered(void)


    File : :file:`state_labels.irp.f`

    .. code:: fortran

        subroutine pg_state_warn_filtered(label, n)


    One line per state whose known irrep the filter rejects:
    "symmetry: irrep filter: state <k> is <name>, not in irrep_states".
    This is the honest message when the requested irreps are not among the
    N_st kept states: raise n_states, or drop irrep_states.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`pg_irrep_allowed`
       * :c:data:`pg_n_col`

 
.. c:function:: void pg_unrank(void)


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        integer function pg_unrank(dim, n, rank)


    Bit mask of the ``rank``-th ``n``-subset of {0,...,dim-1}. ``dim`` <= 30
    (a symmetry block is an irrep space and has d_a <= 5 elements; the
    guard is in :c:func:`pg_elem_me`).

