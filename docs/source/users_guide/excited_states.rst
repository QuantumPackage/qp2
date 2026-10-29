.. _excited_states:

==============
Excited states
==============

Several electronic states are computed at once by setting
:option:`determinants n_states` to their number. The wave function then
has one column per state, the selection criterion is the largest
contribution over the states, and the program prints the energy of each
state together with the excitation energies relative to the first one.

.. code:: bash

    qp set determinants n_states 3
    qp run fci | tee hcn.fci.exc.out


Spin states
-----------

With :option:`determinants s2_eig` set to |true|, which is the default,
the wave functions are constrained to be eigenfunctions of |S^2| with the
value :option:`determinants expected_s2`. When this option is not set,
the expected value is :math:`S(S+1)` for the spin multiplicity of the
molecule, so a singlet calculation returns the lowest singlet states and
a doublet calculation the lowest doublet states. Asking for states of
different multiplicities in the same calculation requires
:option:`determinants s2_eig` set to |false|, and the wave functions may
then contain a mixture of spin states.

The ``only_expected_s2`` option, |true| by default, discards the
diagonalization vectors whose |S^2| deviates from ``expected_s2`` by more
than a hard-coded tolerance of the order of 0.5. Both filters are needed
in practice: the first one keeps the Hamiltonian block that is
diagonalized spin-pure, and the second one removes the vectors that the
iterative diagonalization produces from numerical noise.

.. code:: bash

    qp set determinants s2_eig true
    qp set determinants expected_s2 0.0


Guess wave functions
--------------------

Starting a multi-state |CIPSI| calculation from the Hartree-Fock
determinant works for the lowest state, but the excited states then need
many iterations to build their leading configurations. A |CIS| or
|CISD| calculation on the same |MO| space provides much better starting
vectors, and its cost is negligible for a small active space:

.. code:: bash

    qp set determinants n_states 3
    qp run cis
    qp set determinants read_wf true
    qp run fci | tee hcn.fci.exc.out

The :ref:`cis` program keeps the lowest :option:`determinants n_states`
states of the single-excitation space, and :ref:`cisd` does the same in
the space of single and double excitations. Both save their states in the
|EZFIO| database, and :option:`determinants read_wf` tells :ref:`fci` to
use them instead of the Hartree-Fock determinant. For HCN in 6-31G with
three states requested, the |CIS| guess gives

.. list-table::
   :header-rows: 1

   * - State
     - Energy (au)
     - Excitation (au)
     - Excitation (eV)
   * - 1
     - -92.82785666
     - 0.0
     - 0.0
   * - 2
     - -92.56236388
     - 0.26549278
     - 7.224
   * - 3
     - -92.54146594
     - 0.28639072
     - 7.793

The |CIPSI| run converges these energies and prints, for each state, a
block with the variational energy difference and the difference including
the |PT2| correction:

.. code:: text

    Variational Energy difference (au | eV)
    Delta E =   0.34252281292364728        9.3205241211619043
    Delta E =   0.34894534812463007        9.4952902739586094
    -----
    Variational + perturbative Energy difference (au | eV)
    Delta E =   0.34139629372938884        9.2898699605428909
    Delta E =   0.35037228589277447        9.5341192435500499

The excitation energies of a truncated calculation are not converged
yet: like the total energy, they must be checked against the
:option:`determinants n_det_max` and |PT2| criteria of :doc:`cipsi`.


Following the states
--------------------

The iterative diagonalization can swap two states from one iteration to
the next when their energies are close. The ``state_following`` option in
the :file:`davidson_keywords` section reorders the states so that each
one has the largest overlap with the state of the previous iteration:

.. code:: bash

    qp set davidson_keywords state_following true

It requires ``only_expected_s2`` to be |true| and at least as many
determinants as states, and the program stops with an explicit message
when these conditions are not met.

The ``n_states_diag`` option, also in :file:`davidson_keywords` and equal
to 4 by default, sets the number of vectors computed by the Davidson
algorithm to :option:`determinants n_states` times its value. Increasing
it gives the algorithm more room to converge the excited states when the
diagonalization is slow.


Natural orbitals
----------------

A state-averaged calculation benefits from |MOs| that describe all the
states equally. The :ref:`save_natorb` program computes the natural
orbitals of the wave function and replaces the |MOs| in the database with
them. For several states, the one-body density matrix is the average of
the density matrices of the states, weighted by
:option:`determinants state_average_weight`, which is a vector of
:option:`determinants n_states` numbers, equal by default:

.. code:: bash

    qp run save_natorb

The program invalidates the integrals stored in the database, so the next
program recomputes them in the new orbital basis. Running :ref:`scf` again
after :ref:`save_natorb` and :ref:`cis` is a common way to make sure that
the Hartree-Fock solution is the lowest one in the natural orbital basis.

|qp| does not use spatial symmetry, so states of different spatial
symmetries appear together in the same diagonalization. A state-averaged
calculation over states of different symmetries is therefore possible,
and the only separation that the program enforces is the spin one.

.. seealso::

    The documentation of the :ref:`scf`, :ref:`cis`, :ref:`cisd`,
    :ref:`save_natorb` and :ref:`fci` programs, and the
    :doc:`natural_orbitals` page.
