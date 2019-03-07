.. _excited_states:

==============
Excited states
==============

It is possible to run excited states calculations with the quantum
package. To do this, set :option:`determinants n_states` to the number
of requested states. The selection criterion will be the maximum of the
selection criteria for each state. If the Davidson diagonalization has
difficulties to converge, increase the :option:`davidson n_states_diag`
value.

When computing multiple states, it is good to have the
:option:`determinants s2_eig` flag |true|. This will force the Davidson
algorithm to choose only vectors with a value of |S^2| equal to
:option:`determinants expected_s2`. Otherwise, different spin states
will come out in the diagonalization.

The |qp| doesn't take account of the symmetry. Due to numerical noise,
excited states of different symmetries may enter in the calculation.
Note that it is possible to make state-average calculation of states
with different symmetries and/or different spin multiplicities.

To include excited states of all possible symmetries, a simple trick is
to run a preliminary multi-state |CIS| calculation using the :ref:`CIS`
program, and then running the selected |FCI| restarting from the |CIS|
states, setting :option:`determinants read_wf` to |true|.

Usually, it is good practice to use state-averaged natural |MOs| so that
all states have |MOs| of comparable quality. This allows for a faster
convergence of excitation energies.


.. seealso:: 

    The documentation of the :ref:`scf`, :ref:`cis` and
    :ref:`fci` programs.

