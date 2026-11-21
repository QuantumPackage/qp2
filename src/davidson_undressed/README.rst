==================
davidson_undressed
==================

Davidson diagonalization with no dressing.

For the programs that diagonalize the undressed Hamiltonian, this
module provides the dressing machinery with null dressing vectors:
:c:data:`dressing_column_h`, :c:data:`dressing_column_s` and
:c:data:`dressing_delta` are all zero, so that the Davidson algorithm
of the :ref:`module_davidson_dressed` module applies to the bare
Hamiltonian.
