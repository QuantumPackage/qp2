====
mrci
====

|CIPSI| algorithm in the multi-reference |CI| space (|CAS| + singles
and doubles).

The :ref:`mrci` program performs a |CIPSI|-like selected |CI| using
a stochastic scheme for both the selection of the important Slater
determinants and the computation of the |PT2| correction, for the
:option:`determinants n_states` lowest states of the variational
space. The inactive, active and virtual |MOs| are set with the
:ref:`qp_set_mo_class` program.

The boolean options `do_ddci` and `do_only_1h1p` restrict the space
of excitations: removing the purely inactive doubles, or keeping only
the one hole/one particle excitations.

.. seealso::

    The documentation of the :ref:`module_cipsi` module.
