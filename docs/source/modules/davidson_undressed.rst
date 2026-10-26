.. _module_davidson_undressed: 
 
.. program:: davidson_undressed 
 
.. default-role:: option 
 
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
 
 
 
Providers 
--------- 
 
.. c:var:: double dressing_column_h


    File : :file:`davidson_undressed/null_dressing_vector.irp.f`

    .. code:: fortran

        double precision, allocatable	:: dressing_column_h	(N_det,N_states)
        double precision, allocatable	:: dressing_column_s	(N_det,N_states)
        double precision, allocatable	:: dressing_delta	(N_det,N_states)


    Null dressing vectors

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ci_electronic_energy`

 
.. c:var:: double dressing_column_s


    File : :file:`davidson_undressed/null_dressing_vector.irp.f`

    .. code:: fortran

        double precision, allocatable	:: dressing_column_h	(N_det,N_states)
        double precision, allocatable	:: dressing_column_s	(N_det,N_states)
        double precision, allocatable	:: dressing_delta	(N_det,N_states)


    Null dressing vectors

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ci_electronic_energy`

 
.. c:var:: double dressing_delta


    File : :file:`davidson_undressed/null_dressing_vector.irp.f`

    .. code:: fortran

        double precision, allocatable	:: dressing_column_h	(N_det,N_states)
        double precision, allocatable	:: dressing_column_s	(N_det,N_states)
        double precision, allocatable	:: dressing_delta	(N_det,N_states)


    Null dressing vectors

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ci_electronic_energy`

