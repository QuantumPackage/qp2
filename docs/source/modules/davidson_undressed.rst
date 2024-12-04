.. _module_davidson_undressed: 
 
.. program:: davidson_undressed 
 
.. default-role:: option 
 
==================
davidson_undressed
==================

Module for main files Davidson's algorithm with no dressing.

 
 
 
Providers 
--------- 
 
.. c:var:: dressing_column_h


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

 
.. c:var:: dressing_column_s


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

 
.. c:var:: dressing_delta


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

