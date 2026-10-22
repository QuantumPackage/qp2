.. _module_iterations: 
 
.. program:: iterations 
 
.. default-role:: option 
 
==========
iterations
==========

Iteration counters and energy history, for the extrapolation to the
|FCI| limit.

:c:data:`n_iter` and :c:data:`n_iter_max` count the iterations
(:c:func:`increment_n_iter`), and the energies of the iterations are
stored in :c:data:`energy_iterations` (the variational energies) and
:c:data:`pt2_iterations` (the |PT2| energies).
:c:data:`extrapolated_energy` is the value extrapolated to the |FCI|
limit, computed and printed by :c:func:`print_extrapolated_energy`,
and the summary of the calculation is printed by
:c:func:`print_summary` (or :c:func:`print_summary_tc` for the
transcorrelated variant).
 
 
 
Providers 
--------- 
 
.. c:var:: energy_iterations


    File : :file:`iterations/iterations.irp.f`

    .. code:: fortran

        double precision, allocatable	:: energy_iterations	(n_states,N_iter_max)
        double precision, allocatable	:: pt2_iterations	(n_states,N_iter_max)
        double precision, allocatable	:: extrapolated_energy	(N_iter_max,N_states)


    The energy at each iteration for the extrapolations

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_iter_max`
       * :c:data:`n_states`


 
.. c:var:: extrapolated_energy


    File : :file:`iterations/iterations.irp.f`

    .. code:: fortran

        double precision, allocatable	:: energy_iterations	(n_states,N_iter_max)
        double precision, allocatable	:: pt2_iterations	(n_states,N_iter_max)
        double precision, allocatable	:: extrapolated_energy	(N_iter_max,N_states)


    The energy at each iteration for the extrapolations

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_iter_max`
       * :c:data:`n_states`


 
.. c:var:: n_iter


    File : :file:`iterations/iterations.irp.f`

    .. code:: fortran

        integer	:: n_iter	


    Number of CIPSI iterations


 
.. c:var:: n_iter_max


    File : :file:`iterations/iterations.irp.f`

    .. code:: fortran

        integer	:: n_iter_max	


    Max number of iterations for extrapolations

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`energy_iterations`

 
.. c:var:: pt2_iterations


    File : :file:`iterations/iterations.irp.f`

    .. code:: fortran

        double precision, allocatable	:: energy_iterations	(n_states,N_iter_max)
        double precision, allocatable	:: pt2_iterations	(n_states,N_iter_max)
        double precision, allocatable	:: extrapolated_energy	(N_iter_max,N_states)


    The energy at each iteration for the extrapolations

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_iter_max`
       * :c:data:`n_states`


 
 
Subroutines / functions 
----------------------- 
 
.. c:macro:: increment_n_iter:


    File : :file:`iterations/iterations.irp.f`

    .. code:: fortran

        subroutine increment_n_iter(e, pt2_data)


    Does what is necessary to increment n_iter

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`energy_iterations`
       * :c:data:`n_det`
       * :c:data:`n_iter`
       * :c:data:`n_iter_max`
       * :c:data:`n_states`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_cipsi`
       * :c:func:`run_stochastic_cipsi`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`extrapolate_data`

 
.. c:macro:: print_extrapolated_energy:


    File : :file:`iterations/print_extrapolation.irp.f`

    .. code:: fortran

        subroutine print_extrapolated_energy


    Print the extrapolated energy in the output

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`energy_iterations`
       * :c:data:`n_det`
       * :c:data:`n_iter`
       * :c:data:`n_states`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_cipsi`
       * :c:func:`run_stochastic_cipsi`

 
.. c:macro:: print_summary:


    File : :file:`iterations/print_summary.irp.f`

    .. code:: fortran

        subroutine print_summary(e_,pt2_data,pt2_data_err,n_det_,n_configuration_,n_st,s2_)


    Print the extrapolated energy in the output

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`csf_based`
       * :c:data:`do_pt2`
       * :c:data:`n_csf`
       * :c:data:`s2_eig`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_cipsi`
       * :c:func:`run_stochastic_cipsi`

 
.. c:macro:: print_summary_tc:


    File : :file:`iterations/summary_tc.irp.f`

    .. code:: fortran

        subroutine print_summary_tc(e_,pt2_data,pt2_data_err,n_det_,n_configuration_,n_st,s2_)


    Print the extrapolated energy in the output

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`do_pt2`
       * :c:data:`n_csf`
       * :c:data:`only_expected_s2`
       * :c:data:`s2_eig`

