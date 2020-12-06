.. _module_iterations: 
 
.. program:: iterations 
 
.. default-role:: option 
 
==========
iterations
==========

Module which saves the computed energies for an extrapolation to
the |FCI| limit.
 
 
 
EZFIO parameters 
---------------- 
 
.. option:: n_iter
 
    Number of saved iterations
 
    Default: 1
 
.. option:: n_det_iterations
 
    Number of determinants at each iteration
 
 
.. option:: energy_iterations
 
    The variational energy at each iteration
 
 
.. option:: pt2_iterations
 
    The |PT2| correction at each iteration
 
 
 
Providers 
--------- 
 
.. c:var:: extrapolated_energy


    File : :file:`iterations/iterations.irp.f`

    .. code:: fortran

        double precision, allocatable	:: extrapolated_energy	(N_iter,N_states)


    Extrapolated energy, using E_var = f(PT2) where PT2=0

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`energy_iterations`
       * :c:data:`n_det`
       * :c:data:`n_iter`
       * :c:data:`n_states`
       * :c:data:`pt2_iterations`


 
.. c:var:: n_iter


    File : :file:`iterations/io.irp.f`

    .. code:: fortran

        integer	:: n_iter	


    number of iterations

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`mpi_master`
       * :c:data:`n_states`
       * :c:data:`output_wall_time_0`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`extrapolated_energy`

 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: print_extrapolated_energy:


    File : :file:`iterations/print_extrapolation.irp.f`

    .. code:: fortran

        subroutine print_extrapolated_energy


    Print the extrapolated energy in the output

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`extrapolated_energy`
       * :c:data:`n_det`
       * :c:data:`n_iter`
       * :c:data:`n_states`
       * :c:data:`pt2_iterations`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_cipsi`
       * :c:func:`run_stochastic_cipsi`

 
.. c:function:: print_summary:


    File : :file:`iterations/print_summary.irp.f`

    .. code:: fortran

        subroutine print_summary(e_,pt2_data,pt2_data_err,n_det_,n_occ_pattern_,n_st,s2_)


    Print the extrapolated energy in the output

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`do_pt2`
       * :c:data:`s2_eig`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_cipsi`
       * :c:func:`run_stochastic_cipsi`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`print_energy_components`

 
.. c:function:: save_iterations:


    File : :file:`iterations/iterations.irp.f`

    .. code:: fortran

        subroutine save_iterations(e_, pt2_,n_)


    Update the energy in the EZFIO file.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`energy_iterations`
       * :c:data:`n_det_iterations`
       * :c:data:`n_iter`
       * :c:data:`n_states`
       * :c:data:`pt2_iterations`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_cipsi`
       * :c:func:`run_stochastic_cipsi`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_set_iterations_energy_iterations`
       * :c:func:`ezfio_set_iterations_n_det_iterations`
       * :c:func:`ezfio_set_iterations_n_iter`
       * :c:func:`ezfio_set_iterations_pt2_iterations`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`n_iter`

