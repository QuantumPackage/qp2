.. _module_iterations_tc: 
 
.. program:: iterations_tc 
 
.. default-role:: option 
 
iterations_tc
=============

No README.rst in this module; the providers and
subroutines below are documented from the irpf90 man pages. 
 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: get_extrapolated_energy:


    File : :file:`iterations_tc/iterations.irp.f`

    .. code:: fortran

        subroutine get_extrapolated_energy(Niter,ept2,pt1,extrap_energy)



    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`extrapolate_data`

 
.. c:function:: save_iterations:


    File : :file:`iterations_tc/iterations.irp.f`

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

