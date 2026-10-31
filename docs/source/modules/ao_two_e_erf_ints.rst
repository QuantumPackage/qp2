.. _module_ao_two_e_erf_ints: 
 
.. program:: ao_two_e_erf_ints 
 
.. default-role:: option 
 
ao_two_e_erf_ints
=================

No README.rst in this module; the providers and
subroutines below are documented from the irpf90 man pages. 
 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: ao_two_e_integrals_erf_in_map_collector:


    File : :file:`ao_two_e_erf_ints/integrals_erf_in_map_slave.irp.f`

    .. code:: fortran

        subroutine ao_two_e_integrals_erf_in_map_collector(zmq_socket_pull)


    Collects results from the AO integral calculation

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_erf_map`
       * :c:data:`ao_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_integrals_erf_in_map`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`end_zmq_to_qp_run_socket`
       * :c:func:`insert_into_ao_integrals_erf_map`

 
.. c:function:: ao_two_e_integrals_erf_in_map_slave:


    File : :file:`ao_two_e_erf_ints/integrals_erf_in_map_slave.irp.f`

    .. code:: fortran

        subroutine ao_two_e_integrals_erf_in_map_slave(thread,iproc)


    Computes a buffer of integrals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ao_two_e_integrals_erf_in_map_slave_inproc`
       * :c:func:`ao_two_e_integrals_erf_in_map_slave_tcp`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`compute_ao_integrals_erf_jl`
       * :c:func:`end_zmq_push_socket`
       * :c:func:`end_zmq_to_qp_run_socket`
       * :c:func:`push_integrals`

 
.. c:function:: ao_two_e_integrals_erf_in_map_slave_inproc:


    File : :file:`ao_two_e_erf_ints/integrals_erf_in_map_slave.irp.f`

    .. code:: fortran

        subroutine ao_two_e_integrals_erf_in_map_slave_inproc(i)


    Computes a buffer of integrals. i is the ID of the current thread.

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_integrals_erf_in_map`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ao_two_e_integrals_erf_in_map_slave`

 
.. c:function:: ao_two_e_integrals_erf_in_map_slave_tcp:


    File : :file:`ao_two_e_erf_ints/integrals_erf_in_map_slave.irp.f`

    .. code:: fortran

        subroutine ao_two_e_integrals_erf_in_map_slave_tcp(i)


    Computes a buffer of integrals. i is the ID of the current thread.

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ao_two_e_integrals_erf_in_map_slave`

