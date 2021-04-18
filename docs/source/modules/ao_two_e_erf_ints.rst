.. _module_ao_two_e_erf_ints: 
 
.. program:: ao_two_e_erf_ints 
 
.. default-role:: option 
 
======================
ao_two_e_erf_ints
======================

Here, all two-electron integrals (:math:`erf(\mu r_{12})/r_{12}`) are computed.
As they have 4 indices and many are zero, they are stored in a map, as defined
in :file:`utils/map_module.f90`.

The main parameter of this module is :option:`ao_two_e_erf_ints mu_erf` which is the range-separation parameter.

To fetch an |AO| integral, use the
`get_ao_two_e_integral_erf(i,j,k,l,ao_integrals_erf_map)` function.


The conventions are:
* For |AO| integrals : (ij|kl) = (11|22) = <ik|jl> = <12|12>



 
 
 
EZFIO parameters 
---------------- 
 
.. option:: io_ao_two_e_integrals_erf
 
    Read/Write |AO| integrals with the long range interaction from/to disk [ Write | Read | None ]
 
    Default: None
 
.. option:: mu_erf
 
    cutting of the interaction in the range separated model
 
    Default: 0.5
 
 
Providers 
--------- 
 
.. c:var:: ao_integrals_erf_cache


    File : :file:`ao_two_e_erf_ints/map_integrals_erf.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_integrals_erf_cache	(0:64*64*64*64)


    Cache of |AO| integrals for fast access

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_erf_cache_min`
       * :c:data:`ao_integrals_erf_map`
       * :c:data:`ao_two_e_integrals_erf_in_map`


 
.. c:var:: ao_integrals_erf_cache_max


    File : :file:`ao_two_e_erf_ints/map_integrals_erf.irp.f`

    .. code:: fortran

        integer	:: ao_integrals_erf_cache_min	
        integer	:: ao_integrals_erf_cache_max	


    Min and max values of the AOs for which the integrals are in the cache

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_erf_cache`

 
.. c:var:: ao_integrals_erf_cache_min


    File : :file:`ao_two_e_erf_ints/map_integrals_erf.irp.f`

    .. code:: fortran

        integer	:: ao_integrals_erf_cache_min	
        integer	:: ao_integrals_erf_cache_max	


    Min and max values of the AOs for which the integrals are in the cache

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_erf_cache`

 
.. c:var:: ao_integrals_erf_map


    File : :file:`ao_two_e_erf_ints/map_integrals_erf.irp.f`

    .. code:: fortran

        type(map_type)	:: ao_integrals_erf_map	


    |AO| integrals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_erf_cache`
       * :c:data:`ao_two_e_integrals_erf_in_map`
       * :c:data:`mo_two_e_int_erf_jj_from_ao`

 
.. c:var:: ao_two_e_integral_erf_schwartz


    File : :file:`ao_two_e_erf_ints/providers_ao_erf.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_two_e_integral_erf_schwartz	(ao_num,ao_num)


    Needed to compute Schwartz inequalities

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_expo_ordered_transp`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`mu_erf`
       * :c:data:`n_pt_max_integrals`
       * :c:data:`nucl_coord`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_two_e_int_erf_jj_from_ao`

 
.. c:var:: ao_two_e_integrals_erf_in_map


    File : :file:`ao_two_e_erf_ints/providers_ao_erf.irp.f`

    .. code:: fortran

        logical	:: ao_two_e_integrals_erf_in_map	


    Map of Atomic integrals
       i(r1) j(r2) 1/r12 k(r1) l(r2)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_expo_ordered_transp`
       * :c:data:`ao_integrals_erf_map`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`ezfio_filename`
       * :c:data:`io_ao_two_e_integrals_erf`
       * :c:data:`mu_erf`
       * :c:data:`n_pt_max_integrals`
       * :c:data:`nproc`
       * :c:data:`nucl_coord`
       * :c:data:`read_ao_two_e_integrals_erf`
       * :c:data:`zmq_context`
       * :c:data:`zmq_socket_pull_tcp_address`
       * :c:data:`zmq_state`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_erf_cache`
       * :c:data:`mo_two_e_int_erf_jj_from_ao`
       * :c:data:`mo_two_e_integrals_erf_in_map`

 
.. c:function:: general_primitive_integral_erf:


    File : :file:`ao_two_e_erf_ints/two_e_integrals_erf.irp.f`

    .. code:: fortran

        double precision function general_primitive_integral_erf(dim,            &
      P_new,P_center,fact_p,p,p_inv,iorder_p,                        &
      Q_new,Q_center,fact_q,q,q_inv,iorder_q)


    Computes the integral <pq|rs> where p,q,r,s are Gaussian primitives

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mu_erf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`add_poly_multiply`
       * :c:func:`give_polynom_mult_center_x`
       * :c:func:`multiply_poly`

 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: ao_two_e_integral_erf:


    File : :file:`ao_two_e_erf_ints/two_e_integrals_erf.irp.f`

    .. code:: fortran

        double precision function ao_two_e_integral_erf(i,j,k,l)


    integral of the AO basis <ik|jl> or (ij|kl)
       i(r1) j(r1) 1/r12 k(r2) l(r2)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_expo_ordered_transp`
       * :c:data:`ao_nucl`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`mu_erf`
       * :c:data:`n_pt_max_integrals`
       * :c:data:`nucl_coord`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_explicit_poly_and_gaussian`

 
.. c:function:: ao_two_e_integral_schwartz_accel_erf:


    File : :file:`ao_two_e_erf_ints/two_e_integrals_erf.irp.f`

    .. code:: fortran

        double precision function ao_two_e_integral_schwartz_accel_erf(i,j,k,l)


    integral of the AO basis <ik|jl> or (ij|kl)
       i(r1) j(r1) 1/r12 k(r2) l(r2)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_expo_ordered_transp`
       * :c:data:`ao_integrals_threshold`
       * :c:data:`ao_nucl`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`n_pt_max_integrals`
       * :c:data:`nucl_coord`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_explicit_poly_and_gaussian`

 
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

 
.. c:function:: clear_ao_erf_map:


    File : :file:`ao_two_e_erf_ints/map_integrals_erf.irp.f`

    .. code:: fortran

        subroutine clear_ao_erf_map


    Frees the memory of the |AO| map

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_erf_map`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`map_deinit`

 
.. c:function:: compute_ao_integrals_erf_jl:


    File : :file:`ao_two_e_erf_ints/two_e_integrals_erf.irp.f`

    .. code:: fortran

        subroutine compute_ao_integrals_erf_jl(j,l,n_integrals,buffer_i,buffer_value)


    Parallel client for AO integrals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_threshold`
       * :c:data:`ao_num`
       * :c:data:`ao_two_e_integral_erf_schwartz`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ao_two_e_integrals_erf_in_map_slave`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`two_e_integrals_index`

 
.. c:function:: compute_ao_two_e_integrals_erf:


    File : :file:`ao_two_e_erf_ints/two_e_integrals_erf.irp.f`

    .. code:: fortran

        subroutine compute_ao_two_e_integrals_erf(j,k,l,sze,buffer_value)


    Compute AO 1/r12 integrals for all i and fixed j,k,l

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_two_e_integral_erf_schwartz`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_two_e_int_erf_jj_from_ao`

 
.. c:function:: dump_ao_integrals_erf:


    File : :file:`ao_two_e_erf_ints/map_integrals_erf.irp.f`

    .. code:: fortran

        subroutine dump_ao_integrals_erf(filename)


    Save to disk the |AO| erf integrals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_erf_map`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_set_work_empty`

 
.. c:function:: eri_erf:


    File : :file:`ao_two_e_erf_ints/two_e_integrals_erf.irp.f`

    .. code:: fortran

        double precision function ERI_erf(alpha,beta,delta,gama,a_x,b_x,c_x,d_x,a_y,b_y,c_y,d_y,a_z,b_z,c_z,d_z)


    Atomic primtive two-electron integral between the 4 primitives :
    
    * primitive 1 : $x_1^{a_x} y_1^{a_y} z_1^{a_z} \exp(-\alpha * r1^2)$
    * primitive 2 : $x_1^{b_x} y_1^{b_y} z_1^{b_z} \exp(- \beta * r1^2)$
    * primitive 3 : $x_2^{c_x} y_2^{c_y} z_2^{c_z} \exp(-\delta * r2^2)$
    * primitive 4 : $x_2^{d_x} y_2^{d_y} z_2^{d_z} \exp(-\gamma * r2^2)$
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mu_erf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`integrale_new_erf`

 
.. c:function:: get_ao_erf_map_size:


    File : :file:`ao_two_e_erf_ints/map_integrals_erf.irp.f`

    .. code:: fortran

        function get_ao_erf_map_size()


    Returns the number of elements in the |AO| map

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_erf_map`

 
.. c:function:: get_ao_two_e_integral_erf:


    File : :file:`ao_two_e_erf_ints/map_integrals_erf.irp.f`

    .. code:: fortran

        double precision function get_ao_two_e_integral_erf(i,j,k,l,map) result(result)


    Gets one |AO| two-electron integral from the |AO| map

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_erf_cache`
       * :c:data:`ao_integrals_erf_cache_min`
       * :c:data:`ao_integrals_threshold`
       * :c:data:`ao_two_e_integral_erf_schwartz`
       * :c:data:`ao_two_e_integrals_erf_in_map`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`map_get`
       * :c:func:`two_e_integrals_index`

 
.. c:function:: get_ao_two_e_integrals_erf:


    File : :file:`ao_two_e_erf_ints/map_integrals_erf.irp.f`

    .. code:: fortran

        subroutine get_ao_two_e_integrals_erf(j,k,l,sze,out_val)


    Gets multiple |AO| two-electron integral from the |AO| map .
    All i are retrieved for j,k,l fixed.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_erf_map`
       * :c:data:`ao_integrals_threshold`
       * :c:data:`ao_two_e_integrals_erf_in_map`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`add_integrals_to_map_erf`

 
.. c:function:: get_ao_two_e_integrals_erf_non_zero:


    File : :file:`ao_two_e_erf_ints/map_integrals_erf.irp.f`

    .. code:: fortran

        subroutine get_ao_two_e_integrals_erf_non_zero(j,k,l,sze,out_val,out_val_index,non_zero_int)


    Gets multiple |AO| two-electron integrals from the |AO| map .
    All non-zero i are retrieved for j,k,l fixed.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_erf_map`
       * :c:data:`ao_integrals_threshold`
       * :c:data:`ao_two_e_integral_erf_schwartz`
       * :c:data:`ao_two_e_integrals_erf_in_map`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_two_e_int_erf_jj_from_ao`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`map_get`
       * :c:func:`two_e_integrals_index`

 
.. c:function:: insert_into_ao_integrals_erf_map:


    File : :file:`ao_two_e_erf_ints/map_integrals_erf.irp.f`

    .. code:: fortran

        subroutine insert_into_ao_integrals_erf_map(n_integrals,buffer_i, buffer_values)


    Create new entry into |AO| map

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_erf_map`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ao_two_e_integrals_erf_in_map_collector`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`map_append`

 
.. c:function:: integrale_new_erf:


    File : :file:`ao_two_e_erf_ints/two_e_integrals_erf.irp.f`

    .. code:: fortran

        subroutine integrale_new_erf(I_f,a_x,b_x,c_x,d_x,a_y,b_y,c_y,d_y,a_z,b_z,c_z,d_z,p,q,n_pt)


    Calculate the integral of the polynomial :
    
    $I_x1(a_x+b_x, c_x+d_x,p,q) \, I_x1(a_y+b_y, c_y+d_y,p,q) \, I_x1(a_z+b_z, c_z+d_z,p,q)$
    
    between $( 0 ; 1)$

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`gauleg_t2`
       * :c:data:`mu_erf`
       * :c:data:`n_pt_max_integrals`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`eri_erf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_new`

 
.. c:function:: load_ao_integrals_erf:


    File : :file:`ao_two_e_erf_ints/map_integrals_erf.irp.f`

    .. code:: fortran

        integer function load_ao_integrals_erf(filename)


    Read from disk the |AO| erf integrals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_erf_map`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`cache_map_reallocate`
       * :c:func:`map_deinit`
       * :c:func:`map_sort`

 
.. c:function:: save_erf_two_e_integrals_ao:


    File : :file:`ao_two_e_erf_ints/routines_save_integrals_erf.irp.f`

    .. code:: fortran

        subroutine save_erf_two_e_integrals_ao



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_erf_map`
       * :c:data:`ao_two_e_integrals_erf_in_map`
       * :c:data:`ezfio_filename`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`routine`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_set_ao_two_e_erf_ints_io_ao_two_e_integrals_erf`
       * :c:func:`ezfio_set_work_empty`
       * :c:func:`map_save_to_disk`

 
.. c:function:: save_erf_two_e_ints_ao_into_ints_ao:


    File : :file:`ao_two_e_erf_ints/routines_save_integrals_erf.irp.f`

    .. code:: fortran

        subroutine save_erf_two_e_ints_ao_into_ints_ao



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_erf_map`
       * :c:data:`ao_two_e_integrals_erf_in_map`
       * :c:data:`ezfio_filename`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_set_ao_two_e_ints_io_ao_two_e_integrals`
       * :c:func:`ezfio_set_work_empty`
       * :c:func:`map_save_to_disk`

