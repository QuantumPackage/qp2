.. _module_ao_two_e_ints: 
 
.. program:: ao_two_e_ints 
 
.. default-role:: option 
 
==================
ao_two_e_ints
==================

Here, all two-electron integrals (:math:`1/r_{12}`) are computed.
As they have 4 indices and many are zero, they are stored in a map, as defined
in :file:`utils/map_module.f90`.

To fetch an |AO| integral, use the
`get_ao_two_e_integral(i,j,k,l,ao_integrals_map)` function.


The conventions are:
* For |AO| integrals : (ij|kl) = (11|22) = <ik|jl> = <12|12>



 
 
 
EZFIO parameters 
---------------- 
 
.. option:: io_ao_two_e_integrals
 
    Read/Write |AO| integrals from/to disk [ Write | Read | None ]
 
    Default: None
 
.. option:: io_ao_cholesky
 
    Read/Write |AO| Cholesky integrals from/to disk [ Write | Read | None ]
 
    Default: None
 
.. option:: ao_integrals_threshold
 
    If | (pq|rs) | < `ao_integrals_threshold` then (pq|rs) is zero
 
    Default: 1.e-15
 
.. option:: ao_cholesky_threshold
 
    If | (ii|jj) | < `ao_cholesky_threshold` then (ii|jj) is zero
 
    Default: 1.e-12
 
.. option:: do_ao_cholesky
 
    Perform Cholesky decomposition of AO integrals
 
    Default: True
 
.. option:: io_ao_two_e_integrals_erf
 
    Read/Write |AO| erf integrals from/to disk [ Write | Read | None ]
 
    Default: None
 
.. option:: use_only_lr
 
    If true, use only the long range part of the two-electron integrals instead of 1/r12
 
    Default: False
 
 
Providers 
--------- 
 
.. c:var:: ao_2e_cgtos_schwartz


    File : :file:`ao_two_e_ints/two_e_coul_integrals_cgtos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_2e_cgtos_schwartz	(ao_num,ao_num)


    Needed to compute Schwartz inequalities

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_cgtos_norm_ord_transp`
       * :c:data:`ao_expo_cgtos_ord_transp`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`n_pt_max_integrals`
       * :c:data:`nucl_coord`
       * :c:data:`use_pw`


 
.. c:var:: ao_2e_cosgtos_schwartz


    File : :file:`ao_two_e_ints/two_e_Coul_integrals_cosgtos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_2e_cosgtos_schwartz	(ao_num,ao_num)


    Needed to compute Schwartz inequalities

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_norm_ord_transp_cosgtos`
       * :c:data:`ao_expo_ord_transp_cosgtos`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`n_pt_max_integrals`
       * :c:data:`nucl_coord`


 
.. c:var:: ao_integrals_cache


    File : :file:`ao_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_integrals_cache	(0:64*64*64*64)


    Cache of AO integrals for fast access

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_cache_min`
       * :c:data:`ao_integrals_map`
       * :c:data:`ao_two_e_integrals_in_map`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cholesky_ao_num`

 
.. c:var:: ao_integrals_cache_max


    File : :file:`ao_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        integer	:: ao_integrals_cache_min	
        integer	:: ao_integrals_cache_max	


    Min and max values of the AOs for which the integrals are in the cache

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_cache`
       * :c:data:`ao_integrals_cache_periodic`
       * :c:data:`cholesky_ao_num`

 
.. c:var:: ao_integrals_cache_min


    File : :file:`ao_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        integer	:: ao_integrals_cache_min	
        integer	:: ao_integrals_cache_max	


    Min and max values of the AOs for which the integrals are in the cache

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_cache`
       * :c:data:`ao_integrals_cache_periodic`
       * :c:data:`cholesky_ao_num`

 
.. c:var:: ao_integrals_cache_periodic


    File : :file:`ao_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        complex*16, allocatable	:: ao_integrals_cache_periodic	(0:64*64*64*64)


    Cache of AO integrals for fast access

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_cache_min`
       * :c:data:`ao_integrals_map`
       * :c:data:`ao_two_e_integrals_in_map`


 
.. c:var:: ao_integrals_erf_cache


    File : :file:`ao_two_e_ints/map_integrals_erf.irp.f`

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


    File : :file:`ao_two_e_ints/map_integrals_erf.irp.f`

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


    File : :file:`ao_two_e_ints/map_integrals_erf.irp.f`

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


    File : :file:`ao_two_e_ints/map_integrals_erf.irp.f`

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

 
.. c:var:: ao_integrals_map


    File : :file:`ao_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        type(map_type)	:: ao_integrals_map	


    AO integrals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_cache`
       * :c:data:`ao_integrals_cache_periodic`
       * :c:data:`ao_two_e_integral_alpha`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`cholesky_ao_num`

 
.. c:var:: ao_two_e_integral_erf_schwartz


    File : :file:`ao_two_e_ints/providers_ao_erf.irp.f`

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

 
.. c:var:: ao_two_e_integral_schwartz


    File : :file:`ao_two_e_ints/two_e_integrals.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_two_e_integral_schwartz	(ao_num,ao_num)


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
       * :c:data:`n_pt_max_integrals`
       * :c:data:`nucl_coord`
       * :c:data:`use_cgtos`
       * :c:data:`use_only_lr`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_integral_alpha`
       * :c:data:`cholesky_ao_num`

 
.. c:var:: ao_two_e_integrals_erf_in_map


    File : :file:`ao_two_e_ints/providers_ao_erf.irp.f`

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

 
.. c:var:: ao_two_e_integrals_in_map


    File : :file:`ao_two_e_ints/two_e_integrals.irp.f`

    .. code:: fortran

        logical	:: ao_two_e_integrals_in_map	


    Map of Atomic integrals
       i(r1) j(r2) 1/r12 k(r1) l(r2)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_expo_ordered_transp`
       * :c:data:`ao_integrals_map`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`ezfio_filename`
       * :c:data:`io_ao_two_e_integrals`
       * :c:data:`mpi_master`
       * :c:data:`n_pt_max_integrals`
       * :c:data:`nproc`
       * :c:data:`nucl_coord`
       * :c:data:`read_ao_two_e_integrals`
       * :c:data:`use_cgtos`
       * :c:data:`use_only_lr`
       * :c:data:`zmq_context`
       * :c:data:`zmq_socket_pull_tcp_address`
       * :c:data:`zmq_state`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_cache`
       * :c:data:`ao_integrals_cache_periodic`
       * :c:data:`ao_two_e_integral_alpha`
       * :c:data:`cholesky_ao_num`
       * :c:data:`mo_two_e_integrals_erf_in_map`
       * :c:data:`mo_two_e_integrals_in_map`

 
.. c:var:: cholesky_ao


    File : :file:`ao_two_e_ints/cholesky.irp.f`

    .. code:: fortran

        integer	:: cholesky_ao_num	
        double precision, allocatable	:: cholesky_ao	(ao_num,ao_num,1)


    Cholesky vectors in AO basis: (ik|a):
    <ij|kl> = (ik|jl) = sum_a (ik|a).(a|jl)
    
    Last dimension of cholesky_ao is cholesky_ao_num
    
    https://mogp-emulator.readthedocs.io/en/latest/methods/proc/ProcPivotedCholesky.html
    
    https://doi.org/10.1016/j.apnum.2011.10.001 : Page 4, Algorithm 1
    
    https://www.diva-portal.org/smash/get/diva2:396223/FULLTEXT01.pdf

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_cholesky_threshold`
       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_expo_ordered_transp`
       * :c:data:`ao_integrals_cache`
       * :c:data:`ao_integrals_cache_min`
       * :c:data:`ao_integrals_map`
       * :c:data:`ao_integrals_threshold`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_overlap_abs`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`ao_two_e_integral_schwartz`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`do_direct_integrals`
       * :c:data:`elec_num`
       * :c:data:`ezfio_work_dir`
       * :c:data:`is_periodic`
       * :c:data:`n_pt_max_integrals`
       * :c:data:`nproc`
       * :c:data:`nucl_coord`
       * :c:data:`qp_max_mem`
       * :c:data:`read_ao_cholesky`
       * :c:data:`read_ao_two_e_integrals`
       * :c:data:`use_cgtos`
       * :c:data:`use_only_lr`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_integral_alpha_chol`
       * :c:data:`cholesky_ao_transp`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`

 
.. c:var:: cholesky_ao_num


    File : :file:`ao_two_e_ints/cholesky.irp.f`

    .. code:: fortran

        integer	:: cholesky_ao_num	
        double precision, allocatable	:: cholesky_ao	(ao_num,ao_num,1)


    Cholesky vectors in AO basis: (ik|a):
    <ij|kl> = (ik|jl) = sum_a (ik|a).(a|jl)
    
    Last dimension of cholesky_ao is cholesky_ao_num
    
    https://mogp-emulator.readthedocs.io/en/latest/methods/proc/ProcPivotedCholesky.html
    
    https://doi.org/10.1016/j.apnum.2011.10.001 : Page 4, Algorithm 1
    
    https://www.diva-portal.org/smash/get/diva2:396223/FULLTEXT01.pdf

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_cholesky_threshold`
       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_expo_ordered_transp`
       * :c:data:`ao_integrals_cache`
       * :c:data:`ao_integrals_cache_min`
       * :c:data:`ao_integrals_map`
       * :c:data:`ao_integrals_threshold`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_overlap_abs`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`ao_two_e_integral_schwartz`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`do_direct_integrals`
       * :c:data:`elec_num`
       * :c:data:`ezfio_work_dir`
       * :c:data:`is_periodic`
       * :c:data:`n_pt_max_integrals`
       * :c:data:`nproc`
       * :c:data:`nucl_coord`
       * :c:data:`qp_max_mem`
       * :c:data:`read_ao_cholesky`
       * :c:data:`read_ao_two_e_integrals`
       * :c:data:`use_cgtos`
       * :c:data:`use_only_lr`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_integral_alpha_chol`
       * :c:data:`cholesky_ao_transp`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`

 
.. c:var:: cholesky_ao_transp


    File : :file:`ao_two_e_ints/cholesky.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cholesky_ao_transp	(cholesky_ao_num,ao_num,ao_num)


    Transposed of the Cholesky vectors in AO basis set

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`cholesky_ao_num`


 
.. c:var:: do_direct_integrals


    File : :file:`ao_two_e_ints/direct.irp.f`

    .. code:: fortran

        logical	:: do_direct_integrals	


    Compute integrals on the fly

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`do_ao_cholesky`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_integral_alpha`
       * :c:data:`cholesky_ao_num`
       * :c:data:`mo_two_e_int_erf_jj_from_ao`

 
.. c:var:: gauleg_t2


    File : :file:`ao_two_e_ints/gauss_legendre.irp.f`

    .. code:: fortran

        double precision, allocatable	:: gauleg_t2	(n_pt_max_integrals,n_pt_max_integrals/2)
        double precision, allocatable	:: gauleg_w	(n_pt_max_integrals,n_pt_max_integrals/2)


    t_w(i,1,k) = w(i)
    t_w(i,2,k) = t(i)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_pt_max_integrals`


 
.. c:var:: gauleg_w


    File : :file:`ao_two_e_ints/gauss_legendre.irp.f`

    .. code:: fortran

        double precision, allocatable	:: gauleg_t2	(n_pt_max_integrals,n_pt_max_integrals/2)
        double precision, allocatable	:: gauleg_w	(n_pt_max_integrals,n_pt_max_integrals/2)


    t_w(i,1,k) = w(i)
    t_w(i,2,k) = t(i)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_pt_max_integrals`


 
.. c:function:: general_primitive_integral:


    File : :file:`ao_two_e_ints/two_e_integrals.irp.f`

    .. code:: fortran

        double precision function general_primitive_integral(dim,            &
      P_new,P_center,fact_p,p,p_inv,iorder_p,                        &
      Q_new,Q_center,fact_q,q,q_inv,iorder_q)


    Computes the integral <pq|rs> where p,q,r,s are Gaussian primitives

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`add_poly_multiply`
       * :c:func:`give_polynom_mult_center_x`

 
.. c:function:: general_primitive_integral_cgtos:


    File : :file:`ao_two_e_ints/two_e_coul_integrals_cgtos.irp.f`

    .. code:: fortran

        complex*16 function general_primitive_integral_cgtos(dim, P_new, P_center, fact_p, p, p_inv, iorder_p, &
                                                          Q_new, Q_center, fact_q, q, q_inv, iorder_q)


    
    Computes the integral <pq|rs> where p,q,r,s are cos-cGTOS primitives
    

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`add_cpoly_multiply`
       * :c:func:`give_cpolynom_mult_center_x`
       * :c:func:`multiply_cpoly`

 
.. c:function:: general_primitive_integral_cosgtos:


    File : :file:`ao_two_e_ints/two_e_Coul_integrals_cosgtos.irp.f`

    .. code:: fortran

        complex*16 function general_primitive_integral_cosgtos( dim, P_new, P_center, fact_p, p, p_inv, iorder_p &
                                                           , Q_new, Q_center, fact_q, q, q_inv, iorder_q )


    
    Computes the integral <pq|rs> where p,q,r,s are cos-cGTOS primitives
    

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`add_cpoly_multiply`
       * :c:func:`give_cpolynom_mult_center_x`
       * :c:func:`multiply_cpoly`

 
.. c:function:: general_primitive_integral_erf:


    File : :file:`ao_two_e_ints/two_e_integrals_erf.irp.f`

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

 
.. c:function:: give_cpolynom_mult_center_x:


    File : :file:`ao_two_e_ints/two_e_coul_integrals_cgtos.irp.f`

    .. code:: fortran

        subroutine give_cpolynom_mult_center_x(P_center, Q_center, a_x, d_x, p, q, n_pt_in, &
                                       pq_inv, pq_inv_2, p10_1, p01_1, p10_2, p01_2, d, n_pt_out)


    subroutine that returns the explicit polynom in term of the "t"
    variable of the following polynoms :
    
    $I_{x_1}(a_x,d_x,p,q) \, I_{x_1}(a_y,d_y,p,q) \ I_{x_1}(a_z,d_z,p,q)$

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`general_primitive_integral_cgtos`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_pol_mult_cgtos`

 
.. c:function:: i_x1_new:


    File : :file:`ao_two_e_ints/two_e_integrals.irp.f`

    .. code:: fortran

        recursive subroutine I_x1_new(a,c,B_10,B_01,B_00,res,n_pt)


    recursive function involved in the two-electron integral

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_pt_max_integrals`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_new`
       * :c:func:`i_x2_new`
       * :c:func:`integrale_new`
       * :c:func:`integrale_new_erf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_new`
       * :c:func:`i_x2_new`

 
.. c:function:: i_x1_new_cgtos:


    File : :file:`ao_two_e_ints/two_e_coul_integrals_cgtos.irp.f`

    .. code:: fortran

        recursive subroutine I_x1_new_cgtos(a, c, B_10, B_01, B_00, res, n_pt)


    recursive function involved in the two-electron integral

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_pt_max_integrals`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_new_cgtos`
       * :c:func:`i_x2_new_cgtos`
       * :c:func:`integrale_new_cgtos`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_new_cgtos`
       * :c:func:`i_x2_new_cgtos`

 
.. c:function:: i_x1_new_cosgtos:


    File : :file:`ao_two_e_ints/two_e_Coul_integrals_cosgtos.irp.f`

    .. code:: fortran

        recursive subroutine I_x1_new_cosgtos(a, c, B_10, B_01, B_00, res, n_pt)


    recursive function involved in the two-electron integral

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_pt_max_integrals`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_new_cosgtos`
       * :c:func:`i_x2_new_cosgtos`
       * :c:func:`integrale_new_cosgtos`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_new_cosgtos`
       * :c:func:`i_x2_new_cosgtos`

 
.. c:function:: i_x1_pol_mult_a1:


    File : :file:`ao_two_e_ints/two_e_integrals.irp.f`

    .. code:: fortran

        recursive subroutine I_x1_pol_mult_a1(c,B_10,B_01,B_00,C_00,D_00,d,nd,n_pt_in)


    Recursive function involved in the two-electron integral

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_pol_mult`
       * :c:func:`i_x1_pol_mult_a2`
       * :c:func:`i_x1_pol_mult_recurs`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_x2_pol_mult`

 
.. c:function:: i_x1_pol_mult_a1_cgtos:


    File : :file:`ao_two_e_ints/two_e_coul_integrals_cgtos.irp.f`

    .. code:: fortran

        recursive subroutine I_x1_pol_mult_a1_cgtos(c,B_10,B_01,B_00,C_00,D_00,d,nd,n_pt_in)


    Recursive function involved in the two-electron integral

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_pol_mult_a2_cgtos`
       * :c:func:`i_x1_pol_mult_cgtos`
       * :c:func:`i_x1_pol_mult_recurs_cgtos`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_x2_pol_mult_cgtos`
       * :c:func:`multiply_cpoly`

 
.. c:function:: i_x1_pol_mult_a1_cosgtos:


    File : :file:`ao_two_e_ints/two_e_Coul_integrals_cosgtos.irp.f`

    .. code:: fortran

        recursive subroutine I_x1_pol_mult_a1_cosgtos(c,B_10,B_01,B_00,C_00,D_00,d,nd,n_pt_in)


    Recursive function involved in the two-electron integral

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_pol_mult_a2_cosgtos`
       * :c:func:`i_x1_pol_mult_cosgtos`
       * :c:func:`i_x1_pol_mult_recurs_cosgtos`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_x2_pol_mult_cosgtos`
       * :c:func:`multiply_cpoly`

 
.. c:function:: i_x1_pol_mult_a2:


    File : :file:`ao_two_e_ints/two_e_integrals.irp.f`

    .. code:: fortran

        recursive subroutine I_x1_pol_mult_a2(c,B_10,B_01,B_00,C_00,D_00,d,nd,n_pt_in)


    Recursive function involved in the two-electron integral

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_pol_mult`
       * :c:func:`i_x1_pol_mult_recurs`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_pol_mult_a1`
       * :c:func:`i_x2_pol_mult`

 
.. c:function:: i_x1_pol_mult_a2_cgtos:


    File : :file:`ao_two_e_ints/two_e_coul_integrals_cgtos.irp.f`

    .. code:: fortran

        recursive subroutine I_x1_pol_mult_a2_cgtos(c, B_10, B_01, B_00, C_00, D_00, d, nd, n_pt_in)


    Recursive function involved in the two-electron integral

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_pol_mult_cgtos`
       * :c:func:`i_x1_pol_mult_recurs_cgtos`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_pol_mult_a1_cgtos`
       * :c:func:`i_x2_pol_mult_cgtos`
       * :c:func:`multiply_cpoly`

 
.. c:function:: i_x1_pol_mult_a2_cosgtos:


    File : :file:`ao_two_e_ints/two_e_Coul_integrals_cosgtos.irp.f`

    .. code:: fortran

        recursive subroutine I_x1_pol_mult_a2_cosgtos(c, B_10, B_01, B_00, C_00, D_00, d, nd, n_pt_in)


    Recursive function involved in the two-electron integral

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_pol_mult_cosgtos`
       * :c:func:`i_x1_pol_mult_recurs_cosgtos`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_pol_mult_a1_cosgtos`
       * :c:func:`i_x2_pol_mult_cosgtos`
       * :c:func:`multiply_cpoly`

 
.. c:function:: i_x1_pol_mult_recurs:


    File : :file:`ao_two_e_ints/two_e_integrals.irp.f`

    .. code:: fortran

        recursive subroutine I_x1_pol_mult_recurs(a,c,B_10,B_01,B_00,C_00,D_00,d,nd,n_pt_in)


    Recursive function involved in the two-electron integral

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_pol_mult`
       * :c:func:`i_x1_pol_mult_recurs`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_pol_mult_a1`
       * :c:func:`i_x1_pol_mult_a2`
       * :c:func:`i_x1_pol_mult_recurs`

 
.. c:function:: i_x1_pol_mult_recurs_cgtos:


    File : :file:`ao_two_e_ints/two_e_coul_integrals_cgtos.irp.f`

    .. code:: fortran

        recursive subroutine I_x1_pol_mult_recurs_cgtos(a, c, B_10, B_01, B_00, C_00, D_00, d, nd, n_pt_in)


    Recursive function involved in the two-electron integral

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_pol_mult_cgtos`
       * :c:func:`i_x1_pol_mult_recurs_cgtos`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_pol_mult_a1_cgtos`
       * :c:func:`i_x1_pol_mult_a2_cgtos`
       * :c:func:`i_x1_pol_mult_recurs_cgtos`
       * :c:func:`multiply_cpoly`

 
.. c:function:: i_x1_pol_mult_recurs_cosgtos:


    File : :file:`ao_two_e_ints/two_e_Coul_integrals_cosgtos.irp.f`

    .. code:: fortran

        recursive subroutine I_x1_pol_mult_recurs_cosgtos(a, c, B_10, B_01, B_00, C_00, D_00, d, nd, n_pt_in)


    Recursive function involved in the two-electron integral

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_pol_mult_cosgtos`
       * :c:func:`i_x1_pol_mult_recurs_cosgtos`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_pol_mult_a1_cosgtos`
       * :c:func:`i_x1_pol_mult_a2_cosgtos`
       * :c:func:`i_x1_pol_mult_recurs_cosgtos`
       * :c:func:`multiply_cpoly`

 
.. c:function:: i_x2_new:


    File : :file:`ao_two_e_ints/two_e_integrals.irp.f`

    .. code:: fortran

        recursive subroutine I_x2_new(c,B_10,B_01,B_00,res,n_pt)


    recursive function involved in the two-electron integral

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_pt_max_integrals`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_new`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_new`

 
.. c:function:: i_x2_new_cgtos:


    File : :file:`ao_two_e_ints/two_e_coul_integrals_cgtos.irp.f`

    .. code:: fortran

        recursive subroutine I_x2_new_cgtos(c, B_10, B_01, B_00, res, n_pt)


    recursive function involved in the two-electron integral

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_pt_max_integrals`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_new_cgtos`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_new_cgtos`

 
.. c:function:: i_x2_new_cosgtos:


    File : :file:`ao_two_e_ints/two_e_Coul_integrals_cosgtos.irp.f`

    .. code:: fortran

        recursive subroutine I_x2_new_cosgtos(c, B_10, B_01, B_00, res, n_pt)


    recursive function involved in the two-electron integral

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_pt_max_integrals`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_new_cosgtos`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_new_cosgtos`

 
.. c:function:: i_x2_pol_mult:


    File : :file:`ao_two_e_ints/two_e_integrals.irp.f`

    .. code:: fortran

        recursive subroutine I_x2_pol_mult(c,B_10,B_01,B_00,C_00,D_00,d,nd,dim)


    Recursive function involved in the two-electron integral

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_pol_mult`
       * :c:func:`i_x1_pol_mult_a1`
       * :c:func:`i_x1_pol_mult_a2`
       * :c:func:`i_x2_pol_mult`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_x2_pol_mult`

 
.. c:function:: i_x2_pol_mult_cgtos:


    File : :file:`ao_two_e_ints/two_e_coul_integrals_cgtos.irp.f`

    .. code:: fortran

        recursive subroutine I_x2_pol_mult_cgtos(c, B_10, B_01, B_00, C_00, D_00, d, nd, dim)


    Recursive function involved in the two-electron integral

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_pol_mult_a1_cgtos`
       * :c:func:`i_x1_pol_mult_a2_cgtos`
       * :c:func:`i_x1_pol_mult_cgtos`
       * :c:func:`i_x2_pol_mult_cgtos`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_x2_pol_mult_cgtos`
       * :c:func:`multiply_cpoly`

 
.. c:function:: i_x2_pol_mult_cosgtos:


    File : :file:`ao_two_e_ints/two_e_Coul_integrals_cosgtos.irp.f`

    .. code:: fortran

        recursive subroutine I_x2_pol_mult_cosgtos(c, B_10, B_01, B_00, C_00, D_00, d, nd, dim)


    Recursive function involved in the two-electron integral

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_pol_mult_a1_cosgtos`
       * :c:func:`i_x1_pol_mult_a2_cosgtos`
       * :c:func:`i_x1_pol_mult_cosgtos`
       * :c:func:`i_x2_pol_mult_cosgtos`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_x2_pol_mult_cosgtos`
       * :c:func:`multiply_cpoly`

 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: ao_2e_cgtos_schwartz_accel:


    File : :file:`ao_two_e_ints/two_e_coul_integrals_cgtos.irp.f`

    .. code:: fortran

        double precision function ao_2e_cgtos_schwartz_accel(i, j, k, l)


    integral of the AO basis <ik|jl> or (ij|kl)
       i(r1) j(r1) 1/r12 k(r2) l(r2)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_cgtos_norm_ord_transp`
       * :c:data:`ao_expo_cgtos_ord_transp`
       * :c:data:`ao_integrals_threshold`
       * :c:data:`ao_nucl`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`n_pt_max_integrals`
       * :c:data:`nucl_coord`
       * :c:data:`use_pw`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_explicit_cpoly_and_cgaussian`

 
.. c:function:: ao_2e_cosgtos_schwartz_accel:


    File : :file:`ao_two_e_ints/two_e_Coul_integrals_cosgtos.irp.f`

    .. code:: fortran

        double precision function ao_2e_cosgtos_schwartz_accel(i, j, k, l)


    integral of the AO basis <ik|jl> or (ij|kl)
       i(r1) j(r1) 1/r12 k(r2) l(r2)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_norm_ord_transp_cosgtos`
       * :c:data:`ao_expo_ord_transp_cosgtos`
       * :c:data:`ao_integrals_threshold`
       * :c:data:`ao_nucl`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`n_pt_max_integrals`
       * :c:data:`nucl_coord`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_explicit_cpoly_and_cgaussian`

 
.. c:function:: ao_idx2_sq:


    File : :file:`ao_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        subroutine ao_idx2_sq(i,j,ij)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`two_e_integrals_index_2fold`

 
.. c:function:: ao_idx2_sq_rev:


    File : :file:`ao_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        subroutine ao_idx2_sq_rev(i,k,ik)


    reverse square compound index

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`two_e_integrals_index_reverse_2fold`

 
.. c:function:: ao_idx2_tri_key:


    File : :file:`ao_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        subroutine ao_idx2_tri_key(i,j,ij)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`two_e_integrals_index_2fold`

 
.. c:function:: ao_idx2_tri_rev_key:


    File : :file:`ao_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        subroutine ao_idx2_tri_rev_key(i,k,ik)


    return i<=k

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`two_e_integrals_index_reverse_2fold`

 
.. c:function:: ao_l4:


    File : :file:`ao_two_e_ints/two_e_integrals.irp.f`

    .. code:: fortran

        integer function ao_l4(i,j,k,l)


    Computes the product of l values of i,j,k,and l

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_l`

 
.. c:function:: ao_two_e_integral:


    File : :file:`ao_two_e_ints/two_e_integrals.irp.f`

    .. code:: fortran

        double precision function ao_two_e_integral(i, j, k, l)


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
       * :c:data:`n_pt_max_integrals`
       * :c:data:`nucl_coord`
       * :c:data:`use_cgtos`
       * :c:data:`use_only_lr`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_explicit_poly_and_gaussian`

 
.. c:function:: ao_two_e_integral_cgtos:


    File : :file:`ao_two_e_ints/two_e_coul_integrals_cgtos.irp.f`

    .. code:: fortran

        double precision function ao_two_e_integral_cgtos(i, j, k, l)


    integral of the AO basis <ik|jl> or (ij|kl)
       i(r1) j(r1) 1/r12 k(r2) l(r2)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_cgtos_norm_ord_transp`
       * :c:data:`ao_expo_cgtos_ord_transp`
       * :c:data:`ao_nucl`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`n_pt_max_integrals`
       * :c:data:`nucl_coord`
       * :c:data:`use_pw`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_explicit_cpoly_and_cgaussian`

 
.. c:function:: ao_two_e_integral_cosgtos:


    File : :file:`ao_two_e_ints/two_e_Coul_integrals_cosgtos.irp.f`

    .. code:: fortran

        double precision function ao_two_e_integral_cosgtos(i, j, k, l)


    integral of the AO basis <ik|jl> or (ij|kl)
       i(r1) j(r1) 1/r12 k(r2) l(r2)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_norm_ord_transp_cosgtos`
       * :c:data:`ao_expo_ord_transp_cosgtos`
       * :c:data:`ao_nucl`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`n_pt_max_integrals`
       * :c:data:`nucl_coord`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_explicit_cpoly_and_cgaussian`

 
.. c:function:: ao_two_e_integral_erf:


    File : :file:`ao_two_e_ints/two_e_integrals_erf.irp.f`

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

 
.. c:function:: ao_two_e_integral_schwartz_accel:


    File : :file:`ao_two_e_ints/two_e_integrals.irp.f`

    .. code:: fortran

        double precision function ao_two_e_integral_schwartz_accel(i,j,k,l)


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

 
.. c:function:: ao_two_e_integral_schwartz_accel_erf:


    File : :file:`ao_two_e_ints/two_e_integrals_erf.irp.f`

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

 
.. c:function:: ao_two_e_integral_zero:


    File : :file:`ao_two_e_ints/screening.irp.f`

    .. code:: fortran

        logical function ao_two_e_integral_zero(i,j,k,l)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_threshold`
       * :c:data:`ao_overlap_abs`
       * :c:data:`ao_two_e_integral_schwartz`
       * :c:data:`is_periodic`
       * :c:data:`read_ao_two_e_integrals`
       * :c:data:`use_cgtos`

 
.. c:function:: ao_two_e_integrals_erf_in_map_collector:


    File : :file:`ao_two_e_ints/integrals_erf_in_map_slave.irp.f`

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


    File : :file:`ao_two_e_ints/integrals_erf_in_map_slave.irp.f`

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


    File : :file:`ao_two_e_ints/integrals_erf_in_map_slave.irp.f`

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


    File : :file:`ao_two_e_ints/integrals_erf_in_map_slave.irp.f`

    .. code:: fortran

        subroutine ao_two_e_integrals_erf_in_map_slave_tcp(i)


    Computes a buffer of integrals. i is the ID of the current thread.

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ao_two_e_integrals_erf_in_map_slave`

 
.. c:function:: ao_two_e_integrals_in_map_collector:


    File : :file:`ao_two_e_ints/integrals_in_map_slave.irp.f`

    .. code:: fortran

        subroutine ao_two_e_integrals_in_map_collector(zmq_socket_pull)


    Collects results from the AO integral calculation

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_map`
       * :c:data:`ao_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_integrals_in_map`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`end_zmq_to_qp_run_socket`
       * :c:func:`insert_into_ao_integrals_map`

 
.. c:function:: ao_two_e_integrals_in_map_slave:


    File : :file:`ao_two_e_ints/integrals_in_map_slave.irp.f`

    .. code:: fortran

        subroutine ao_two_e_integrals_in_map_slave(thread,iproc)


    Computes a buffer of integrals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ao_two_e_integrals_in_map_slave_inproc`
       * :c:func:`ao_two_e_integrals_in_map_slave_tcp`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`compute_ao_integrals_jl`
       * :c:func:`end_zmq_push_socket`
       * :c:func:`end_zmq_to_qp_run_socket`
       * :c:func:`push_integrals`
       * :c:func:`sscanf_dd`

 
.. c:function:: ao_two_e_integrals_in_map_slave_inproc:


    File : :file:`ao_two_e_ints/integrals_in_map_slave.irp.f`

    .. code:: fortran

        subroutine ao_two_e_integrals_in_map_slave_inproc(i)


    Computes a buffer of integrals. i is the ID of the current thread.

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_integrals_in_map`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ao_two_e_integrals_in_map_slave`

 
.. c:function:: ao_two_e_integrals_in_map_slave_tcp:


    File : :file:`ao_two_e_ints/integrals_in_map_slave.irp.f`

    .. code:: fortran

        subroutine ao_two_e_integrals_in_map_slave_tcp(i)


    Computes a buffer of integrals. i is the ID of the current thread.

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ao_two_e_integrals_in_map_slave`

 
.. c:function:: clear_ao_erf_map:


    File : :file:`ao_two_e_ints/map_integrals_erf.irp.f`

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

 
.. c:function:: clear_ao_map:


    File : :file:`ao_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        subroutine clear_ao_map


    Frees the memory of the AO map

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_map`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`map_deinit`

 
.. c:function:: compute_ao_integrals_erf_jl:


    File : :file:`ao_two_e_ints/two_e_integrals_erf.irp.f`

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

 
.. c:function:: compute_ao_integrals_jl:


    File : :file:`ao_two_e_ints/two_e_integrals.irp.f`

    .. code:: fortran

        subroutine compute_ao_integrals_jl(j,l,n_integrals,buffer_i,buffer_value)


    Parallel client for AO integrals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_threshold`
       * :c:data:`ao_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ao_two_e_integrals_in_map_slave`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`two_e_integrals_index`

 
.. c:function:: compute_ao_two_e_integrals:


    File : :file:`ao_two_e_ints/two_e_integrals.irp.f`

    .. code:: fortran

        subroutine compute_ao_two_e_integrals(j,k,l,sze,buffer_value)


    Compute AO 1/r12 integrals for all i and fixed j,k,l

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`

 
.. c:function:: compute_ao_two_e_integrals_erf:


    File : :file:`ao_two_e_ints/two_e_integrals_erf.irp.f`

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

 
.. c:function:: do_schwartz_accel:


    File : :file:`ao_two_e_ints/two_e_integrals.irp.f`

    .. code:: fortran

        logical function do_schwartz_accel(i,j,k,l)


    If true, use Schwatrz to accelerate direct integral calculation

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_prim_num`
       * :c:data:`do_ao_cholesky`

 
.. c:function:: dump_ao_integrals_erf:


    File : :file:`ao_two_e_ints/map_integrals_erf.irp.f`

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

 
.. c:function:: eri:


    File : :file:`ao_two_e_ints/two_e_integrals.irp.f`

    .. code:: fortran

        double precision function ERI(alpha,beta,delta,gama,a_x,b_x,c_x,d_x,a_y,b_y,c_y,d_y,a_z,b_z,c_z,d_z)


    ATOMIC PRIMTIVE two-electron integral between the 4 primitives ::
           primitive_1 = x1**(a_x) y1**(a_y) z1**(a_z) exp(-alpha * r1**2)
           primitive_2 = x1**(b_x) y1**(b_y) z1**(b_z) exp(- beta * r1**2)
           primitive_3 = x2**(c_x) y2**(c_y) z2**(c_z) exp(-delta * r2**2)
           primitive_4 = x2**(d_x) y2**(d_y) z2**(d_z) exp(- gama * r2**2)

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`integrale_new`

 
.. c:function:: eri_cgtos:


    File : :file:`ao_two_e_ints/two_e_coul_integrals_cgtos.irp.f`

    .. code:: fortran

        complex*16 function ERI_cgtos(alpha, beta, delta, gama, a_x, b_x, c_x, d_x, a_y, b_y, c_y, d_y, a_z, b_z, c_z, d_z)


    ATOMIC PRIMTIVE two-electron integral between the 4 primitives ::
           primitive_1 = x1**(a_x) y1**(a_y) z1**(a_z) exp(-alpha * r1**2)
           primitive_2 = x1**(b_x) y1**(b_y) z1**(b_z) exp(- beta * r1**2)
           primitive_3 = x2**(c_x) y2**(c_y) z2**(c_z) exp(-delta * r2**2)
           primitive_4 = x2**(d_x) y2**(d_y) z2**(d_z) exp(- gama * r2**2)

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`integrale_new_cgtos`

 
.. c:function:: eri_cosgtos:


    File : :file:`ao_two_e_ints/two_e_Coul_integrals_cosgtos.irp.f`

    .. code:: fortran

        complex*16 function ERI_cosgtos(alpha, beta, delta, gama, a_x, b_x, c_x, d_x, a_y, b_y, c_y, d_y, a_z, b_z, c_z, d_z)


    ATOMIC PRIMTIVE two-electron integral between the 4 primitives ::
           primitive_1 = x1**(a_x) y1**(a_y) z1**(a_z) exp(-alpha * r1**2)
           primitive_2 = x1**(b_x) y1**(b_y) z1**(b_z) exp(- beta * r1**2)
           primitive_3 = x2**(c_x) y2**(c_y) z2**(c_z) exp(-delta * r2**2)
           primitive_4 = x2**(d_x) y2**(d_y) z2**(d_z) exp(- gama * r2**2)

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`integrale_new_cosgtos`

 
.. c:function:: eri_erf:


    File : :file:`ao_two_e_ints/two_e_integrals_erf.irp.f`

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

 
.. c:function:: gauleg:


    File : :file:`ao_two_e_ints/gauss_legendre.irp.f`

    .. code:: fortran

        subroutine gauleg(x1,x2,x,w,n)


    Gauss-Legendre

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`gauleg_t2`

 
.. c:function:: get_ao_erf_map_size:


    File : :file:`ao_two_e_ints/map_integrals_erf.irp.f`

    .. code:: fortran

        function get_ao_erf_map_size()


    Returns the number of elements in the |AO| map

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_erf_map`

 
.. c:function:: get_ao_integ_chol:


    File : :file:`ao_two_e_ints/cholesky.irp.f`

    .. code:: fortran

        double precision function get_ao_integ_chol(i,j,k,l)


    CHOLESKY representation of the integral of the AO basis <ik|jl> or (ij|kl)
       i(r1) j(r1) 1/r12 k(r2) l(r2)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cholesky_ao_num`
       * :c:data:`cholesky_ao_transp`

 
.. c:function:: get_ao_map_size:


    File : :file:`ao_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        function get_ao_map_size()


    Returns the number of elements in the AO map

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_map`

 
.. c:function:: get_ao_two_e_integral:


    File : :file:`ao_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        double precision function get_ao_two_e_integral(i, j, k, l, map) result(result)


    Gets one AO bi-electronic integral from the AO map in PHYSICIST NOTATION
    
    <1:k, 2:l |1:i, 2:j>

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_cache`
       * :c:data:`ao_integrals_cache_min`
       * :c:data:`ao_two_e_integrals_in_map`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`map_get`
       * :c:func:`two_e_integrals_index`

 
.. c:function:: get_ao_two_e_integral_erf:


    File : :file:`ao_two_e_ints/map_integrals_erf.irp.f`

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

 
.. c:function:: get_ao_two_e_integral_periodic:


    File : :file:`ao_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        complex*16 function get_ao_two_e_integral_periodic(i,j,k,l,map) result(result)


    Gets one AO bi-electronic integral from the AO map

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_cache_min`
       * :c:data:`ao_integrals_cache_periodic`
       * :c:data:`ao_integrals_map`
       * :c:data:`ao_two_e_integrals_in_map`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`map_get`
       * :c:func:`two_e_integrals_index_2fold`

 
.. c:function:: get_ao_two_e_integrals:


    File : :file:`ao_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        subroutine get_ao_two_e_integrals(j,k,l,sze,out_val)


    Gets multiple AO bi-electronic integral from the AO map .
    All i are retrieved for j,k,l fixed.
    physicist convention : <ij|kl>

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_map`
       * :c:data:`ao_two_e_integrals_in_map`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`add_integrals_to_map`
       * :c:func:`four_idx_dgemm`

 
.. c:function:: get_ao_two_e_integrals_erf:


    File : :file:`ao_two_e_ints/map_integrals_erf.irp.f`

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
       * :c:func:`four_idx_dgemm_erf`

 
.. c:function:: get_ao_two_e_integrals_erf_non_zero:


    File : :file:`ao_two_e_ints/map_integrals_erf.irp.f`

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

 
.. c:function:: get_ao_two_e_integrals_non_zero:


    File : :file:`ao_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        subroutine get_ao_two_e_integrals_non_zero(j,k,l,sze,out_val,out_val_index,non_zero_int)


    Gets multiple AO bi-electronic integral from the AO map .
    All non-zero i are retrieved for j,k,l fixed.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_map`
       * :c:data:`ao_integrals_threshold`
       * :c:data:`ao_two_e_integrals_in_map`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`map_get`
       * :c:func:`two_e_integrals_index`

 
.. c:function:: get_ao_two_e_integrals_non_zero_jl:


    File : :file:`ao_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        subroutine get_ao_two_e_integrals_non_zero_jl(j,l,thresh,sze_max,sze,out_val,out_val_index,non_zero_int)


    Gets multiple AO bi-electronic integral from the AO map .
    All non-zero i are retrieved for j,k,l fixed.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_map`
       * :c:data:`ao_two_e_integrals_in_map`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`map_get`
       * :c:func:`two_e_integrals_index`

 
.. c:function:: get_ao_two_e_integrals_non_zero_jl_from_list:


    File : :file:`ao_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        subroutine get_ao_two_e_integrals_non_zero_jl_from_list(j,l,thresh,list,n_list,sze_max,out_val,out_val_index,non_zero_int)


    Gets multiple AO two-electron integrals from the AO map .
    All non-zero i are retrieved for j,k,l fixed.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_map`
       * :c:data:`ao_two_e_integrals_in_map`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`map_get`
       * :c:func:`two_e_integrals_index`

 
.. c:function:: get_ao_two_e_integrals_periodic:


    File : :file:`ao_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        subroutine get_ao_two_e_integrals_periodic(j,k,l,sze,out_val)


    Gets multiple AO bi-electronic integral from the AO map .
    All i are retrieved for j,k,l fixed.
    physicist convention : <ij|kl>

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_map`
       * :c:data:`ao_two_e_integrals_in_map`

 
.. c:function:: give_polynom_mult_center_x:


    File : :file:`ao_two_e_ints/two_e_integrals.irp.f`

    .. code:: fortran

        subroutine give_polynom_mult_center_x(P_center,Q_center,a_x,d_x,p,q,n_pt_in,pq_inv,pq_inv_2,p10_1,p01_1,p10_2,p01_2,d,n_pt_out)


    subroutine that returns the explicit polynom in term of the "t"
    variable of the following polynomw :
    
    $I_{x_1}(a_x,d_x,p,q) \, I_{x_1}(a_y,d_y,p,q) \ I_{x_1}(a_z,d_z,p,q)$

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`general_primitive_integral`
       * :c:func:`general_primitive_integral_erf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_pol_mult`

 
.. c:function:: i_x1_pol_mult:


    File : :file:`ao_two_e_ints/two_e_integrals.irp.f`

    .. code:: fortran

        subroutine I_x1_pol_mult(a,c,B_10,B_01,B_00,C_00,D_00,d,nd,n_pt_in)


    Recursive function involved in the two-electron integral

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`give_polynom_mult_center_x`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_pol_mult_a1`
       * :c:func:`i_x1_pol_mult_a2`
       * :c:func:`i_x1_pol_mult_recurs`
       * :c:func:`i_x2_pol_mult`

 
.. c:function:: i_x1_pol_mult_cgtos:


    File : :file:`ao_two_e_ints/two_e_coul_integrals_cgtos.irp.f`

    .. code:: fortran

        subroutine I_x1_pol_mult_cgtos(a, c, B_10, B_01, B_00, C_00, D_00, d, nd, n_pt_in)


    Recursive function involved in the two-electron integral

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`give_cpolynom_mult_center_x`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_pol_mult_a1_cgtos`
       * :c:func:`i_x1_pol_mult_a2_cgtos`
       * :c:func:`i_x1_pol_mult_recurs_cgtos`
       * :c:func:`i_x2_pol_mult_cgtos`

 
.. c:function:: i_x1_pol_mult_cosgtos:


    File : :file:`ao_two_e_ints/two_e_Coul_integrals_cosgtos.irp.f`

    .. code:: fortran

        subroutine I_x1_pol_mult_cosgtos(a, c, B_10, B_01, B_00, C_00, D_00, d, nd, n_pt_in)


    Recursive function involved in the two-electron integral

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`give_cpolynom_mult_center_x`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_pol_mult_a1_cosgtos`
       * :c:func:`i_x1_pol_mult_a2_cosgtos`
       * :c:func:`i_x1_pol_mult_recurs_cosgtos`
       * :c:func:`i_x2_pol_mult_cosgtos`

 
.. c:function:: idx2_tri_int:


    File : :file:`ao_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        subroutine idx2_tri_int(i,j,ij)



 
.. c:function:: idx2_tri_rev_int:


    File : :file:`ao_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        subroutine idx2_tri_rev_int(i,k,ik)


    return i<=k

 
.. c:function:: insert_into_ao_integrals_erf_map:


    File : :file:`ao_two_e_ints/map_integrals_erf.irp.f`

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

 
.. c:function:: insert_into_ao_integrals_map:


    File : :file:`ao_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        subroutine insert_into_ao_integrals_map(n_integrals,buffer_i, buffer_values)


    Create new entry into AO map

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_map`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ao_two_e_integrals_in_map_collector`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`map_append`

 
.. c:function:: integrale_new:


    File : :file:`ao_two_e_ints/two_e_integrals.irp.f`

    .. code:: fortran

        subroutine integrale_new(I_f,a_x,b_x,c_x,d_x,a_y,b_y,c_y,d_y,a_z,b_z,c_z,d_z,p,q,n_pt)


    Calculates the integral of the polynomial :
    
    $I_{x_1}(a_x+b_x,c_x+d_x,p,q) \, I_{x_1}(a_y+b_y,c_y+d_y,p,q) \, I_{x_1}(a_z+b_z,c_z+d_z,p,q)$
    in $( 0 ; 1)$

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`gauleg_t2`
       * :c:data:`n_pt_max_integrals`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`eri`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_new`

 
.. c:function:: integrale_new_cgtos:


    File : :file:`ao_two_e_ints/two_e_coul_integrals_cgtos.irp.f`

    .. code:: fortran

        subroutine integrale_new_cgtos(I_f, a_x, b_x, c_x, d_x, a_y, b_y, c_y, d_y, a_z, b_z, c_z, d_z, p, q, n_pt)


    Calculates the integral of the polynomial :
    
    $I_{x_1}(a_x+b_x, c_x+d_x, p, q) \, I_{x_1}(a_y+b_y, c_y+d_y, p, q) \, I_{x_1}(a_z+b_z, c_z+d_z, p, q)$
    in $( 0 ; 1)$

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`gauleg_t2`
       * :c:data:`n_pt_max_integrals`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`eri_cgtos`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_new_cgtos`

 
.. c:function:: integrale_new_cosgtos:


    File : :file:`ao_two_e_ints/two_e_Coul_integrals_cosgtos.irp.f`

    .. code:: fortran

        subroutine integrale_new_cosgtos(I_f, a_x, b_x, c_x, d_x, a_y, b_y, c_y, d_y, a_z, b_z, c_z, d_z, p, q, n_pt)


    Calculates the integral of the polynomial :
    
    $I_{x_1}(a_x+b_x, c_x+d_x, p, q) \, I_{x_1}(a_y+b_y, c_y+d_y, p, q) \, I_{x_1}(a_z+b_z, c_z+d_z, p, q)$
    in $( 0 ; 1)$

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`gauleg_t2`
       * :c:data:`n_pt_max_integrals`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`eri_cosgtos`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_new_cosgtos`

 
.. c:function:: integrale_new_erf:


    File : :file:`ao_two_e_ints/two_e_integrals_erf.irp.f`

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


    File : :file:`ao_two_e_ints/map_integrals_erf.irp.f`

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

 
.. c:function:: multiply_poly_c2_inline_2e:


    File : :file:`ao_two_e_ints/two_e_integrals.irp.f`

    .. code:: fortran

        subroutine multiply_poly_c2_inline_2e(b,nb,c,d,nd)


    Multiply two polynomials
    D(t) =! D(t) +( B(t)*C(t))

 
.. c:function:: multiply_poly_local:


    File : :file:`ao_two_e_ints/two_e_integrals.irp.f`

    .. code:: fortran

        subroutine multiply_poly_local(b,nb,c,nc,d,nd)


    Multiply two polynomials
    D(t) =! D(t) +( B(t)*C(t))

 
.. c:function:: n_pt_sup:


    File : :file:`ao_two_e_ints/two_e_integrals.irp.f`

    .. code:: fortran

        integer function n_pt_sup(a_x,b_x,c_x,d_x,a_y,b_y,c_y,d_y,a_z,b_z,c_z,d_z)


    Returns the upper boundary of the degree of the polynomial involved in the
    two-electron integral :
    
    $I_x(a_x,b_x,c_x,d_x) \, I_y(a_y,b_y,c_y,d_y) \, I_z(a_z,b_z,c_z,d_z)$

 
.. c:function:: push_integrals:


    File : :file:`ao_two_e_ints/integrals_in_map_slave.irp.f`

    .. code:: fortran

        subroutine push_integrals(zmq_socket_push, n_integrals, buffer_i, buffer_value, task_id)


    Push integrals in the push socket

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ao_two_e_integrals_erf_in_map_slave`
       * :c:func:`ao_two_e_integrals_in_map_slave`

 
.. c:function:: save_erf_two_e_integrals_ao:


    File : :file:`ao_two_e_ints/routines_save_integrals_erf.irp.f`

    .. code:: fortran

        subroutine save_erf_two_e_integrals_ao



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_erf_map`
       * :c:data:`ao_two_e_integrals_erf_in_map`
       * :c:data:`ezfio_filename`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_set_ao_two_e_ints_io_ao_two_e_integrals_erf`
       * :c:func:`ezfio_set_work_empty`
       * :c:func:`map_save_to_disk`

 
.. c:function:: save_erf_two_e_ints_ao_into_ints_ao:


    File : :file:`ao_two_e_ints/routines_save_integrals_erf.irp.f`

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

 
.. c:function:: two_e_integrals_index:


    File : :file:`ao_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        subroutine two_e_integrals_index(i,j,k,l,i1)


    Gives a unique index for i,j,k,l using permtuation symmetry.
    i <-> k, j <-> l, and (i,k) <-> (j,l) for non-periodic systems

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_cache`
       * :c:data:`ao_integrals_erf_cache`
       * :c:data:`ao_integrals_erf_map`
       * :c:data:`ao_integrals_map`
       * :c:func:`compute_ao_integrals_erf_jl`
       * :c:func:`compute_ao_integrals_jl`
       * :c:func:`get_ao_two_e_integral`
       * :c:func:`get_ao_two_e_integral_erf`
       * :c:func:`get_ao_two_e_integrals_erf_non_zero`
       * :c:func:`get_ao_two_e_integrals_non_zero`
       * :c:func:`get_ao_two_e_integrals_non_zero_jl`
       * :c:func:`get_ao_two_e_integrals_non_zero_jl_from_list`
       * :c:func:`get_mo_two_e_integral_erf`
       * :c:func:`get_mo_two_e_integrals_erf`
       * :c:func:`get_mo_two_e_integrals_erf_coulomb_ii`
       * :c:func:`get_mo_two_e_integrals_erf_exch_ii`
       * :c:func:`get_mo_two_e_integrals_erf_i1j1`
       * :c:func:`get_mo_two_e_integrals_erf_ij`
       * :c:func:`get_two_e_integral`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_integrals_erf_cache`
       * :c:data:`mo_integrals_erf_map`
       * :c:data:`mo_integrals_map`

 
.. c:function:: two_e_integrals_index_2fold:


    File : :file:`ao_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        subroutine two_e_integrals_index_2fold(i,j,k,l,i1)



    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_cache_periodic`
       * :c:func:`get_ao_two_e_integral_periodic`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ao_idx2_sq`
       * :c:func:`ao_idx2_tri_key`

 
.. c:function:: two_e_integrals_index_reverse:


    File : :file:`ao_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        subroutine two_e_integrals_index_reverse(i,j,k,l,i1)


    Computes the 4 indices $i,j,k,l$ from a unique index $i_1$.
    For 2 indices $i,j$ and $i \le j$, we have
    $p = i(i-1)/2 + j$.
    The key point is that because $j < i$,
    $i(i-1)/2 < p \le i(i+1)/2$. So $i$ can be found by solving
    $i^2 - i - 2p=0$. One obtains $i=1 + \sqrt{1+8p}/2$
    and $j = p - i(i-1)/2$.
    This rule is applied 3 times. First for the symmetry of the
    pairs (i,k) and (j,l), and then for the symmetry within each pair.

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_integral_alpha`

 
.. c:function:: two_e_integrals_index_reverse_2fold:


    File : :file:`ao_two_e_ints/map_integrals.irp.f`

    .. code:: fortran

        subroutine two_e_integrals_index_reverse_2fold(i,j,k,l,i1)



    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ao_idx2_sq_rev`
       * :c:func:`ao_idx2_tri_rev_key`

