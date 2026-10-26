.. _module_csf: 
 
.. program:: csf 
 
.. default-role:: option 
 
csf
===

No README.rst in this module; the providers and
subroutines below are documented from the irpf90 man pages. 
 
 
Programs 
-------- 
 
 * :ref:`print_csf_map` 
 
Providers 
--------- 
 
.. c:var:: int cfg_nsomo_max


    File : :file:`configurations.irp.f`

    .. code:: fortran

        integer, allocatable	:: cfg_seniority_index	(0:elec_num)
        integer	:: cfg_nsomo_max	


    Returns the index in psi_configuration of the first cfg with
    the requested seniority
    
    cfg_nsomo_max : Max number of SOMO in the current wave function

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_num`
       * :c:data:`n_int`
       * :c:data:`psi_configuration`


 
.. c:var:: int cfg_seniority_index


    File : :file:`configurations.irp.f`

    .. code:: fortran

        integer, allocatable	:: cfg_seniority_index	(0:elec_num)
        integer	:: cfg_nsomo_max	


    Returns the index in psi_configuration of the first cfg with
    the requested seniority
    
    cfg_nsomo_max : Max number of SOMO in the current wave function

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_num`
       * :c:data:`n_int`
       * :c:data:`psi_configuration`


 
.. c:var:: int csf_unique_vector_index


    File : :file:`csf.irp.f`

    .. code:: fortran

        integer	:: n_csf	
        integer, allocatable	:: psi_configuration_to_psi_csf	(2,N_configuration)
        integer, allocatable	:: psi_configuration_n_csf	(N_configuration)
        integer, allocatable	:: csf_unique_vector_index	(N_det)
        double precision, allocatable	:: csf_unique_vector_value	(csf_unique_vector_value_size)


    To get the det-csf transformation matrix, call get_det_csf_transformation

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`csf_unique_vector_value_size`
       * :c:data:`expected_s2`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_configuration`
       * :c:data:`psi_configuration_to_psi_det`
       * :c:data:`psi_det`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_csf_coef`

 
.. c:var:: double csf_unique_vector_value


    File : :file:`csf.irp.f`

    .. code:: fortran

        integer	:: n_csf	
        integer, allocatable	:: psi_configuration_to_psi_csf	(2,N_configuration)
        integer, allocatable	:: psi_configuration_n_csf	(N_configuration)
        integer, allocatable	:: csf_unique_vector_index	(N_det)
        double precision, allocatable	:: csf_unique_vector_value	(csf_unique_vector_value_size)


    To get the det-csf transformation matrix, call get_det_csf_transformation

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`csf_unique_vector_value_size`
       * :c:data:`expected_s2`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_configuration`
       * :c:data:`psi_configuration_to_psi_det`
       * :c:data:`psi_det`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_csf_coef`

 
.. c:var:: int csf_unique_vector_value_size


    File : :file:`csf.irp.f`

    .. code:: fortran

        integer	:: csf_unique_vector_value_size	


    Is reset in the provider of csf_unique_vector_value

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`n_csf`

 
.. c:var:: int det_to_configuration


    File : :file:`configurations.irp.f`

    .. code:: fortran

        integer, allocatable	:: det_to_configuration	(N_det)


    Returns the index of the configuration for each determinant

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_configuration`
       * :c:data:`psi_det`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`dominant_cfg`
       * :c:data:`psi_configuration_hii`
       * :c:data:`psi_configuration_to_psi_det`
       * :c:data:`weight_configuration`
       * :c:data:`weight_configuration_average`

 
.. c:var:: int dominant_cfg


    File : :file:`configurations.irp.f`

    .. code:: fortran

        integer, allocatable	:: dominant_cfg	(N_states)


    Configuration of the determinants with the largest weight, for each state

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`det_to_configuration`
       * :c:data:`dominant_det`
       * :c:data:`n_det`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`dominant_dets_of_cfgs`
       * :c:data:`n_dominant_dets_of_cfgs`

 
.. c:var:: int dominant_dets_of_cfgs


    File : :file:`configurations.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: dominant_dets_of_cfgs	(N_int,2,N_dominant_dets_of_cfgs)


    Configuration of the determinants with the largest weight, for each state

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`binom_int`
       * :c:data:`dominant_cfg`
       * :c:data:`elec_alpha_num`
       * :c:data:`n_dominant_dets_of_cfgs`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_configuration`


 
.. c:var:: int n_configuration


    File : :file:`configurations.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_configuration	(N_int,2,psi_det_size)
        integer	:: n_configuration	


    Array of the configurations present in the wave function.
    
    psi_configuration(:,1,j) = j-th configuration of the wave function : represents all the single occupations
    
    psi_configuration(:,2,j) = j-th configuration of the wave function : represents all the double occupations
    
    The occ patterns are sorted by :c:func:`configuration_search_key`

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cfg_seniority_index`
       * :c:data:`det_to_configuration`
       * :c:data:`dominant_dets_of_cfgs`
       * :c:data:`n_csf`
       * :c:data:`n_dominant_dets_of_cfgs`
       * :c:data:`n_elec_alpha_for_psi_configuration`
       * :c:data:`psi_configuration_hii`
       * :c:data:`psi_configuration_sorted`
       * :c:data:`psi_configuration_to_psi_det`
       * :c:data:`psi_csf_coef`
       * :c:data:`weight_configuration`
       * :c:data:`weight_configuration_average`

 
.. c:var:: int n_csf


    File : :file:`csf.irp.f`

    .. code:: fortran

        integer	:: n_csf	
        integer, allocatable	:: psi_configuration_to_psi_csf	(2,N_configuration)
        integer, allocatable	:: psi_configuration_n_csf	(N_configuration)
        integer, allocatable	:: csf_unique_vector_index	(N_det)
        double precision, allocatable	:: csf_unique_vector_value	(csf_unique_vector_value_size)


    To get the det-csf transformation matrix, call get_det_csf_transformation

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`csf_unique_vector_value_size`
       * :c:data:`expected_s2`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_configuration`
       * :c:data:`psi_configuration_to_psi_det`
       * :c:data:`psi_det`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_csf_coef`

 
.. c:var:: int n_det_per_config_max


    File : :file:`configurations.irp.f`

    .. code:: fortran

        integer, allocatable	:: psi_configuration_to_psi_det	(2,N_configuration)
        integer, allocatable	:: psi_configuration_n_det	(N_configuration)
        integer, allocatable	:: psi_configuration_to_psi_det_data	(N_det)
        integer	:: n_det_per_config_max	


    psi_configuration_to_psi_det_data(k) -> i : i is the index of the
    determinant in psi_det.
    
    psi_configuration_to_psi_det(1:2,k) gives the first and last index of the
    determinants of configuration k in array psi_configuration_to_psi_det_data.
    
    n_det_per_config_max: Maximum number of determinants per configuration

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`binom_int`
       * :c:data:`det_to_configuration`
       * :c:data:`elec_alpha_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_configuration`
       * :c:data:`psi_det`
       * :c:data:`qp_max_mem`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`n_csf`
       * :c:data:`n_elec_alpha_for_psi_configuration`
       * :c:data:`psi_csf_coef`

 
.. c:var:: int n_dominant_dets_of_cfgs


    File : :file:`configurations.irp.f`

    .. code:: fortran

        integer	:: n_dominant_dets_of_cfgs	


    Number of determinants in all dominant determinants

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`binom_int`
       * :c:data:`dominant_cfg`
       * :c:data:`elec_alpha_num`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_configuration`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`dominant_dets_of_cfgs`

 
.. c:var:: int n_elec_alpha_for_psi_configuration


    File : :file:`configurations.irp.f`

    .. code:: fortran

        integer, allocatable	:: n_elec_alpha_for_psi_configuration	(N_configuration)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`psi_configuration`
       * :c:data:`psi_configuration_to_psi_det`
       * :c:data:`psi_det`


 
.. c:var:: int psi_configuration


    File : :file:`configurations.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_configuration	(N_int,2,psi_det_size)
        integer	:: n_configuration	


    Array of the configurations present in the wave function.
    
    psi_configuration(:,1,j) = j-th configuration of the wave function : represents all the single occupations
    
    psi_configuration(:,2,j) = j-th configuration of the wave function : represents all the double occupations
    
    The occ patterns are sorted by :c:func:`configuration_search_key`

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cfg_seniority_index`
       * :c:data:`det_to_configuration`
       * :c:data:`dominant_dets_of_cfgs`
       * :c:data:`n_csf`
       * :c:data:`n_dominant_dets_of_cfgs`
       * :c:data:`n_elec_alpha_for_psi_configuration`
       * :c:data:`psi_configuration_hii`
       * :c:data:`psi_configuration_sorted`
       * :c:data:`psi_configuration_to_psi_det`
       * :c:data:`psi_csf_coef`
       * :c:data:`weight_configuration`
       * :c:data:`weight_configuration_average`

 
.. c:var:: double psi_configuration_hii


    File : :file:`configurations.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_configuration_hii	(N_configuration)


    :math:`\langle I|H|I \rangle`  where :math:`|I\rangle`  is a configuration.
    This is the minimum :math:`H_{ii}` , where the :math:`|i\rangle`  are the
    determinants of :math:`|I\rangle` .

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`det_to_configuration`
       * :c:data:`n_det`
       * :c:data:`psi_configuration`
       * :c:data:`psi_det_hii`


 
.. c:var:: int psi_configuration_n_csf


    File : :file:`csf.irp.f`

    .. code:: fortran

        integer	:: n_csf	
        integer, allocatable	:: psi_configuration_to_psi_csf	(2,N_configuration)
        integer, allocatable	:: psi_configuration_n_csf	(N_configuration)
        integer, allocatable	:: csf_unique_vector_index	(N_det)
        double precision, allocatable	:: csf_unique_vector_value	(csf_unique_vector_value_size)


    To get the det-csf transformation matrix, call get_det_csf_transformation

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`csf_unique_vector_value_size`
       * :c:data:`expected_s2`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_configuration`
       * :c:data:`psi_configuration_to_psi_det`
       * :c:data:`psi_det`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_csf_coef`

 
.. c:var:: int psi_configuration_n_det


    File : :file:`configurations.irp.f`

    .. code:: fortran

        integer, allocatable	:: psi_configuration_to_psi_det	(2,N_configuration)
        integer, allocatable	:: psi_configuration_n_det	(N_configuration)
        integer, allocatable	:: psi_configuration_to_psi_det_data	(N_det)
        integer	:: n_det_per_config_max	


    psi_configuration_to_psi_det_data(k) -> i : i is the index of the
    determinant in psi_det.
    
    psi_configuration_to_psi_det(1:2,k) gives the first and last index of the
    determinants of configuration k in array psi_configuration_to_psi_det_data.
    
    n_det_per_config_max: Maximum number of determinants per configuration

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`binom_int`
       * :c:data:`det_to_configuration`
       * :c:data:`elec_alpha_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_configuration`
       * :c:data:`psi_det`
       * :c:data:`qp_max_mem`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`n_csf`
       * :c:data:`n_elec_alpha_for_psi_configuration`
       * :c:data:`psi_csf_coef`

 
.. c:var:: int psi_configuration_sorted


    File : :file:`configurations.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_configuration_sorted	(N_int,2,N_configuration)
        double precision, allocatable	:: weight_configuration_average_sorted	(N_configuration)
        integer, allocatable	:: psi_configuration_sorted_order	(N_configuration)
        integer, allocatable	:: psi_configuration_sorted_order_reverse	(N_configuration)


    Configurations sorted by weight

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`psi_configuration`
       * :c:data:`weight_configuration_average`


 
.. c:var:: int psi_configuration_sorted_order


    File : :file:`configurations.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_configuration_sorted	(N_int,2,N_configuration)
        double precision, allocatable	:: weight_configuration_average_sorted	(N_configuration)
        integer, allocatable	:: psi_configuration_sorted_order	(N_configuration)
        integer, allocatable	:: psi_configuration_sorted_order_reverse	(N_configuration)


    Configurations sorted by weight

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`psi_configuration`
       * :c:data:`weight_configuration_average`


 
.. c:var:: int psi_configuration_sorted_order_reverse


    File : :file:`configurations.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_configuration_sorted	(N_int,2,N_configuration)
        double precision, allocatable	:: weight_configuration_average_sorted	(N_configuration)
        integer, allocatable	:: psi_configuration_sorted_order	(N_configuration)
        integer, allocatable	:: psi_configuration_sorted_order_reverse	(N_configuration)


    Configurations sorted by weight

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`psi_configuration`
       * :c:data:`weight_configuration_average`


 
.. c:var:: int psi_configuration_to_psi_csf


    File : :file:`csf.irp.f`

    .. code:: fortran

        integer	:: n_csf	
        integer, allocatable	:: psi_configuration_to_psi_csf	(2,N_configuration)
        integer, allocatable	:: psi_configuration_n_csf	(N_configuration)
        integer, allocatable	:: csf_unique_vector_index	(N_det)
        double precision, allocatable	:: csf_unique_vector_value	(csf_unique_vector_value_size)


    To get the det-csf transformation matrix, call get_det_csf_transformation

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`csf_unique_vector_value_size`
       * :c:data:`expected_s2`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_configuration`
       * :c:data:`psi_configuration_to_psi_det`
       * :c:data:`psi_det`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_csf_coef`

 
.. c:var:: int psi_configuration_to_psi_det


    File : :file:`configurations.irp.f`

    .. code:: fortran

        integer, allocatable	:: psi_configuration_to_psi_det	(2,N_configuration)
        integer, allocatable	:: psi_configuration_n_det	(N_configuration)
        integer, allocatable	:: psi_configuration_to_psi_det_data	(N_det)
        integer	:: n_det_per_config_max	


    psi_configuration_to_psi_det_data(k) -> i : i is the index of the
    determinant in psi_det.
    
    psi_configuration_to_psi_det(1:2,k) gives the first and last index of the
    determinants of configuration k in array psi_configuration_to_psi_det_data.
    
    n_det_per_config_max: Maximum number of determinants per configuration

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`binom_int`
       * :c:data:`det_to_configuration`
       * :c:data:`elec_alpha_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_configuration`
       * :c:data:`psi_det`
       * :c:data:`qp_max_mem`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`n_csf`
       * :c:data:`n_elec_alpha_for_psi_configuration`
       * :c:data:`psi_csf_coef`

 
.. c:var:: int psi_configuration_to_psi_det_data


    File : :file:`configurations.irp.f`

    .. code:: fortran

        integer, allocatable	:: psi_configuration_to_psi_det	(2,N_configuration)
        integer, allocatable	:: psi_configuration_n_det	(N_configuration)
        integer, allocatable	:: psi_configuration_to_psi_det_data	(N_det)
        integer	:: n_det_per_config_max	


    psi_configuration_to_psi_det_data(k) -> i : i is the index of the
    determinant in psi_det.
    
    psi_configuration_to_psi_det(1:2,k) gives the first and last index of the
    determinants of configuration k in array psi_configuration_to_psi_det_data.
    
    n_det_per_config_max: Maximum number of determinants per configuration

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`binom_int`
       * :c:data:`det_to_configuration`
       * :c:data:`elec_alpha_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_configuration`
       * :c:data:`psi_det`
       * :c:data:`qp_max_mem`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`n_csf`
       * :c:data:`n_elec_alpha_for_psi_configuration`
       * :c:data:`psi_csf_coef`

 
.. c:var:: double psi_csf_coef


    File : :file:`csf.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_csf_coef	(N_csf,N_states)


    Wafe function in CSF basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_csf`
       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`psi_coef`
       * :c:data:`psi_configuration`
       * :c:data:`psi_configuration_to_psi_det`


 
.. c:var:: int spin_multiplicity


    File : :file:`configurations.irp.f`

    .. code:: fortran

        integer	:: spin_multiplicity	


    n_alpha - n_beta + 1

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`


 
.. c:var:: double weight_configuration


    File : :file:`configurations.irp.f`

    .. code:: fortran

        double precision, allocatable	:: weight_configuration	(N_configuration,N_states)


    Weight of the configurations in the wave function

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`det_to_configuration`
       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`psi_coef`
       * :c:data:`psi_configuration`


 
.. c:var:: double weight_configuration_average


    File : :file:`configurations.irp.f`

    .. code:: fortran

        double precision, allocatable	:: weight_configuration_average	(N_configuration)


    State-average weight of the configurations in the wave function

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`det_to_configuration`
       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`psi_coef`
       * :c:data:`psi_configuration`
       * :c:data:`state_average_weight`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_configuration_sorted`

 
.. c:var:: double weight_configuration_average_sorted


    File : :file:`configurations.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_configuration_sorted	(N_int,2,N_configuration)
        double precision, allocatable	:: weight_configuration_average_sorted	(N_configuration)
        integer, allocatable	:: psi_configuration_sorted_order	(N_configuration)
        integer, allocatable	:: psi_configuration_sorted_order_reverse	(N_configuration)


    Configurations sorted by weight

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`psi_configuration`
       * :c:data:`weight_configuration_average`


 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: void binary_search_cfg(void)


    File : :file:`configurations.irp.f`

    .. code:: fortran

        subroutine binary_search_cfg(cfgInp,addcfg)


    Documentation for binary_search
    
    Does a binary search to find
    the address of a configuration in a list of
    configurations.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`psi_configuration`

 
.. c:function:: void configuration_of_det(void)


    File : :file:`configurations.irp.f`

    .. code:: fortran

        subroutine configuration_of_det(d,o,Nint)


    Transforms a determinant to a configuration
    
    occ(:,1) : Single occupations
    
    occ(:,2) : Double occupations
    

 
.. c:function:: void configuration_search_key(void)


    File : :file:`connected_to_ref.irp.f`

    .. code:: fortran

        integer*8 function configuration_search_key(cfg,Nint)


    Returns an integer*8 corresponding to a determinant index for searching.
    The left-most 8 bits contain the number of open shells+1. This ensures that the CSF
    are packed with the same seniority.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`

 
.. c:function:: void configuration_to_dets(void)


    File : :file:`configurations.irp.f`

    .. code:: fortran

        subroutine configuration_to_dets(o,d,sze,n_alpha,Nint)


    Generate all possible determinants for a given configuration
    
    Input :
       o   : configuration : (doubly occupied, singly occupied)
       sze : Number of produced determinants, computed by `configuration_to_dets_size`
       n_alpha : Number of $\alpha$ electrons
       Nint    : N_int
    
    Output:
       d : determinants
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`binom_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`dominant_dets_of_cfgs`
       * :c:func:`make_s2_eigenfunction`
       * :c:data:`psi_configuration_to_psi_det`

 
.. c:function:: void configuration_to_dets_size(void)


    File : :file:`configurations.irp.f`

    .. code:: fortran

        subroutine configuration_to_dets_size(o,sze,n_alpha,Nint)


    Number of possible determinants for a given configuration

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`binom_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`make_s2_eigenfunction`
       * :c:data:`n_dominant_dets_of_cfgs`

 
.. c:function:: void convertwffromcsftodet(void)


    File : :file:`csf.irp.f`

    .. code:: fortran

        subroutine convertWFfromCSFtoDET(N_st,psi_coef_csf_in, psi_coef_det_out)


    Converts a wave functions from determinant space to CSF space

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_csf`
       * :c:data:`n_det`
       * :c:data:`psi_configuration`
       * :c:data:`psi_configuration_to_psi_det`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`get_det_csf_transformation`

 
.. c:function:: void convertwffromdettocsf(void)


    File : :file:`csf.irp.f`

    .. code:: fortran

        subroutine convertWFfromDETtoCSF(N_st,psi_coef_det_in, psi_coef_csf_out)


    Converts a wave functions from determinant space to CSF space

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_csf`
       * :c:data:`n_det`
       * :c:data:`psi_configuration`
       * :c:data:`psi_configuration_to_psi_det`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_csf_coef`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`get_det_csf_transformation`

 
.. c:function:: void get_det_csf_transformation(void)


    File : :file:`csf.irp.f`

    .. code:: fortran

        subroutine get_det_csf_transformation(matrix, sze, i_cfg)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_csf`
       * :c:data:`psi_configuration_to_psi_det`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`convertwffromcsftodet`
       * :c:func:`convertwffromdettocsf`
       * :c:func:`print_csf_map`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`qp_bug`

 
.. c:function:: void make_s2_eigenfunction(void)


    File : :file:`configurations.irp.f`

    .. code:: fortran

        subroutine make_s2_eigenfunction



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`psi_configuration`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`psi_coef`
       * :c:data:`psi_configuration`
       * :c:data:`psi_det`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`configuration_to_dets`
       * :c:func:`configuration_to_dets_size`
       * :c:func:`copy_h_apply_buffer_to_wf`
       * :c:func:`fill_h_apply_buffer_no_selection`
       * :c:func:`write_int`
       * :c:func:`write_time`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`psi_configuration`
       * :c:data:`n_det`
       * :c:data:`c0_weight`
       * :c:data:`psi_coef`
       * :c:data:`psi_det_sorted_bit`
       * :c:data:`psi_det_sorted_bit`
       * :c:data:`psi_configuration`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`
       * :c:data:`psi_det_sorted_bit`

