.. _module_csf: 
 
.. program:: csf 
 
.. default-role:: option 
 
csf
===

No README.rst in this module; the providers and
subroutines below are documented from the irpf90 man pages. 
 
 
Providers 
--------- 
 
.. c:var:: aijpqcontainer


    File : :file:`csf/sigma_vector.irp.f`

    .. code:: fortran

        real*8, allocatable	:: aijpqcontainer	(0:NSOMOMax,0:NSOMOMax,4,NSOMOMax,NSOMOMax,NBFMax,NBFMax)


    Documentation for AIJpqMatrixList
    The prototype matrix containing the <I|E_{pq}|J>
    matrices for each I,J somo pair and orb ids.
    
    Due to the different types of excitations which
    include DOMOs and VMOs two prototype DOMOs and two
    prototype VMOs are needed. Therefore
    the 4th and 5th dimensions are NSOMOMax+4 and NSOMOMax+4
    respectively.
    
    The type of excitations are ordered as follows:
    Type 1 - SOMO -> SOMO
    Type 2 - DOMO -> VMO
    Type 3 - SOMO -> VMO
    Type 4 - DOMO -> SOMO

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`aijpqmatrixdimslist`
       * :c:data:`nsomomax`


 
.. c:var:: aijpqmatrixdimslist


    File : :file:`csf/sigma_vector.irp.f`

    .. code:: fortran

        integer, allocatable	:: aijpqmatrixdimslist	(0:NSOMOMax,0:NSOMOMax,4,NSOMOMax,NSOMOMax,2)
        integer	:: rowsmax	
        integer	:: colsmax	


    Documentation for AIJpqMatrixList
    The prototype matrix containing the <I|E_{pq}|J>
    matrices for each I,J somo pair and orb ids.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`nsomomax`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`aijpqcontainer`

 
.. c:var:: int cfg_nsomo_max


    File : :file:`csf/configurations.irp.f`

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


    File : :file:`csf/configurations.irp.f`

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


 
.. c:var:: colsmax


    File : :file:`csf/sigma_vector.irp.f`

    .. code:: fortran

        integer, allocatable	:: aijpqmatrixdimslist	(0:NSOMOMax,0:NSOMOMax,4,NSOMOMax,NSOMOMax,2)
        integer	:: rowsmax	
        integer	:: colsmax	


    Documentation for AIJpqMatrixList
    The prototype matrix containing the <I|E_{pq}|J>
    matrices for each I,J somo pair and orb ids.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`nsomomax`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`aijpqcontainer`

 
.. c:var:: int csf_unique_vector_index


    File : :file:`csf/csf.irp.f`

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

       * :c:data:`ci_electronic_energy`
       * :c:data:`psi_csf_coef`

 
.. c:var:: double csf_unique_vector_value


    File : :file:`csf/csf.irp.f`

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

       * :c:data:`ci_electronic_energy`
       * :c:data:`psi_csf_coef`

 
.. c:var:: int csf_unique_vector_value_size


    File : :file:`csf/csf.irp.f`

    .. code:: fortran

        integer	:: csf_unique_vector_value_size	


    Is reset in the provider of csf_unique_vector_value

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`n_csf`

 
.. c:var:: int det_to_configuration


    File : :file:`csf/configurations.irp.f`

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

 
.. c:var:: dettocsftransformationmatrix


    File : :file:`csf/sigma_vector.irp.f`

    .. code:: fortran

        real*8, allocatable	:: dettocsftransformationmatrix	(0:NSOMOMax,NBFMax,maxDetDimPerBF)
        real*8, allocatable	:: psi_coef_config	(n_CSF)
        integer, allocatable	:: psi_config_data	(N_configuration,2)


    Documentation for DetToCSFTransformationMatrix
    Provides the matrix of transformatons for the
    conversion between determinant to CSF basis (in BFs)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`binom`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`n_int`
       * :c:data:`nsomomax`
       * :c:data:`psi_coef`
       * :c:data:`psi_configuration`
       * :c:data:`psi_configuration_to_psi_det`
       * :c:data:`psi_det`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_csf_coef`

 
.. c:var:: int dominant_cfg


    File : :file:`csf/configurations.irp.f`

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


    File : :file:`csf/configurations.irp.f`

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


 
.. c:var:: maxdetdimperbf


    File : :file:`csf/sigma_vector.irp.f`

    .. code:: fortran

        integer	:: nsomomax	
        integer	:: ncsfmax	
        integer*8	:: nmo	
        integer	:: nbfmax	
        integer	:: n_csf	
        integer	:: maxdetdimperbf	


    Documentation for NSOMOMax
    The maximum number of SOMOs for the current calculation.
    required for the calculation of prototype arrays.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`binom`
       * :c:data:`cfg_seniority_index`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`elec_num`
       * :c:data:`n_act_orb`
       * :c:data:`psi_configuration`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`aijpqcontainer`
       * :c:data:`aijpqmatrixdimslist`
       * :c:data:`ci_electronic_energy`
       * :c:data:`dettocsftransformationmatrix`
       * :c:data:`psi_csf_coef`

 
.. c:var:: int n_configuration


    File : :file:`csf/configurations.irp.f`

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


    File : :file:`csf/csf.irp.f`

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

       * :c:data:`ci_electronic_energy`
       * :c:data:`psi_csf_coef`

 
.. c:var:: int n_det_per_config_max


    File : :file:`csf/configurations.irp.f`

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


    File : :file:`csf/configurations.irp.f`

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


    File : :file:`csf/configurations.irp.f`

    .. code:: fortran

        integer, allocatable	:: n_elec_alpha_for_psi_configuration	(N_configuration)


    Number of alpha electrons of configuration i (1..N_configuration),
    per its name. The code takes the determinant of psi_det at index
    psi_configuration_to_psi_det(2,i) and sums, over orbitals, the
    popcnt of the AND of its alpha bit strings with the configuration's
    singly occupied mask.
    The loop writes det_alpha(N_int) but should write
    det_alpha(l) = iand(det_tmp(l,1),psi_configuration(l,1,i)) -
    bug (author, 2026-10-31).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`psi_configuration`
       * :c:data:`psi_configuration_to_psi_det`
       * :c:data:`psi_det`


 
.. c:var:: nbfmax


    File : :file:`csf/sigma_vector.irp.f`

    .. code:: fortran

        integer	:: nsomomax	
        integer	:: ncsfmax	
        integer*8	:: nmo	
        integer	:: nbfmax	
        integer	:: n_csf	
        integer	:: maxdetdimperbf	


    Documentation for NSOMOMax
    The maximum number of SOMOs for the current calculation.
    required for the calculation of prototype arrays.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`binom`
       * :c:data:`cfg_seniority_index`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`elec_num`
       * :c:data:`n_act_orb`
       * :c:data:`psi_configuration`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`aijpqcontainer`
       * :c:data:`aijpqmatrixdimslist`
       * :c:data:`ci_electronic_energy`
       * :c:data:`dettocsftransformationmatrix`
       * :c:data:`psi_csf_coef`

 
.. c:var:: ncsfmax


    File : :file:`csf/sigma_vector.irp.f`

    .. code:: fortran

        integer	:: nsomomax	
        integer	:: ncsfmax	
        integer*8	:: nmo	
        integer	:: nbfmax	
        integer	:: n_csf	
        integer	:: maxdetdimperbf	


    Documentation for NSOMOMax
    The maximum number of SOMOs for the current calculation.
    required for the calculation of prototype arrays.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`binom`
       * :c:data:`cfg_seniority_index`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`elec_num`
       * :c:data:`n_act_orb`
       * :c:data:`psi_configuration`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`aijpqcontainer`
       * :c:data:`aijpqmatrixdimslist`
       * :c:data:`ci_electronic_energy`
       * :c:data:`dettocsftransformationmatrix`
       * :c:data:`psi_csf_coef`

 
.. c:var:: nmo


    File : :file:`csf/sigma_vector.irp.f`

    .. code:: fortran

        integer	:: nsomomax	
        integer	:: ncsfmax	
        integer*8	:: nmo	
        integer	:: nbfmax	
        integer	:: n_csf	
        integer	:: maxdetdimperbf	


    Documentation for NSOMOMax
    The maximum number of SOMOs for the current calculation.
    required for the calculation of prototype arrays.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`binom`
       * :c:data:`cfg_seniority_index`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`elec_num`
       * :c:data:`n_act_orb`
       * :c:data:`psi_configuration`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`aijpqcontainer`
       * :c:data:`aijpqmatrixdimslist`
       * :c:data:`ci_electronic_energy`
       * :c:data:`dettocsftransformationmatrix`
       * :c:data:`psi_csf_coef`

 
.. c:var:: nsomomax


    File : :file:`csf/sigma_vector.irp.f`

    .. code:: fortran

        integer	:: nsomomax	
        integer	:: ncsfmax	
        integer*8	:: nmo	
        integer	:: nbfmax	
        integer	:: n_csf	
        integer	:: maxdetdimperbf	


    Documentation for NSOMOMax
    The maximum number of SOMOs for the current calculation.
    required for the calculation of prototype arrays.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`binom`
       * :c:data:`cfg_seniority_index`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`elec_num`
       * :c:data:`n_act_orb`
       * :c:data:`psi_configuration`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`aijpqcontainer`
       * :c:data:`aijpqmatrixdimslist`
       * :c:data:`ci_electronic_energy`
       * :c:data:`dettocsftransformationmatrix`
       * :c:data:`psi_csf_coef`

 
.. c:var:: psi_coef_config


    File : :file:`csf/sigma_vector.irp.f`

    .. code:: fortran

        real*8, allocatable	:: dettocsftransformationmatrix	(0:NSOMOMax,NBFMax,maxDetDimPerBF)
        real*8, allocatable	:: psi_coef_config	(n_CSF)
        integer, allocatable	:: psi_config_data	(N_configuration,2)


    Documentation for DetToCSFTransformationMatrix
    Provides the matrix of transformatons for the
    conversion between determinant to CSF basis (in BFs)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`binom`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`n_int`
       * :c:data:`nsomomax`
       * :c:data:`psi_coef`
       * :c:data:`psi_configuration`
       * :c:data:`psi_configuration_to_psi_det`
       * :c:data:`psi_det`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_csf_coef`

 
.. c:var:: psi_config_data


    File : :file:`csf/sigma_vector.irp.f`

    .. code:: fortran

        real*8, allocatable	:: dettocsftransformationmatrix	(0:NSOMOMax,NBFMax,maxDetDimPerBF)
        real*8, allocatable	:: psi_coef_config	(n_CSF)
        integer, allocatable	:: psi_config_data	(N_configuration,2)


    Documentation for DetToCSFTransformationMatrix
    Provides the matrix of transformatons for the
    conversion between determinant to CSF basis (in BFs)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`binom`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`n_int`
       * :c:data:`nsomomax`
       * :c:data:`psi_coef`
       * :c:data:`psi_configuration`
       * :c:data:`psi_configuration_to_psi_det`
       * :c:data:`psi_det`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_csf_coef`

 
.. c:var:: int psi_configuration


    File : :file:`csf/configurations.irp.f`

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


    File : :file:`csf/configurations.irp.f`

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


    File : :file:`csf/csf.irp.f`

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

       * :c:data:`ci_electronic_energy`
       * :c:data:`psi_csf_coef`

 
.. c:var:: int psi_configuration_n_det


    File : :file:`csf/configurations.irp.f`

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


    File : :file:`csf/configurations.irp.f`

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


    File : :file:`csf/configurations.irp.f`

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


    File : :file:`csf/configurations.irp.f`

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


    File : :file:`csf/csf.irp.f`

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

       * :c:data:`ci_electronic_energy`
       * :c:data:`psi_csf_coef`

 
.. c:var:: int psi_configuration_to_psi_det


    File : :file:`csf/configurations.irp.f`

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


    File : :file:`csf/configurations.irp.f`

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


    File : :file:`csf/csf.irp.f`

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


 
.. c:var:: rowsmax


    File : :file:`csf/sigma_vector.irp.f`

    .. code:: fortran

        integer, allocatable	:: aijpqmatrixdimslist	(0:NSOMOMax,0:NSOMOMax,4,NSOMOMax,NSOMOMax,2)
        integer	:: rowsmax	
        integer	:: colsmax	


    Documentation for AIJpqMatrixList
    The prototype matrix containing the <I|E_{pq}|J>
    matrices for each I,J somo pair and orb ids.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`nsomomax`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`aijpqcontainer`

 
.. c:var:: int spin_multiplicity


    File : :file:`csf/configurations.irp.f`

    .. code:: fortran

        integer	:: spin_multiplicity	


    n_alpha - n_beta + 1

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`


 
.. c:var:: double weight_configuration


    File : :file:`csf/configurations.irp.f`

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


    File : :file:`csf/configurations.irp.f`

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


    File : :file:`csf/configurations.irp.f`

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


    File : :file:`csf/configurations.irp.f`

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


    File : :file:`csf/configurations.irp.f`

    .. code:: fortran

        subroutine configuration_of_det(d,o,Nint)


    Transforms a determinant to a configuration
    
    occ(:,1) : Single occupations
    
    occ(:,2) : Double occupations
    

 
.. c:function:: void configuration_search_key(void)


    File : :file:`csf/connected_to_ref.irp.f`

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


    File : :file:`csf/configurations.irp.f`

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


    File : :file:`csf/configurations.irp.f`

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

 
.. c:function:: configuration_to_dets_tree_addressing:


    File : :file:`csf/configurations.irp.f`

    .. code:: fortran

        subroutine configuration_to_dets_tree_addressing(o,d,sze,n_alpha,Nint)


    Generate all possible determinants for a given configuration
    
    This function preserves the tree addressing i.e.
    the time-reversal determinants are at the opposite ends
    and not one after the other as in the parent function.
    
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

       * :c:data:`psi_configuration_to_psi_det`

 
.. c:function:: convertorbidstomodelspaceids:


    File : :file:`csf/configuration_CI_sigma_helpers.irp.f`

    .. code:: fortran

        subroutine convertOrbIdsToModelSpaceIds(Ialpha, Jcfg, p, q, extype, pmodel, qmodel)


    This function converts the orbital ids
    in real space to those used in model space
    in order to identify the matrices required
    for the calculation of MEs.
    
    The type of excitations are ordered as follows:
    Type 1 - SOMO -> SOMO
    Type 2 - DOMO -> VMO
    Type 3 - SOMO -> VMO
    Type 4 - DOMO -> SOMO

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

 
.. c:function:: void convertwffromcsftodet(void)


    File : :file:`csf/csf.irp.f`

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

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`davidson_diag_csf_hjj`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`get_det_csf_transformation`

 
.. c:function:: void convertwffromdettocsf(void)


    File : :file:`csf/csf.irp.f`

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

       * :c:func:`davidson_diag_csf_hjj`
       * :c:data:`psi_csf_coef`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`get_det_csf_transformation`

 
.. c:function:: do_single_excitation_cfg:


    File : :file:`csf/create_excitations.irp.f`

    .. code:: fortran

        subroutine do_single_excitation_cfg(key_in,key_out,i_hole,i_particle,ok)


    Applies the single excitation operator to a configuration
    If the excitation is possible, ok is True

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`generate_all_singles_cfg`

 
.. c:function:: do_single_excitation_cfg_with_type:


    File : :file:`csf/create_excitations.irp.f`

    .. code:: fortran

        subroutine do_single_excitation_cfg_with_type(key_in,key_out,i_hole,i_particle,ex_type,ok)


    Applies the single excitation operator to a configuration
    Returns the type of excitation in ex_type
    where the following convention is used
    1 = (SOMO -> SOMO) 1 change in Nsomo
    2 = (DOMO ->  VMO) 1 change in Nsomo
    3 = (SOMO ->  VMO) 0 change in Nsomo
    4 = (DOMO -> SOMO) 0 change in Nsomo
    If the excitation is possible, ok is True

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`generate_all_singles_cfg_with_type`

 
.. c:function:: generate_all_singles_cfg:


    File : :file:`csf/create_excitations.irp.f`

    .. code:: fortran

        subroutine generate_all_singles_cfg(cfg,singles,n_singles,Nint)


    Generate all single excitation wrt a configuration
    
    n_singles : on input, max number of singles :
      elec_alpha_num * (mo_num - elec_beta_num)
      on output, number of generated singles

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`do_single_excitation_cfg`

 
.. c:function:: generate_all_singles_cfg_with_type:


    File : :file:`csf/create_excitations.irp.f`

    .. code:: fortran

        subroutine generate_all_singles_cfg_with_type(cfgInp,singles,idxs_singles,pq_singles,ex_type_singles,n_singles,Nint)


    Generate all single excitation wrt a configuration
    
    n_singles : on input, max number of singles :
      elec_alpha_num * (mo_num - elec_beta_num)
      on output, number of generated singles
    ex_type_singles : on output contains type of excitations  :
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_act_orb`
       * :c:data:`n_core_orb`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`binary_search_cfg`
       * :c:func:`do_single_excitation_cfg_with_type`

 
.. c:function:: void get_det_csf_transformation(void)


    File : :file:`csf/csf.irp.f`

    .. code:: fortran

        subroutine get_det_csf_transformation(matrix, sze, i_cfg)


    Fills matrix with the det-csf transformation matrix of
    configuration i_cfg: column j (j = 1..ncsf) is the j-th CSF of the
    configuration as a vector of ndet components, in the determinant
    basis, the determinants being ordered as in
    psi_configuration_to_psi_det_data. For a single-determinant
    configuration, matrix(1,1) = 1. sze must be >= ndet and >= ncsf.

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

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`qp_bug`

 
.. c:function:: get_phase_qp_to_cfg:


    File : :file:`csf/sigma_vector.irp.f`

    .. code:: fortran

        subroutine get_phase_qp_to_cfg(Ialpha, Ibeta, phaseout)


    Documentation for get_phase_qp_to_cfg
    
    This function converts from (aaaa)(bbbb)
    notation to (ab)(ab)(ab)(ab)
    notation.
    The cfgCI code works in (ab)(ab)(ab)(ab)
    notation throughout.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`convertwffromcsftodet`
       * :c:func:`convertwffromdettocsf`
       * :c:data:`dettocsftransformationmatrix`

 
.. c:function:: getnsomo:


    File : :file:`csf/configuration_CI_sigma_helpers.irp.f`

    .. code:: fortran

        function getNSOMO(Icfg) result(NSOMO)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

 
.. c:function:: logabsgamma:


    File : :file:`csf/sigma_vector.irp.f`

    .. code:: fortran

        real*8 function logabsgamma(x)



 
.. c:function:: void make_s2_eigenfunction(void)


    File : :file:`csf/configurations.irp.f`

    .. code:: fortran

        subroutine make_s2_eigenfunction


    Adds to the wavefunction the determinants missing from each of its
    configurations: all determinants of each configuration in
    psi_configuration are generated, and the missing ones are added
    with zero coefficient through the H-apply buffer. Used by CIPSI
    when s2_eig is true, so that the eigenstates of H can also be
    eigenstates of S^2.

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

