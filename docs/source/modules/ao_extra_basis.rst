.. _module_ao_extra_basis: 
 
.. program:: ao_extra_basis 
 
.. default-role:: option 
 
===========
extra_basis
===========

Plugin to handle an extra basis, which is attached to the extra_nuclei.
It is essentially a duplication of all important quantities (coefficients, exponents and so on) of the usual |AO| basis. 

Check in the directory "tuto" for a simple example of how to create a fictious system "B" attached independently to a system "A"

Another interesting feature is the possibility to fit any basis made at most with "p" functions onto a purely "s" basis. 
This is done with the various scripts here: 
 
 - qp_fit_1s_basis : script that creates an |EZFIO| folder corresponding to an .xyz file and a basis fitted with only "s" functions
 - qp_add_extra_fit_system : script that takes as input an |EZFIO| folder and an .xyz file and add an extra_basis and extra_nuclei with a 1s fitted basis
 
Ex: 
qp_add_extra_fit_system LiH.ezfio/ h2o.xyz # takes the EZFIO folder "LiH.ezfio" and creates all necessary additional basis and nuclei based on h2o.xyz, but only with 1s functions. 

 
 
 
EZFIO parameters 
---------------- 
 
.. option:: ao_extra_basis
 
    Name of the |ao_extra| basis set
 
 
.. option:: ao_extra_only_1s
 
    If |true|, you know that the additional AO basis is built only with 1s functions
 
    Default: true
 
.. option:: ao_extra_num
 
    Number of |ao_extras|
 
 
.. option:: ao_extra_prim_num
 
    Number of primitives per |ao_extra|
 
 
.. option:: ao_extra_prim_num_max
 
    Maximum number of primitives
 
    Default: =maxval(ao_extra_basis.ao_extra_prim_num)
 
.. option:: ao_extra_nucl
 
    Index of the nucleus on which the |ao_extra| is centered
 
 
.. option:: ao_extra_power
 
    Powers of x, y and z for each |ao_extra|
 
 
.. option:: ao_extra_coef
 
    Primitive coefficients, read from input. Those should not be used directly, as the MOs are expressed on the basis of **normalized** ao_extras.
 
 
.. option:: ao_extra_expo
 
    Exponents for each primitive of each |ao_extra|
 
 
.. option:: ao_extra_md5
 
    MD5 key, specific of the |ao_extra| basis
 
 
.. option:: ao_extra_cartesian
 
    If |true|, use |ao_extras| in Cartesian coordinates (6d,10f,...)
 
    Default: false
 
.. option:: ao_extra_normalized
 
    Use normalized basis functions
 
    Default: true
 
.. option:: primitives_normalized_extra
 
    Use normalized primitive functions
 
    Default: true
 
.. option:: ao_extra_expo_im
 
    imag part for Exponents for each primitive of each cGTOs |ao_extra|
 
 
.. option:: ao_extra_expo_pw
 
    plane wave part for each primitive GTOs |ao_extra|
 
 
.. option:: ao_extra_expo_phase
 
    phase shift for each primitive GTOs |ao_extra|
 
 
.. option:: ao_extra_one_e_dm
 
    reduced density matrix on the ao extra basis
 
 
.. option:: ao_extra_center
 
    shift with which the atoms are shifted to mimick p functions
 
 
 
Programs 
-------- 
 
 * :ref:`extra_basis` 
 * :ref:`fit_1s_basis` 
 * :ref:`write_extra_grid` 
 
Providers 
--------- 
 
.. c:var:: double ao_extra_center


    File : :file:`prov_fit_1s.irp.f`

    .. code:: fortran

        double precision	:: ao_extra_center	



    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`new_nucl_coord_1s_transp`

 
.. c:var:: double ao_extra_center_1s


    File : :file:`density_extra.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_extra_center_1s	(3,ao_extra_num)


    Original position of each extra AO

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_nucl`
       * :c:data:`ao_extra_num`
       * :c:data:`extra_nucl_coord`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_center_1s_dist`

 
.. c:var:: double ao_extra_center_1s_dist


    File : :file:`density_extra.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_extra_center_1s_dist	(ao_extra_num,ao_extra_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_center_1s`
       * :c:data:`ao_extra_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`e_pq_ao_extra`

 
.. c:var:: double ao_extra_coef_normalization_factor


    File : :file:`aos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_extra_coef_normalized	(ao_extra_num,ao_extra_prim_num_max)
        double precision, allocatable	:: ao_extra_coef_normalization_factor	(ao_extra_num)


    Coefficients including the |ao_extra| normalization

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_coef`
       * :c:data:`ao_extra_expo`
       * :c:data:`ao_extra_normalized`
       * :c:data:`ao_extra_num`
       * :c:data:`ao_extra_power`
       * :c:data:`ao_extra_prim_num`
       * :c:data:`ao_extra_prim_num_max`
       * :c:data:`primitives_normalized_extra`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_coef_normalized_ordered`
       * :c:data:`effective_ao_extra_dm`

 
.. c:var:: double ao_extra_coef_normalized


    File : :file:`aos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_extra_coef_normalized	(ao_extra_num,ao_extra_prim_num_max)
        double precision, allocatable	:: ao_extra_coef_normalization_factor	(ao_extra_num)


    Coefficients including the |ao_extra| normalization

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_coef`
       * :c:data:`ao_extra_expo`
       * :c:data:`ao_extra_normalized`
       * :c:data:`ao_extra_num`
       * :c:data:`ao_extra_power`
       * :c:data:`ao_extra_prim_num`
       * :c:data:`ao_extra_prim_num_max`
       * :c:data:`primitives_normalized_extra`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_coef_normalized_ordered`
       * :c:data:`effective_ao_extra_dm`

 
.. c:var:: double ao_extra_coef_normalized_ordered


    File : :file:`aos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_extra_coef_normalized_ordered	(ao_extra_num,ao_extra_prim_num_max)
        double precision, allocatable	:: ao_extra_expo_ordered	(ao_extra_num,ao_extra_prim_num_max)


    Sorted primitives to accelerate 4 index |MO| transformation

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_coef_normalized`
       * :c:data:`ao_extra_expo`
       * :c:data:`ao_extra_num`
       * :c:data:`ao_extra_prim_num`
       * :c:data:`ao_extra_prim_num_max`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_coef_normalized_ordered_transp`
       * :c:data:`ao_extra_expo_ordered_transp`

 
.. c:var:: double ao_extra_coef_normalized_ordered_transp


    File : :file:`aos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_extra_coef_normalized_ordered_transp	(ao_extra_prim_num_max,ao_extra_num)


    Transposed :c:data:`ao_extra_coef_normalized_ordered`

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_coef_normalized_ordered`
       * :c:data:`ao_extra_num`
       * :c:data:`ao_extra_prim_num_max`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_coef_normalized_ordered_transp_per_nucl`

 
.. c:var:: double ao_extra_coef_normalized_ordered_transp_per_nucl


    File : :file:`aos_transp.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_extra_coef_normalized_ordered_transp_per_nucl	(ao_extra_prim_num_max,N_ao_extras_max,nucl_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_coef_normalized_ordered_transp`
       * :c:data:`ao_extra_prim_num`
       * :c:data:`ao_extra_prim_num_max`
       * :c:data:`nucl_ao_extras_transposed`
       * :c:data:`nucl_n_ao_extras`
       * :c:data:`nucl_num`


 
.. c:var:: double ao_extra_expo_ordered


    File : :file:`aos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_extra_coef_normalized_ordered	(ao_extra_num,ao_extra_prim_num_max)
        double precision, allocatable	:: ao_extra_expo_ordered	(ao_extra_num,ao_extra_prim_num_max)


    Sorted primitives to accelerate 4 index |MO| transformation

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_coef_normalized`
       * :c:data:`ao_extra_expo`
       * :c:data:`ao_extra_num`
       * :c:data:`ao_extra_prim_num`
       * :c:data:`ao_extra_prim_num_max`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_coef_normalized_ordered_transp`
       * :c:data:`ao_extra_expo_ordered_transp`

 
.. c:var:: double ao_extra_expo_ordered_transp


    File : :file:`aos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_extra_expo_ordered_transp	(ao_extra_prim_num_max,ao_extra_num)


    Transposed :c:data:`ao_extra_expo_ordered`

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_coef_normalized_ordered`
       * :c:data:`ao_extra_num`
       * :c:data:`ao_extra_prim_num_max`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_expo_ordered_transp_per_nucl`

 
.. c:var:: double ao_extra_expo_ordered_transp_per_nucl


    File : :file:`aos_transp.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_extra_expo_ordered_transp_per_nucl	(ao_extra_prim_num_max,N_ao_extras_max,nucl_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_expo_ordered_transp`
       * :c:data:`ao_extra_prim_num`
       * :c:data:`ao_extra_prim_num_max`
       * :c:data:`nucl_ao_extras_transposed`
       * :c:data:`nucl_n_ao_extras`
       * :c:data:`nucl_num`


 
.. c:var:: int ao_extra_first_of_shell


    File : :file:`aos.irp.f`

    .. code:: fortran

        integer, allocatable	:: ao_extra_first_of_shell	(shell_num)


    Index of the shell to which the ao_extra corresponds

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`shell_ang_mom`
       * :c:data:`shell_num`


 
.. c:var:: int ao_extra_l


    File : :file:`aos.irp.f`

    .. code:: fortran

        integer, allocatable	:: ao_extra_l	(ao_extra_num)
        integer	:: ao_extra_l_max	
        character*(128), allocatable	:: ao_extra_l_char	(ao_extra_num)


    :math:`l` value of the |ao_extra|: :math`a+b+c` in :math:`x^a y^b z^c`

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_num`
       * :c:data:`ao_extra_power`
       * :c:data:`l_to_character`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_l_char_space`
       * :c:data:`nucl_list_shell_ao_extras`

 
.. c:var:: char ao_extra_l_char


    File : :file:`aos.irp.f`

    .. code:: fortran

        integer, allocatable	:: ao_extra_l	(ao_extra_num)
        integer	:: ao_extra_l_max	
        character*(128), allocatable	:: ao_extra_l_char	(ao_extra_num)


    :math:`l` value of the |ao_extra|: :math`a+b+c` in :math:`x^a y^b z^c`

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_num`
       * :c:data:`ao_extra_power`
       * :c:data:`l_to_character`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_l_char_space`
       * :c:data:`nucl_list_shell_ao_extras`

 
.. c:var:: char ao_extra_l_char_space


    File : :file:`aos.irp.f`

    .. code:: fortran

        character*(4), allocatable	:: ao_extra_l_char_space	(ao_extra_num)


    Converts an l value to a string

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_l`
       * :c:data:`ao_extra_num`
       * :c:data:`ao_extra_power`


 
.. c:var:: int ao_extra_l_max


    File : :file:`aos.irp.f`

    .. code:: fortran

        integer, allocatable	:: ao_extra_l	(ao_extra_num)
        integer	:: ao_extra_l_max	
        character*(128), allocatable	:: ao_extra_l_char	(ao_extra_num)


    :math:`l` value of the |ao_extra|: :math`a+b+c` in :math:`x^a y^b z^c`

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_num`
       * :c:data:`ao_extra_power`
       * :c:data:`l_to_character`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_l_char_space`
       * :c:data:`nucl_list_shell_ao_extras`

 
.. c:var:: double ao_extra_one_e_dm_at_extra_r


    File : :file:`aos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_extra_one_e_dm_at_extra_r	(n_points_extra_final_grid)


    ao_extra_one_e_dm_at_extra_r(i) = extra density on the extra grid points

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_num`
       * :c:data:`ao_extra_one_e_dm`
       * :c:data:`final_grid_points_extra`
       * :c:data:`n_points_extra_final_grid`


 
.. c:var:: int ao_extra_power_ordered_transp_per_nucl


    File : :file:`aos_transp.irp.f`

    .. code:: fortran

        integer, allocatable	:: ao_extra_power_ordered_transp_per_nucl	(3,N_ao_extras_max,nucl_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_power`
       * :c:data:`nucl_ao_extras_transposed`
       * :c:data:`nucl_n_ao_extras`
       * :c:data:`nucl_num`


 
.. c:var:: int ao_extra_prim_num_max


    File : :file:`aos.irp.f`

    .. code:: fortran

        integer	:: ao_extra_prim_num_max	


    Max number of primitives.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_prim_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_coef`
       * :c:data:`ao_extra_coef_normalized`
       * :c:data:`ao_extra_coef_normalized_ordered`
       * :c:data:`ao_extra_coef_normalized_ordered_transp`
       * :c:data:`ao_extra_coef_normalized_ordered_transp_per_nucl`
       * :c:data:`ao_extra_expo`
       * :c:data:`ao_extra_expo_im`
       * :c:data:`ao_extra_expo_ordered_transp`
       * :c:data:`ao_extra_expo_ordered_transp_per_nucl`
       * :c:data:`ao_extra_expo_phase`
       * :c:data:`ao_extra_expo_pw`

 
.. c:var:: int ao_extra_shell


    File : :file:`aos.irp.f`

    .. code:: fortran

        integer, allocatable	:: ao_extra_shell	(ao_extra_num)


    Index of the shell to which the ao_extra corresponds

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_num`
       * :c:data:`shell_ang_mom`
       * :c:data:`shell_num`


 
.. c:var:: double e_pq_ao_extra


    File : :file:`density_extra.irp.f`

    .. code:: fortran

        double precision, allocatable	:: e_pq_ao_extra	(ao_extra_num,ao_extra_num)


    E_pq_ao_extra = exp(-alpha_p alpha_q/gamma_pq (Q_p - Q_q)^2)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_center_1s_dist`
       * :c:data:`ao_extra_expo`
       * :c:data:`ao_extra_num`
       * :c:data:`gamma_pq_ao_extra`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`effective_ao_extra_dm`

 
.. c:var:: double effective_ao_extra_dm


    File : :file:`density_extra.irp.f`

    .. code:: fortran

        double precision, allocatable	:: effective_ao_extra_dm	(ao_extra_num,ao_extra_num)


    effective density matrix : rho_pq x N_p x N_q x E_pq x (pi/gamma_pq)^3/2
    
    where rho_pq is the usual density matrix
    
    N_p and N_q are the normalization factors associated with the two Gaussians
    
    E_pq is the prefactor resulting from the Gaussian product
    
    gamma_pq = gamma_p + gamm_q

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_coef_normalized`
       * :c:data:`ao_extra_num`
       * :c:data:`ao_extra_one_e_dm`
       * :c:data:`e_pq_ao_extra`
       * :c:data:`gamma_pq_ao_extra`


 
.. c:var:: int extra_fictious_nucl


    File : :file:`prov_fit_1s.irp.f`

    .. code:: fortran

        integer	:: extra_fictious_nucl	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_2p_func_orig`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`new_nucl_coord_1s_transp`
       * :c:data:`new_nucl_label_1s`

 
.. c:var:: int extra_nucl_real_fictious_list_prov


    File : :file:`prov_fit_1s.irp.f`

    .. code:: fortran

        double precision, allocatable	:: new_nucl_coord_1s_transp	(3,new_nucl_num)
        double precision, allocatable	:: new_nucl_charge_1s	(new_nucl_num)
        integer, allocatable	:: extra_nucl_real_fictious_list_prov	(extra_fictious_nucl)


    the real atoms are located in the first nucl_num entries
    
    then the fictious atoms are located after

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_center`
       * :c:data:`ao_l`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`extra_fictious_nucl`
       * :c:data:`new_nucl_num`
       * :c:data:`nucl_aos_transposed`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_coord_transp`
       * :c:data:`nucl_n_aos`
       * :c:data:`nucl_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`new_nucl_coord_1s`

 
.. c:var:: double gamma_pq_ao_extra


    File : :file:`density_extra.irp.f`

    .. code:: fortran

        double precision, allocatable	:: gamma_pq_ao_extra	(ao_extra_num,ao_extra_num)
        double precision, allocatable	:: inv_gamma_pq_ao_extra	(ao_extra_num,ao_extra_num)
        double precision, allocatable	:: inv_pi_gamma_pq_3_2_ao_extra	(ao_extra_num,ao_extra_num)
        double precision, allocatable	:: sqrt_gamma_pq_ao_extra	(ao_extra_num,ao_extra_num)


    gamma_pq_ao_extra = gamma_p + gamma_q
    
    inv_gamma_pq_ao_extra = 1/(gamma_p + gamma_q)
    
    sqrt_gamma_pq_ao_extra = sqrt(gamma_p + gamma_q)
    
    WARNING :: VALID ONLY IN THE CASE OF A PURELY 1S BASIS

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_expo`
       * :c:data:`ao_extra_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`e_pq_ao_extra`
       * :c:data:`effective_ao_extra_dm`

 
.. c:var:: double inv_gamma_pq_ao_extra


    File : :file:`density_extra.irp.f`

    .. code:: fortran

        double precision, allocatable	:: gamma_pq_ao_extra	(ao_extra_num,ao_extra_num)
        double precision, allocatable	:: inv_gamma_pq_ao_extra	(ao_extra_num,ao_extra_num)
        double precision, allocatable	:: inv_pi_gamma_pq_3_2_ao_extra	(ao_extra_num,ao_extra_num)
        double precision, allocatable	:: sqrt_gamma_pq_ao_extra	(ao_extra_num,ao_extra_num)


    gamma_pq_ao_extra = gamma_p + gamma_q
    
    inv_gamma_pq_ao_extra = 1/(gamma_p + gamma_q)
    
    sqrt_gamma_pq_ao_extra = sqrt(gamma_p + gamma_q)
    
    WARNING :: VALID ONLY IN THE CASE OF A PURELY 1S BASIS

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_expo`
       * :c:data:`ao_extra_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`e_pq_ao_extra`
       * :c:data:`effective_ao_extra_dm`

 
.. c:var:: double inv_pi_gamma_pq_3_2_ao_extra


    File : :file:`density_extra.irp.f`

    .. code:: fortran

        double precision, allocatable	:: gamma_pq_ao_extra	(ao_extra_num,ao_extra_num)
        double precision, allocatable	:: inv_gamma_pq_ao_extra	(ao_extra_num,ao_extra_num)
        double precision, allocatable	:: inv_pi_gamma_pq_3_2_ao_extra	(ao_extra_num,ao_extra_num)
        double precision, allocatable	:: sqrt_gamma_pq_ao_extra	(ao_extra_num,ao_extra_num)


    gamma_pq_ao_extra = gamma_p + gamma_q
    
    inv_gamma_pq_ao_extra = 1/(gamma_p + gamma_q)
    
    sqrt_gamma_pq_ao_extra = sqrt(gamma_p + gamma_q)
    
    WARNING :: VALID ONLY IN THE CASE OF A PURELY 1S BASIS

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_expo`
       * :c:data:`ao_extra_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`e_pq_ao_extra`
       * :c:data:`effective_ao_extra_dm`

 
.. c:var:: int list_2p_functions


    File : :file:`prov_fit_1s.irp.f`

    .. code:: fortran

        integer, allocatable	:: list_2p_functions	(n_2p_func_orig)


    list of 2p functions in the original basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_l`
       * :c:data:`ao_num`
       * :c:data:`n_2p_func_orig`


 
.. c:var:: int list_fict_nucl


    File : :file:`prov_fit_1s.irp.f`

    .. code:: fortran

        character*(32), allocatable	:: new_nucl_label_1s	(new_nucl_num)
        integer, allocatable	:: list_real_nucl	(nucl_num)
        integer, allocatable	:: list_fict_nucl	(extra_fictious_nucl)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`extra_fictious_nucl`
       * :c:data:`new_nucl_num`
       * :c:data:`nucl_label`
       * :c:data:`nucl_num`


 
.. c:var:: int list_real_nucl


    File : :file:`prov_fit_1s.irp.f`

    .. code:: fortran

        character*(32), allocatable	:: new_nucl_label_1s	(new_nucl_num)
        integer, allocatable	:: list_real_nucl	(nucl_num)
        integer, allocatable	:: list_fict_nucl	(extra_fictious_nucl)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`extra_fictious_nucl`
       * :c:data:`new_nucl_num`
       * :c:data:`nucl_label`
       * :c:data:`nucl_num`


 
.. c:var:: bool lmax_too_big


    File : :file:`prov_fit_1s.irp.f`

    .. code:: fortran

        logical	:: lmax_too_big	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_l`


 
.. c:var:: int n_2p_func_orig


    File : :file:`prov_fit_1s.irp.f`

    .. code:: fortran

        integer	:: n_2p_func_orig	
        integer	:: n_2p_func_tot	


    n_2p_func_orig :: number of 2p functions in the original basis
    
    n_2p_func_tot :: total number of p functions in the fitted basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_l`
       * :c:data:`ao_num`
       * :c:data:`ao_prim_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`extra_fictious_nucl`
       * :c:data:`list_2p_functions`
       * :c:data:`new_nucl_num`

 
.. c:var:: int n_2p_func_tot


    File : :file:`prov_fit_1s.irp.f`

    .. code:: fortran

        integer	:: n_2p_func_orig	
        integer	:: n_2p_func_tot	


    n_2p_func_orig :: number of 2p functions in the original basis
    
    n_2p_func_tot :: total number of p functions in the fitted basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_l`
       * :c:data:`ao_num`
       * :c:data:`ao_prim_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`extra_fictious_nucl`
       * :c:data:`list_2p_functions`
       * :c:data:`new_nucl_num`

 
.. c:var:: int n_ao_extras_max


    File : :file:`aos.irp.f`

    .. code:: fortran

        integer, allocatable	:: nucl_n_ao_extras	(extra_nucl_num)
        integer	:: n_ao_extras_max	


    Number of |ao_extras| per atom

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_nucl`
       * :c:data:`ao_extra_num`
       * :c:data:`extra_nucl_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_coef_normalized_ordered_transp_per_nucl`
       * :c:data:`ao_extra_expo_ordered_transp_per_nucl`
       * :c:data:`ao_extra_power_ordered_transp_per_nucl`
       * :c:data:`nucl_ao_extras`
       * :c:data:`nucl_ao_extras_transposed`
       * :c:data:`nucl_list_shell_ao_extras`

 
.. c:var:: int n_func_tot


    File : :file:`prov_fit_1s.irp.f`

    .. code:: fortran

        integer	:: n_func_tot	


    n_func_tot :: total number of functions in the fitted basis set
    
    returned in an uncontracted way

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_l`
       * :c:data:`ao_num`
       * :c:data:`ao_prim_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`new_ao_coef_1s`
       * :c:data:`new_ao_power_1s`
       * :c:data:`new_ao_prim_num_1s`
       * :c:data:`new_nucl_n_aos`

 
.. c:var:: int n_prim_tot_orig


    File : :file:`prov_fit_1s.irp.f`

    .. code:: fortran

        integer	:: n_prim_tot_orig	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_prim_num`


 
.. c:var:: int n_pt_max_extra_basis_i_x


    File : :file:`dimensions_integrals.irp.f`

    .. code:: fortran

        integer	:: n_pt_max_extra_basis_integrals	
        integer	:: n_pt_max_extra_basis_i_x	


    Number of points used in the numerical integrations.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_power`


 
.. c:var:: int n_pt_max_extra_basis_integrals


    File : :file:`dimensions_integrals.irp.f`

    .. code:: fortran

        integer	:: n_pt_max_extra_basis_integrals	
        integer	:: n_pt_max_extra_basis_i_x	


    Number of points used in the numerical integrations.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_power`


 
.. c:var:: double new_ao_coef_1s


    File : :file:`prov_fit_1s.irp.f`

    .. code:: fortran

        double precision, allocatable	:: new_ao_coef_1s	(n_func_tot)


    Primitive coefficients, read from input. Those should not be used directly, as the MOs are expressed on the basis of **normalized** AOs.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_func_tot`


 
.. c:var:: double new_ao_expo_1s


    File : :file:`prov_fit_1s.irp.f`

    .. code:: fortran

        integer, allocatable	:: new_nucl_n_aos	(new_nucl_num)
        integer, allocatable	:: new_nucl_aos_transposed	(new_n_AOs_max,new_nucl_num)
        double precision, allocatable	:: new_ao_expo_1s	(n_func_tot)
        integer, allocatable	:: new_ao_nucl_1s	(n_func_tot)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_expo`
       * :c:data:`ao_l`
       * :c:data:`ao_prim_num`
       * :c:data:`n_func_tot`
       * :c:data:`new_n_aos_max`
       * :c:data:`new_nucl_num`
       * :c:data:`nucl_aos_transposed`
       * :c:data:`nucl_n_aos`
       * :c:data:`nucl_num`


 
.. c:var:: int new_ao_nucl_1s


    File : :file:`prov_fit_1s.irp.f`

    .. code:: fortran

        integer, allocatable	:: new_nucl_n_aos	(new_nucl_num)
        integer, allocatable	:: new_nucl_aos_transposed	(new_n_AOs_max,new_nucl_num)
        double precision, allocatable	:: new_ao_expo_1s	(n_func_tot)
        integer, allocatable	:: new_ao_nucl_1s	(n_func_tot)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_expo`
       * :c:data:`ao_l`
       * :c:data:`ao_prim_num`
       * :c:data:`n_func_tot`
       * :c:data:`new_n_aos_max`
       * :c:data:`new_nucl_num`
       * :c:data:`nucl_aos_transposed`
       * :c:data:`nucl_n_aos`
       * :c:data:`nucl_num`


 
.. c:var:: int new_ao_power_1s


    File : :file:`prov_fit_1s.irp.f`

    .. code:: fortran

        integer, allocatable	:: new_ao_power_1s	(n_func_tot,3)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_func_tot`


 
.. c:var:: int new_ao_prim_num_1s


    File : :file:`prov_fit_1s.irp.f`

    .. code:: fortran

        integer, allocatable	:: new_ao_prim_num_1s	(n_func_tot)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_func_tot`


 
.. c:var:: int new_n_aos_max


    File : :file:`prov_fit_1s.irp.f`

    .. code:: fortran

        integer	:: new_n_aos_max	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_prim_num_max`
       * :c:data:`nucl_n_aos`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`new_nucl_n_aos`

 
.. c:var:: int new_nucl_aos_transposed


    File : :file:`prov_fit_1s.irp.f`

    .. code:: fortran

        integer, allocatable	:: new_nucl_n_aos	(new_nucl_num)
        integer, allocatable	:: new_nucl_aos_transposed	(new_n_AOs_max,new_nucl_num)
        double precision, allocatable	:: new_ao_expo_1s	(n_func_tot)
        integer, allocatable	:: new_ao_nucl_1s	(n_func_tot)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_expo`
       * :c:data:`ao_l`
       * :c:data:`ao_prim_num`
       * :c:data:`n_func_tot`
       * :c:data:`new_n_aos_max`
       * :c:data:`new_nucl_num`
       * :c:data:`nucl_aos_transposed`
       * :c:data:`nucl_n_aos`
       * :c:data:`nucl_num`


 
.. c:var:: double new_nucl_charge_1s


    File : :file:`prov_fit_1s.irp.f`

    .. code:: fortran

        double precision, allocatable	:: new_nucl_coord_1s_transp	(3,new_nucl_num)
        double precision, allocatable	:: new_nucl_charge_1s	(new_nucl_num)
        integer, allocatable	:: extra_nucl_real_fictious_list_prov	(extra_fictious_nucl)


    the real atoms are located in the first nucl_num entries
    
    then the fictious atoms are located after

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_center`
       * :c:data:`ao_l`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`extra_fictious_nucl`
       * :c:data:`new_nucl_num`
       * :c:data:`nucl_aos_transposed`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_coord_transp`
       * :c:data:`nucl_n_aos`
       * :c:data:`nucl_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`new_nucl_coord_1s`

 
.. c:var:: double new_nucl_coord_1s


    File : :file:`prov_fit_1s.irp.f`

    .. code:: fortran

        double precision, allocatable	:: new_nucl_coord_1s	(new_nucl_num,3)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`new_nucl_coord_1s_transp`
       * :c:data:`new_nucl_num`


 
.. c:var:: double new_nucl_coord_1s_transp


    File : :file:`prov_fit_1s.irp.f`

    .. code:: fortran

        double precision, allocatable	:: new_nucl_coord_1s_transp	(3,new_nucl_num)
        double precision, allocatable	:: new_nucl_charge_1s	(new_nucl_num)
        integer, allocatable	:: extra_nucl_real_fictious_list_prov	(extra_fictious_nucl)


    the real atoms are located in the first nucl_num entries
    
    then the fictious atoms are located after

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_center`
       * :c:data:`ao_l`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`extra_fictious_nucl`
       * :c:data:`new_nucl_num`
       * :c:data:`nucl_aos_transposed`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_coord_transp`
       * :c:data:`nucl_n_aos`
       * :c:data:`nucl_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`new_nucl_coord_1s`

 
.. c:var:: char new_nucl_label_1s


    File : :file:`prov_fit_1s.irp.f`

    .. code:: fortran

        character*(32), allocatable	:: new_nucl_label_1s	(new_nucl_num)
        integer, allocatable	:: list_real_nucl	(nucl_num)
        integer, allocatable	:: list_fict_nucl	(extra_fictious_nucl)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`extra_fictious_nucl`
       * :c:data:`new_nucl_num`
       * :c:data:`nucl_label`
       * :c:data:`nucl_num`


 
.. c:var:: int new_nucl_n_aos


    File : :file:`prov_fit_1s.irp.f`

    .. code:: fortran

        integer, allocatable	:: new_nucl_n_aos	(new_nucl_num)
        integer, allocatable	:: new_nucl_aos_transposed	(new_n_AOs_max,new_nucl_num)
        double precision, allocatable	:: new_ao_expo_1s	(n_func_tot)
        integer, allocatable	:: new_ao_nucl_1s	(n_func_tot)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_expo`
       * :c:data:`ao_l`
       * :c:data:`ao_prim_num`
       * :c:data:`n_func_tot`
       * :c:data:`new_n_aos_max`
       * :c:data:`new_nucl_num`
       * :c:data:`nucl_aos_transposed`
       * :c:data:`nucl_n_aos`
       * :c:data:`nucl_num`


 
.. c:var:: int new_nucl_num


    File : :file:`prov_fit_1s.irp.f`

    .. code:: fortran

        integer	:: new_nucl_num	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_2p_func_orig`
       * :c:data:`nucl_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`new_nucl_coord_1s`
       * :c:data:`new_nucl_coord_1s_transp`
       * :c:data:`new_nucl_label_1s`
       * :c:data:`new_nucl_n_aos`

 
.. c:var:: int nucl_ao_extras


    File : :file:`aos.irp.f`

    .. code:: fortran

        integer, allocatable	:: nucl_ao_extras	(extra_nucl_num,N_ao_extras_max)


    List of |ao_extras| centered on each atom

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_nucl`
       * :c:data:`ao_extra_num`
       * :c:data:`extra_nucl_num`
       * :c:data:`nucl_n_ao_extras`
       * :c:data:`nucl_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_list_shell_ao_extras`

 
.. c:var:: int nucl_ao_extras_transposed


    File : :file:`aos_transp.irp.f`

    .. code:: fortran

        integer, allocatable	:: nucl_ao_extras_transposed	(N_ao_extras_max,nucl_num)


    List of ao_extras attached on each atom

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_nucl`
       * :c:data:`ao_extra_num`
       * :c:data:`nucl_n_ao_extras`
       * :c:data:`nucl_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_coef_normalized_ordered_transp_per_nucl`
       * :c:data:`ao_extra_expo_ordered_transp_per_nucl`
       * :c:data:`ao_extra_power_ordered_transp_per_nucl`

 
.. c:var:: int nucl_list_shell_ao_extras


    File : :file:`aos.irp.f`

    .. code:: fortran

        integer, allocatable	:: nucl_list_shell_ao_extras	(extra_nucl_num,N_ao_extras_max)
        integer, allocatable	:: nucl_num_shell_ao_extras	(nucl_num)


    Index of the shell type |ao_extras| and of the corresponding |ao_extras|
    By convention, for p,d,f and g |ao_extras|, we take the index
    of the |ao_extra| with the the corresponding power in the x axis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_l`
       * :c:data:`ao_extra_power`
       * :c:data:`extra_nucl_num`
       * :c:data:`nucl_ao_extras`
       * :c:data:`nucl_n_ao_extras`
       * :c:data:`nucl_num`


 
.. c:var:: int nucl_n_ao_extras


    File : :file:`aos.irp.f`

    .. code:: fortran

        integer, allocatable	:: nucl_n_ao_extras	(extra_nucl_num)
        integer	:: n_ao_extras_max	


    Number of |ao_extras| per atom

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_nucl`
       * :c:data:`ao_extra_num`
       * :c:data:`extra_nucl_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_coef_normalized_ordered_transp_per_nucl`
       * :c:data:`ao_extra_expo_ordered_transp_per_nucl`
       * :c:data:`ao_extra_power_ordered_transp_per_nucl`
       * :c:data:`nucl_ao_extras`
       * :c:data:`nucl_ao_extras_transposed`
       * :c:data:`nucl_list_shell_ao_extras`

 
.. c:var:: int nucl_num_shell_ao_extras


    File : :file:`aos.irp.f`

    .. code:: fortran

        integer, allocatable	:: nucl_list_shell_ao_extras	(extra_nucl_num,N_ao_extras_max)
        integer, allocatable	:: nucl_num_shell_ao_extras	(nucl_num)


    Index of the shell type |ao_extras| and of the corresponding |ao_extras|
    By convention, for p,d,f and g |ao_extras|, we take the index
    of the |ao_extra| with the the corresponding power in the x axis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_l`
       * :c:data:`ao_extra_power`
       * :c:data:`extra_nucl_num`
       * :c:data:`nucl_ao_extras`
       * :c:data:`nucl_n_ao_extras`
       * :c:data:`nucl_num`


 
.. c:var:: double sqrt_gamma_pq_ao_extra


    File : :file:`density_extra.irp.f`

    .. code:: fortran

        double precision, allocatable	:: gamma_pq_ao_extra	(ao_extra_num,ao_extra_num)
        double precision, allocatable	:: inv_gamma_pq_ao_extra	(ao_extra_num,ao_extra_num)
        double precision, allocatable	:: inv_pi_gamma_pq_3_2_ao_extra	(ao_extra_num,ao_extra_num)
        double precision, allocatable	:: sqrt_gamma_pq_ao_extra	(ao_extra_num,ao_extra_num)


    gamma_pq_ao_extra = gamma_p + gamma_q
    
    inv_gamma_pq_ao_extra = 1/(gamma_p + gamma_q)
    
    sqrt_gamma_pq_ao_extra = sqrt(gamma_p + gamma_q)
    
    WARNING :: VALID ONLY IN THE CASE OF A PURELY 1S BASIS

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_expo`
       * :c:data:`ao_extra_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`e_pq_ao_extra`
       * :c:data:`effective_ao_extra_dm`

 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: void ao_extra_power_index(void)


    File : :file:`aos.irp.f`

    .. code:: fortran

        integer function ao_extra_power_index(nx,ny,nz)


    Unique index given to a triplet of powers:
    
    :math:`\frac{1}{2} (l-n_x) (l-n_x+1) + n_z + 1`

 
.. c:function:: void ao_extra_value(void)


    File : :file:`aos.irp.f`

    .. code:: fortran

        double precision function ao_extra_value(i, r)


    Returns the value of the i-th ao belonging the EXTRA BASIS at point $\textbf{r}$

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_coef_normalized_ordered_transp`
       * :c:data:`ao_extra_expo_ordered_transp`
       * :c:data:`ao_extra_nucl`
       * :c:data:`ao_extra_power`
       * :c:data:`ao_extra_prim_num`
       * :c:data:`extra_nucl_coord`

 
.. c:function:: void extra_density_at_r(void)


    File : :file:`aos.irp.f`

    .. code:: fortran

        double precision function extra_density_at_r(r)


    
    input  : r == r(1) = x and so on
    
    output : density corresponding to the extra system at r
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_num`
       * :c:data:`ao_extra_one_e_dm`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_all_aos_extra_at_r`

 
.. c:function:: void give_all_aos_extra_at_r(void)


    File : :file:`aos.irp.f`

    .. code:: fortran

        subroutine give_all_aos_extra_at_r(r, tmp_array)


    
    input  : r == r(1) = x and so on
    
    output : tmp_array(i) = EXTRA aos(i) evaluated in $\textbf{r}$
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_extra_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`extra_density_at_r`

 
.. c:function:: void r_phi_theta_psi_matrix(void)


    File : :file:`transform_basis.irp.f`

    .. code:: fortran

        subroutine r_phi_theta_psi_matrix(phi,theta,psi,r_mat)


    routine that creates the rotation matrix corresponding to phi,theta,psi
    
    according to conventions in MDFT code

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`rotate_nuclei`

 
.. c:function:: void return_xyz_int(void)


    File : :file:`prov_fit_1s.irp.f`

    .. code:: fortran

        integer function return_xyz_int(power)



 
.. c:function:: void rotate_nuclei(void)


    File : :file:`transform_basis.irp.f`

    .. code:: fortran

        subroutine rotate_nuclei(phi,theta,psi,nucl_centers,n_nucl,nucl_centers_after)


    routine that rotates a set of nuclei according to three axis corresponding to angles phi, theta, psi.

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_ab_prod`
       * :c:func:`r_phi_theta_psi_matrix`

 
.. c:function:: void routine(void)


    File : :file:`write_extra_grid.irp.f`

    .. code:: fortran

        subroutine routine



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points_extra`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`write_extra_grid_in_ezfio`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_set_becke_numerical_grid_io_extra_grid`

