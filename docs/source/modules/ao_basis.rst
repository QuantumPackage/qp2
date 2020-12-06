.. _module_ao_basis: 
 
.. program:: ao_basis 
 
.. default-role:: option 
 
========
ao_basis
========

This module describes the atomic orbitals basis set.

An |AO| :math:`\chi` centered on nucleus A is represented as:

.. math::

   \chi_i({\bf r}) = (x-X_A)^a (y-Y_A)^b (z-Z_A)^c \sum_k c_{ki} e^{-\gamma_{ki} |{\bf r} - {\bf R}_A|^2}


The |AO| coefficients are normalized as:

.. math::

  {\tilde c}_{ki} = \frac{c_{ki}}{ \int \left( (x-X_A)^a (y-Y_A)^b (z-Z_A)^c  e^{-\gamma_{ki} |{\bf r} - {\bf R}_A|^2} \right)^2 dr}


.. warning::

  `ao_coef` contains the |AO| coefficients given in input. These do not
  include the normalization constant of the |AO|. The `ao_coef_normalized`
  provider includes this normalization factor.


The |AOs| are also sorted by increasing exponent to accelerate the calculation of
the two electron integrals.



 
 
 
EZFIO parameters 
---------------- 
 
.. option:: ao_basis
 
    Name of the |AO| basis set
 
 
.. option:: ao_num
 
    Number of |AOs|
 
 
.. option:: ao_prim_num
 
    Number of primitives per |AO|
 
 
.. option:: ao_prim_num_max
 
    Maximum number of primitives
 
    Default: =maxval(ao_basis.ao_prim_num)
 
.. option:: ao_nucl
 
    Index of the nucleus on which the |AO| is centered
 
 
.. option:: ao_power
 
    Powers of x, y and z for each |AO|
 
 
.. option:: ao_coef
 
    Primitive coefficients, read from input. Those should not be used directly, as the MOs are expressed on the basis of **normalized** AOs.
 
 
.. option:: ao_expo
 
    Exponents for each primitive of each |AO|
 
 
.. option:: ao_md5
 
    MD5 key, specific of the |AO| basis
 
 
.. option:: ao_cartesian
 
    If |true|, use |AOs| in Cartesian coordinates (6d,10f,...)
 
    Default: false
 
 
Providers 
--------- 
 
.. c:var:: ao_coef_normalization_factor


    File : :file:`ao_basis/aos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_coef_normalized	(ao_num,ao_prim_num_max)
        double precision, allocatable	:: ao_coef_normalization_factor	(ao_num)


    Coefficients including the |AO| normalization

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef`
       * :c:data:`ao_expo`
       * :c:data:`ao_normalized`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`ao_prim_num_max`
       * :c:data:`primitives_normalized`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalization_libint_factor`
       * :c:data:`ao_coef_normalized_ordered`

 
.. c:var:: ao_coef_normalization_libint_factor


    File : :file:`ao_basis/aos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_coef_normalization_libint_factor	(ao_num)


    |AO| normalization for interfacing with libint

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized`
       * :c:data:`ao_expo`
       * :c:data:`ao_l`
       * :c:data:`ao_num`
       * :c:data:`ao_prim_num`


 
.. c:var:: ao_coef_normalized


    File : :file:`ao_basis/aos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_coef_normalized	(ao_num,ao_prim_num_max)
        double precision, allocatable	:: ao_coef_normalization_factor	(ao_num)


    Coefficients including the |AO| normalization

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef`
       * :c:data:`ao_expo`
       * :c:data:`ao_normalized`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`ao_prim_num_max`
       * :c:data:`primitives_normalized`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalization_libint_factor`
       * :c:data:`ao_coef_normalized_ordered`

 
.. c:var:: ao_coef_normalized_ordered


    File : :file:`ao_basis/aos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_coef_normalized_ordered	(ao_num,ao_prim_num_max)
        double precision, allocatable	:: ao_expo_ordered	(ao_num,ao_prim_num_max)


    Sorted primitives to accelerate 4 index |MO| transformation

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized`
       * :c:data:`ao_expo`
       * :c:data:`ao_num`
       * :c:data:`ao_prim_num`
       * :c:data:`ao_prim_num_max`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_expo_ordered_transp`

 
.. c:var:: ao_coef_normalized_ordered_transp


    File : :file:`ao_basis/aos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_coef_normalized_ordered_transp	(ao_prim_num_max,ao_num)


    Transposed :c:data:`ao_coef_normalized_ordered`

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered`
       * :c:data:`ao_num`
       * :c:data:`ao_prim_num_max`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp_per_nucl`
       * :c:data:`ao_deriv2_x`
       * :c:data:`ao_deriv_1_x`
       * :c:data:`ao_dipole_x`
       * :c:data:`ao_integrals_n_e`
       * :c:data:`ao_integrals_n_e_per_atom`
       * :c:data:`ao_overlap`
       * :c:data:`ao_overlap_abs`
       * :c:data:`ao_pseudo_integrals_local`
       * :c:data:`ao_pseudo_integrals_non_local`
       * :c:data:`ao_spread_x`
       * :c:data:`ao_two_e_integral_schwartz`
       * :c:data:`ao_two_e_integrals_in_map`

 
.. c:var:: ao_coef_normalized_ordered_transp_per_nucl


    File : :file:`ao_basis/aos_transp.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_coef_normalized_ordered_transp_per_nucl	(ao_prim_num_max,N_AOs_max,nucl_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_prim_num`
       * :c:data:`ao_prim_num_max`
       * :c:data:`nucl_aos_transposed`
       * :c:data:`nucl_n_aos`
       * :c:data:`nucl_num`


 
.. c:var:: ao_expo_ordered


    File : :file:`ao_basis/aos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_coef_normalized_ordered	(ao_num,ao_prim_num_max)
        double precision, allocatable	:: ao_expo_ordered	(ao_num,ao_prim_num_max)


    Sorted primitives to accelerate 4 index |MO| transformation

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized`
       * :c:data:`ao_expo`
       * :c:data:`ao_num`
       * :c:data:`ao_prim_num`
       * :c:data:`ao_prim_num_max`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_expo_ordered_transp`

 
.. c:var:: ao_expo_ordered_transp


    File : :file:`ao_basis/aos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_expo_ordered_transp	(ao_prim_num_max,ao_num)


    Transposed :c:data:`ao_expo_ordered`

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered`
       * :c:data:`ao_num`
       * :c:data:`ao_prim_num_max`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_deriv2_x`
       * :c:data:`ao_deriv_1_x`
       * :c:data:`ao_dipole_x`
       * :c:data:`ao_expo_ordered_transp_per_nucl`
       * :c:data:`ao_integrals_n_e`
       * :c:data:`ao_integrals_n_e_per_atom`
       * :c:data:`ao_overlap`
       * :c:data:`ao_overlap_abs`
       * :c:data:`ao_pseudo_integrals_local`
       * :c:data:`ao_pseudo_integrals_non_local`
       * :c:data:`ao_spread_x`
       * :c:data:`ao_two_e_integral_schwartz`
       * :c:data:`ao_two_e_integrals_in_map`

 
.. c:var:: ao_expo_ordered_transp_per_nucl


    File : :file:`ao_basis/aos_transp.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_expo_ordered_transp_per_nucl	(ao_prim_num_max,N_AOs_max,nucl_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_expo_ordered_transp`
       * :c:data:`ao_prim_num`
       * :c:data:`ao_prim_num_max`
       * :c:data:`nucl_aos_transposed`
       * :c:data:`nucl_n_aos`
       * :c:data:`nucl_num`


 
.. c:var:: ao_l


    File : :file:`ao_basis/aos.irp.f`

    .. code:: fortran

        integer, allocatable	:: ao_l	(ao_num)
        integer	:: ao_l_max	
        character*(128), allocatable	:: ao_l_char	(ao_num)


    :math:`l` value of the |AO|: :math`a+b+c` in :math:`x^a y^b z^c`

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`l_to_character`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_cart_to_sphe_coef`
       * :c:data:`ao_coef_normalization_libint_factor`
       * :c:data:`ao_l_char_space`
       * :c:data:`nucl_list_shell_aos`

 
.. c:var:: ao_l_char


    File : :file:`ao_basis/aos.irp.f`

    .. code:: fortran

        integer, allocatable	:: ao_l	(ao_num)
        integer	:: ao_l_max	
        character*(128), allocatable	:: ao_l_char	(ao_num)


    :math:`l` value of the |AO|: :math`a+b+c` in :math:`x^a y^b z^c`

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`l_to_character`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_cart_to_sphe_coef`
       * :c:data:`ao_coef_normalization_libint_factor`
       * :c:data:`ao_l_char_space`
       * :c:data:`nucl_list_shell_aos`

 
.. c:var:: ao_l_char_space


    File : :file:`ao_basis/aos.irp.f`

    .. code:: fortran

        character*(4), allocatable	:: ao_l_char_space	(ao_num)


    Converts an l value to a string

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_l`
       * :c:data:`ao_num`
       * :c:data:`ao_power`


 
.. c:var:: ao_l_max


    File : :file:`ao_basis/aos.irp.f`

    .. code:: fortran

        integer, allocatable	:: ao_l	(ao_num)
        integer	:: ao_l_max	
        character*(128), allocatable	:: ao_l_char	(ao_num)


    :math:`l` value of the |AO|: :math`a+b+c` in :math:`x^a y^b z^c`

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`l_to_character`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_cart_to_sphe_coef`
       * :c:data:`ao_coef_normalization_libint_factor`
       * :c:data:`ao_l_char_space`
       * :c:data:`nucl_list_shell_aos`

 
.. c:var:: ao_power_ordered_transp_per_nucl


    File : :file:`ao_basis/aos_transp.irp.f`

    .. code:: fortran

        integer, allocatable	:: ao_power_ordered_transp_per_nucl	(3,N_AOs_max,nucl_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_power`
       * :c:data:`nucl_aos_transposed`
       * :c:data:`nucl_n_aos`
       * :c:data:`nucl_num`


 
.. c:var:: ao_prim_num_max


    File : :file:`ao_basis/aos.irp.f`

    .. code:: fortran

        integer	:: ao_prim_num_max	


    Max number of primitives.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_prim_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef`
       * :c:data:`ao_coef_normalized`
       * :c:data:`ao_coef_normalized_ordered`
       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_coef_normalized_ordered_transp_per_nucl`
       * :c:data:`ao_expo`
       * :c:data:`ao_expo_ordered_transp`
       * :c:data:`ao_expo_ordered_transp_per_nucl`

 
.. c:var:: cart_to_sphe_0


    File : :file:`ao_basis/spherical_to_cartesian.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cart_to_sphe_0	(1,1)


    Spherical -> Cartesian Transformation matrix for l=0


 
.. c:var:: cart_to_sphe_1


    File : :file:`ao_basis/spherical_to_cartesian.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cart_to_sphe_1	(3,3)


    Spherical -> Cartesian Transformation matrix for l=1

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_cart_to_sphe_coef`

 
.. c:var:: cart_to_sphe_2


    File : :file:`ao_basis/spherical_to_cartesian.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cart_to_sphe_2	(6,5)


    Spherical -> Cartesian Transformation matrix for l=2

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_cart_to_sphe_coef`

 
.. c:var:: cart_to_sphe_3


    File : :file:`ao_basis/spherical_to_cartesian.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cart_to_sphe_3	(10,7)


    Spherical -> Cartesian Transformation matrix for l=3

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_cart_to_sphe_coef`

 
.. c:var:: cart_to_sphe_4


    File : :file:`ao_basis/spherical_to_cartesian.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cart_to_sphe_4	(15,9)


    Spherical -> Cartesian Transformation matrix for l=4

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_cart_to_sphe_coef`

 
.. c:var:: cart_to_sphe_5


    File : :file:`ao_basis/spherical_to_cartesian.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cart_to_sphe_5	(21,11)


    Spherical -> Cartesian Transformation matrix for l=5

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_cart_to_sphe_coef`

 
.. c:var:: cart_to_sphe_6


    File : :file:`ao_basis/spherical_to_cartesian.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cart_to_sphe_6	(28,13)


    Spherical -> Cartesian Transformation matrix for l=6

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_cart_to_sphe_coef`

 
.. c:var:: cart_to_sphe_7


    File : :file:`ao_basis/spherical_to_cartesian.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cart_to_sphe_7	(36,15)


    Spherical -> Cartesian Transformation matrix for l=7

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_cart_to_sphe_coef`

 
.. c:var:: cart_to_sphe_8


    File : :file:`ao_basis/spherical_to_cartesian.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cart_to_sphe_8	(45,17)


    Spherical -> Cartesian Transformation matrix for l=8

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_cart_to_sphe_coef`

 
.. c:var:: cart_to_sphe_9


    File : :file:`ao_basis/spherical_to_cartesian.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cart_to_sphe_9	(55,19)


    Spherical -> Cartesian Transformation matrix for l=9

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_cart_to_sphe_coef`

 
.. c:var:: l_to_character


    File : :file:`ao_basis/aos.irp.f`

    .. code:: fortran

        character*(128), allocatable	:: l_to_character	(0:7)


    Character corresponding to the "l" value of an |AO|

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_l`

 
.. c:var:: n_aos_max


    File : :file:`ao_basis/aos.irp.f`

    .. code:: fortran

        integer, allocatable	:: nucl_n_aos	(nucl_num)
        integer	:: n_aos_max	


    Number of |AOs| per atom

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`nucl_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp_per_nucl`
       * :c:data:`ao_expo_ordered_transp_per_nucl`
       * :c:data:`ao_power_ordered_transp_per_nucl`
       * :c:data:`nucl_aos`
       * :c:data:`nucl_aos_transposed`
       * :c:data:`nucl_list_shell_aos`

 
.. c:var:: n_pt_max_i_x


    File : :file:`ao_basis/dimensions_integrals.irp.f`

    .. code:: fortran

        integer	:: n_pt_max_integrals	
        integer	:: n_pt_max_i_x	


    Number of points used in the numerical integrations.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_power`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_n_e`
       * :c:data:`ao_integrals_n_e_per_atom`
       * :c:data:`ao_two_e_integral_schwartz`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`gauleg_t2`

 
.. c:var:: n_pt_max_integrals


    File : :file:`ao_basis/dimensions_integrals.irp.f`

    .. code:: fortran

        integer	:: n_pt_max_integrals	
        integer	:: n_pt_max_i_x	


    Number of points used in the numerical integrations.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_power`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_n_e`
       * :c:data:`ao_integrals_n_e_per_atom`
       * :c:data:`ao_two_e_integral_schwartz`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`gauleg_t2`

 
.. c:var:: nucl_aos


    File : :file:`ao_basis/aos.irp.f`

    .. code:: fortran

        integer, allocatable	:: nucl_aos	(nucl_num,N_AOs_max)


    List of |AOs| centered on each atom

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`nucl_n_aos`
       * :c:data:`nucl_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_aos_transposed`
       * :c:data:`nucl_list_shell_aos`

 
.. c:var:: nucl_aos_transposed


    File : :file:`ao_basis/aos_transp.irp.f`

    .. code:: fortran

        integer, allocatable	:: nucl_aos_transposed	(N_AOs_max,nucl_num)


    List of AOs attached on each atom

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`nucl_aos`
       * :c:data:`nucl_n_aos`
       * :c:data:`nucl_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp_per_nucl`
       * :c:data:`ao_expo_ordered_transp_per_nucl`
       * :c:data:`ao_power_ordered_transp_per_nucl`

 
.. c:var:: nucl_list_shell_aos


    File : :file:`ao_basis/aos.irp.f`

    .. code:: fortran

        integer, allocatable	:: nucl_list_shell_aos	(nucl_num,N_AOs_max)
        integer, allocatable	:: nucl_num_shell_aos	(nucl_num)


    Index of the shell type |AOs| and of the corresponding |AOs|
    By convention, for p,d,f and g |AOs|, we take the index
    of the |AO| with the the corresponding power in the x axis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_l`
       * :c:data:`ao_power`
       * :c:data:`nucl_aos`
       * :c:data:`nucl_n_aos`
       * :c:data:`nucl_num`


 
.. c:var:: nucl_n_aos


    File : :file:`ao_basis/aos.irp.f`

    .. code:: fortran

        integer, allocatable	:: nucl_n_aos	(nucl_num)
        integer	:: n_aos_max	


    Number of |AOs| per atom

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`nucl_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp_per_nucl`
       * :c:data:`ao_expo_ordered_transp_per_nucl`
       * :c:data:`ao_power_ordered_transp_per_nucl`
       * :c:data:`nucl_aos`
       * :c:data:`nucl_aos_transposed`
       * :c:data:`nucl_list_shell_aos`

 
.. c:var:: nucl_num_shell_aos


    File : :file:`ao_basis/aos.irp.f`

    .. code:: fortran

        integer, allocatable	:: nucl_list_shell_aos	(nucl_num,N_AOs_max)
        integer, allocatable	:: nucl_num_shell_aos	(nucl_num)


    Index of the shell type |AOs| and of the corresponding |AOs|
    By convention, for p,d,f and g |AOs|, we take the index
    of the |AO| with the the corresponding power in the x axis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_l`
       * :c:data:`ao_power`
       * :c:data:`nucl_aos`
       * :c:data:`nucl_n_aos`
       * :c:data:`nucl_num`


 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: ao_power_index:


    File : :file:`ao_basis/aos.irp.f`

    .. code:: fortran

        integer function ao_power_index(nx,ny,nz)


    Unique index given to a triplet of powers:
    
    :math:`\frac{1}{2} (l-n_x) (l-n_x+1) + n_z + 1`

 
.. c:function:: ao_value:


    File : :file:`ao_basis/aos_in_r.irp.f`

    .. code:: fortran

        double precision function ao_value(i,r)


    Returns the value of the i-th ao at point $\textbf{r}$

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_expo_ordered_transp`
       * :c:data:`ao_nucl`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`nucl_coord`

 
.. c:function:: give_all_aos_and_grad_and_lapl_at_r:


    File : :file:`ao_basis/aos_in_r.irp.f`

    .. code:: fortran

        subroutine give_all_aos_and_grad_and_lapl_at_r(r,aos_array,aos_grad_array,aos_lapl_array)


    input  : r(1) ==> r(1) = x, r(2) = y, r(3) = z
    
    output :
    
    * aos_array(i) = ao(i) evaluated at $\textbf{r}$
    * aos_grad_array(1,i) = $\nabla_x$ of the ao(i) evaluated at $\textbf{r}$

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp_per_nucl`
       * :c:data:`ao_expo_ordered_transp_per_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power_ordered_transp_per_nucl`
       * :c:data:`ao_prim_num`
       * :c:data:`nucl_aos_transposed`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_n_aos`
       * :c:data:`nucl_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`give_all_mos_and_grad_and_lapl_at_r`

 
.. c:function:: give_all_aos_and_grad_at_r:


    File : :file:`ao_basis/aos_in_r.irp.f`

    .. code:: fortran

        subroutine give_all_aos_and_grad_at_r(r,aos_array,aos_grad_array)


    input : r(1) ==> r(1) = x, r(2) = y, r(3) = z
    
    output :
    
    * aos_array(i) = ao(i) evaluated at ro
    * aos_grad_array(1,i) = gradient X of the ao(i) evaluated at $\textbf{r}$
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp_per_nucl`
       * :c:data:`ao_expo_ordered_transp_per_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power_ordered_transp_per_nucl`
       * :c:data:`ao_prim_num`
       * :c:data:`nucl_aos_transposed`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_n_aos`
       * :c:data:`nucl_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`give_all_mos_and_grad_at_r`

 
.. c:function:: give_all_aos_at_r:


    File : :file:`ao_basis/aos_in_r.irp.f`

    .. code:: fortran

        subroutine give_all_aos_at_r(r,aos_array)


    input  : r == r(1) = x and so on
    
    output : aos_array(i) = aos(i) evaluated in $\textbf{r}$

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp_per_nucl`
       * :c:data:`ao_expo_ordered_transp_per_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power_ordered_transp_per_nucl`
       * :c:data:`ao_prim_num`
       * :c:data:`nucl_aos_transposed`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_n_aos`
       * :c:data:`nucl_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`give_all_mos_at_r`

 
.. c:function:: primitive_value:


    File : :file:`ao_basis/aos_in_r.irp.f`

    .. code:: fortran

        double precision function primitive_value(i,j,r)


    Returns the value of the j-th primitive of the i-th |AO| at point $\textbf{r}
    **without the coefficient**

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_expo_ordered_transp`
       * :c:data:`ao_nucl`
       * :c:data:`ao_power`
       * :c:data:`nucl_coord`

