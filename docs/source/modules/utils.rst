.. _module_utils: 
 
.. program:: utils 
 
.. default-role:: option 
 
=====
utils
=====

Contains general purpose utilities (sorting, maps, etc).

 
 
 
Providers 
--------- 
 
.. c:var:: binom


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        double precision, allocatable	:: binom	(0:40,0:40)
        double precision, allocatable	:: binom_transp	(0:40,0:40)


    Binomial coefficients

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`binom_int`

 
.. c:var:: binom_int


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        integer*8, allocatable	:: binom_int	(0:40,0:40)
        integer*8, allocatable	:: binom_int_transp	(0:40,0:40)


    Binomial coefficients, as integers*8

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`binom`


 
.. c:var:: binom_int_transp


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        integer*8, allocatable	:: binom_int	(0:40,0:40)
        integer*8, allocatable	:: binom_int_transp	(0:40,0:40)


    Binomial coefficients, as integers*8

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`binom`


 
.. c:var:: binom_transp


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        double precision, allocatable	:: binom	(0:40,0:40)
        double precision, allocatable	:: binom_transp	(0:40,0:40)


    Binomial coefficients

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`binom_int`

 
.. c:var:: degree_max_integration_lebedev


    File : :file:`utils/angular_integration.irp.f`

    .. code:: fortran

        integer	:: degree_max_integration_lebedev	


    integrate correctly a polynom of order "degree_max_integration_lebedev"
    needed for the angular integration according to LEBEDEV formulae

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`n_points_integration_angular_lebedev`
       * :c:data:`theta_angular_integration_lebedev`

 
.. c:function:: dtranspose:


    File : :file:`utils/transpose.irp.f`

    .. code:: fortran

        recursive subroutine dtranspose(A,LDA,B,LDB,d1,d2)


    Transpose input matrix A into output matrix B

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`dtranspose`
       * :c:func:`h_s2_u_0_nstates_openmp`
       * :c:func:`h_s2_u_0_nstates_zmq`
       * :c:func:`h_s2_u_0_two_e_nstates_openmp`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dtranspose`

 
.. c:var:: fact_inv


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fact_inv	(128)


    1/n!


 
.. c:function:: i2radix_sort:


    File : :file:`utils/sort.irp.f_template_644`

    .. code:: fortran

        recursive subroutine i2radix_sort(x,iorder,isize,iradix)


    Sort integer array x(isize) using the radix sort algorithm.
    iorder in input should be (1,2,3,...,isize), and in output
    contains the new order of the elements.
    iradix should be -1 in input.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_mo_two_e_integrals_erf_i1j1`
       * :c:func:`get_mo_two_e_integrals_erf_ij`
       * :c:func:`get_mo_two_e_integrals_i1j1`
       * :c:func:`get_mo_two_e_integrals_ij`
       * :c:func:`i2radix_sort`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i2radix_sort`
       * :c:func:`insertion_i2sort`

 
.. c:function:: i8radix_sort:


    File : :file:`utils/sort.irp.f_template_644`

    .. code:: fortran

        recursive subroutine i8radix_sort(x,iorder,isize,iradix)


    Sort integer array x(isize) using the radix sort algorithm.
    iorder in input should be (1,2,3,...,isize), and in output
    contains the new order of the elements.
    iradix should be -1 in input.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_mo_two_e_integrals_erf_i1j1`
       * :c:func:`get_mo_two_e_integrals_erf_ij`
       * :c:func:`get_mo_two_e_integrals_i1j1`
       * :c:func:`get_mo_two_e_integrals_ij`
       * :c:func:`i8radix_sort`
       * :c:data:`psi_bilinear_matrix_transp_values`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i8radix_sort`
       * :c:func:`insertion_i8sort`

 
.. c:function:: i8radix_sort_big:


    File : :file:`utils/sort.irp.f_template_644`

    .. code:: fortran

        recursive subroutine i8radix_sort_big(x,iorder,isize,iradix)


    Sort integer array x(isize) using the radix sort algorithm.
    iorder in input should be (1,2,3,...,isize), and in output
    contains the new order of the elements.
    iradix should be -1 in input.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`i8radix_sort_big`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i8radix_sort_big`
       * :c:func:`insertion_i8sort_big`

 
.. c:var:: inv_int


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        double precision, allocatable	:: inv_int	(128)


    1/i


 
.. c:function:: iradix_sort:


    File : :file:`utils/sort.irp.f_template_644`

    .. code:: fortran

        recursive subroutine iradix_sort(x,iorder,isize,iradix)


    Sort integer array x(isize) using the radix sort algorithm.
    iorder in input should be (1,2,3,...,isize), and in output
    contains the new order of the elements.
    iradix should be -1 in input.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_mo_two_e_integrals_erf_i1j1`
       * :c:func:`get_mo_two_e_integrals_erf_ij`
       * :c:func:`get_mo_two_e_integrals_i1j1`
       * :c:func:`get_mo_two_e_integrals_ij`
       * :c:func:`iradix_sort`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`insertion_isort`
       * :c:func:`iradix_sort`

 
.. c:function:: iradix_sort_big:


    File : :file:`utils/sort.irp.f_template_644`

    .. code:: fortran

        recursive subroutine iradix_sort_big(x,iorder,isize,iradix)


    Sort integer array x(isize) using the radix sort algorithm.
    iorder in input should be (1,2,3,...,isize), and in output
    contains the new order of the elements.
    iradix should be -1 in input.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`iradix_sort_big`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`insertion_isort_big`
       * :c:func:`iradix_sort_big`

 
.. c:var:: n_points_integration_angular_lebedev


    File : :file:`utils/angular_integration.irp.f`

    .. code:: fortran

        integer	:: n_points_integration_angular_lebedev	


    Number of points needed for the angular integral

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`degree_max_integration_lebedev`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`theta_angular_integration_lebedev`

 
.. c:var:: nproc


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        integer	:: nproc	


    Number of current OpenMP threads

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_integrals_erf_in_map`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`h_apply_buffer_allocated`
       * :c:data:`n_det`
       * :c:data:`nthreads_davidson`
       * :c:data:`nthreads_pt2`

 
.. c:function:: overlap_gaussian_xyz:


    File : :file:`utils/one_e_integration.irp.f`

    .. code:: fortran

        subroutine overlap_gaussian_xyz(A_center,B_center,alpha,beta,power_A,&
      power_B,overlap_x,overlap_y,overlap_z,overlap,dim)


    .. math::
    
       S_x = \int (x-A_x)^{a_x} exp(-\alpha(x-A_x)^2)  (x-B_x)^{b_x} exp(-beta(x-B_x)^2) dx \\
       S = S_x S_y S_z
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalization_libint_factor`
       * :c:data:`ao_coef_normalized`
       * :c:data:`ao_deriv2_x`
       * :c:data:`ao_deriv_1_x`
       * :c:data:`ao_dipole_x`
       * :c:data:`ao_overlap`
       * :c:data:`ao_spread_x`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`gaussian_product_x`
       * :c:func:`give_explicit_poly_and_gaussian`

 
.. c:var:: phi_angular_integration_lebedev


    File : :file:`utils/angular_integration.irp.f`

    .. code:: fortran

        double precision, allocatable	:: theta_angular_integration_lebedev	(n_points_integration_angular_lebedev)
        double precision, allocatable	:: phi_angular_integration_lebedev	(n_points_integration_angular_lebedev)
        double precision, allocatable	:: weights_angular_integration_lebedev	(n_points_integration_angular_lebedev)


    Theta phi values together with the weights values for the angular integration :
    integral [dphi,dtheta] f(x,y,z) = 4 * pi * sum (1<i<n_points_integration_angular_lebedev) f(xi,yi,zi)
    Note that theta and phi are in DEGREES !!

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`degree_max_integration_lebedev`
       * :c:data:`n_points_integration_angular_lebedev`


 
.. c:var:: qp_max_mem


    File : :file:`utils/memory.irp.f`

    .. code:: fortran

        integer	:: qp_max_mem	


    Maximum memory in Gb

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pt2_j`
       * :c:data:`pt2_w`

 
.. c:function:: rec__quicksort:


    File : :file:`utils/sort.irp.f_template_261`

    .. code:: fortran

        recursive subroutine rec__quicksort(x, iorder, isize, first, last, level)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`quick_sort`
       * :c:func:`rec__quicksort`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`rec__quicksort`

 
.. c:function:: rec_d_quicksort:


    File : :file:`utils/sort.irp.f_template_261`

    .. code:: fortran

        recursive subroutine rec_d_quicksort(x, iorder, isize, first, last, level)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`quick_dsort`
       * :c:func:`rec_d_quicksort`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`rec_d_quicksort`

 
.. c:function:: rec_i2_quicksort:


    File : :file:`utils/sort.irp.f_template_261`

    .. code:: fortran

        recursive subroutine rec_i2_quicksort(x, iorder, isize, first, last, level)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`quick_i2sort`
       * :c:func:`rec_i2_quicksort`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`rec_i2_quicksort`

 
.. c:function:: rec_i8_quicksort:


    File : :file:`utils/sort.irp.f_template_261`

    .. code:: fortran

        recursive subroutine rec_i8_quicksort(x, iorder, isize, first, last, level)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`quick_i8sort`
       * :c:func:`rec_i8_quicksort`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`rec_i8_quicksort`

 
.. c:function:: rec_i_quicksort:


    File : :file:`utils/sort.irp.f_template_261`

    .. code:: fortran

        recursive subroutine rec_i_quicksort(x, iorder, isize, first, last, level)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`quick_isort`
       * :c:func:`rec_i_quicksort`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`rec_i_quicksort`

 
.. c:var:: theta_angular_integration_lebedev


    File : :file:`utils/angular_integration.irp.f`

    .. code:: fortran

        double precision, allocatable	:: theta_angular_integration_lebedev	(n_points_integration_angular_lebedev)
        double precision, allocatable	:: phi_angular_integration_lebedev	(n_points_integration_angular_lebedev)
        double precision, allocatable	:: weights_angular_integration_lebedev	(n_points_integration_angular_lebedev)


    Theta phi values together with the weights values for the angular integration :
    integral [dphi,dtheta] f(x,y,z) = 4 * pi * sum (1<i<n_points_integration_angular_lebedev) f(xi,yi,zi)
    Note that theta and phi are in DEGREES !!

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`degree_max_integration_lebedev`
       * :c:data:`n_points_integration_angular_lebedev`


 
.. c:function:: transpose:


    File : :file:`utils/transpose.irp.f`

    .. code:: fortran

        recursive subroutine transpose(A,LDA,B,LDB,d1,d2)


    Transpose input matrix A into output matrix B

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`transpose`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`transpose`

 
.. c:var:: weights_angular_integration_lebedev


    File : :file:`utils/angular_integration.irp.f`

    .. code:: fortran

        double precision, allocatable	:: theta_angular_integration_lebedev	(n_points_integration_angular_lebedev)
        double precision, allocatable	:: phi_angular_integration_lebedev	(n_points_integration_angular_lebedev)
        double precision, allocatable	:: weights_angular_integration_lebedev	(n_points_integration_angular_lebedev)


    Theta phi values together with the weights values for the angular integration :
    integral [dphi,dtheta] f(x,y,z) = 4 * pi * sum (1<i<n_points_integration_angular_lebedev) f(xi,yi,zi)
    Note that theta and phi are in DEGREES !!

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`degree_max_integration_lebedev`
       * :c:data:`n_points_integration_angular_lebedev`


 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: a_coef:


    File : :file:`utils/need.irp.f`

    .. code:: fortran

        double precision function a_coef(n)



 
.. c:function:: add_poly:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        subroutine add_poly(b,nb,c,nc,d,nd)


    Add two polynomials
    D(t) =! D(t) +( B(t)+C(t))

 
.. c:function:: add_poly_multiply:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        subroutine add_poly_multiply(b,nb,cst,d,nd)


    Add a polynomial multiplied by a constant
    D(t) =! D(t) +( cst * B(t))

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`general_primitive_integral`
       * :c:func:`general_primitive_integral_erf`

 
.. c:function:: apply_rotation:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine apply_rotation(A,LDA,R,LDR,B,LDB,m,n)


    Apply the rotation found by find_rotation

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

 
.. c:function:: approx_dble:


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        double precision function approx_dble(a,n)



 
.. c:function:: b_coef:


    File : :file:`utils/need.irp.f`

    .. code:: fortran

        double precision function b_coef(n,u)



 
.. c:function:: binom_func:


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        double precision function binom_func(i,j)


    .. math                       ::
    
      \frac{i!}{j!(i-j)!}
    

 
.. c:function:: check_mem:


    File : :file:`utils/memory.irp.f`

    .. code:: fortran

        subroutine check_mem(rss_in,routine)


    Checks if n gigabytes can be allocated. If not, exit the run.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`qp_max_mem`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`create_selection_buffer`
       * :c:func:`davidson_diag_hjj_sjj`
       * :c:func:`make_selection_buffer_s2`
       * :c:func:`merge_selection_buffers`
       * :c:func:`pt2_collector`
       * :c:data:`pt2_j`
       * :c:data:`pt2_w`
       * :c:func:`remove_duplicates_in_selection_buffer`
       * :c:func:`run_cipsi`
       * :c:func:`run_pt2_slave`
       * :c:func:`run_stochastic_cipsi`
       * :c:func:`select_singles_and_doubles`
       * :c:func:`selection_collector`
       * :c:func:`sort_selection_buffer`
       * :c:func:`testteethbuilding`
       * :c:func:`zmq_pt2`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`resident_memory`

 
.. c:function:: dble_fact:


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        double precision function dble_fact(n)



 
.. c:function:: dble_fact_even:


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        double precision function dble_fact_even(n) result(fact2)


    n!!

 
.. c:function:: dble_fact_odd:


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        double precision function dble_fact_odd(n) result(fact2)


    n!!

 
.. c:function:: dble_logfact:


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        double precision function dble_logfact(n) result(logfact2)


    n!!

 
.. c:function:: ddfact2:


    File : :file:`utils/need.irp.f`

    .. code:: fortran

        double precision function ddfact2(n)



 
.. c:function:: dset_order:


    File : :file:`utils/sort.irp.f_template_347`

    .. code:: fortran

        subroutine dset_order(x,iorder,isize)


    array A has already been sorted, and iorder has contains the new order of
    elements of A. This subroutine changes the order of x to match the new order of A.

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered`
       * :c:func:`h_s2_u_0_nstates_openmp`
       * :c:func:`h_s2_u_0_nstates_zmq`
       * :c:func:`h_s2_u_0_two_e_nstates_openmp`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`

 
.. c:function:: dset_order_big:


    File : :file:`utils/sort.irp.f_template_412`

    .. code:: fortran

        subroutine dset_order_big(x,iorder,isize)


    array A has already been sorted, and iorder has contains the new order of
    elements of A. This subroutine changes the order of x to match the new order of A.
    This is a version for very large arrays where the indices need
    to be in integer*8 format

 
.. c:function:: dsort:


    File : :file:`utils/sort.irp.f_template_293`

    .. code:: fortran

        subroutine dsort(x,iorder,isize)


    Sort array x(isize).
    iorder in input should be (1,2,3,...,isize), and in output
    contains the new order of the elements.

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered`
       * :c:func:`make_selection_buffer_s2`
       * :c:data:`psi_det_sorted`
       * :c:func:`reorder_core_orb`
       * :c:func:`sort_selection_buffer`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`insertion_dsort`
       * :c:func:`quick_dsort`

 
.. c:function:: erf0:


    File : :file:`utils/need.irp.f`

    .. code:: fortran

        double precision function erf0(x)



 
.. c:function:: extrapolate_data:


    File : :file:`utils/extrapolation.irp.f`

    .. code:: fortran

        subroutine extrapolate_data(N_data, data, pt2, output)


    Extrapolate the data to the FCI limit

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`extrapolated_energy`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_pseudo_inverse`

 
.. c:function:: f_integral:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        double precision function F_integral(n,p)


    function that calculates the following integral
    \int_{\-infty}^{+\infty} x^n \exp(-p x^2) dx

 
.. c:function:: fact:


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        double precision function fact(n)


    n!

 
.. c:function:: find_rotation:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine find_rotation(A,LDA,B,m,C,n)


    Find A.C = B

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`get_pseudo_inverse`

 
.. c:function:: gammln:


    File : :file:`utils/need.irp.f`

    .. code:: fortran

        double precision function gammln(xx)



 
.. c:function:: gammp:


    File : :file:`utils/need.irp.f`

    .. code:: fortran

        double precision function gammp(a,x)



    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`gcf`
       * :c:func:`gser`

 
.. c:function:: gaussian_product:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        subroutine gaussian_product(a,xa,b,xb,k,p,xp)


    Gaussian product in 1D.
    e^{-a (x-x_A)^2} e^{-b (x-x_B)^2} = K_{ab}^x e^{-p (x-x_P)^2}

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`give_explicit_poly_and_gaussian`
       * :c:func:`give_explicit_poly_and_gaussian_double`

 
.. c:function:: gaussian_product_x:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        subroutine gaussian_product_x(a,xa,b,xb,k,p,xp)


    Gaussian product in 1D.
    e^{-a (x-x_A)^2} e^{-b (x-x_B)^2} = K_{ab}^x e^{-p (x-x_P)^2}

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`overlap_gaussian_xyz`

 
.. c:function:: gcf:


    File : :file:`utils/need.irp.f`

    .. code:: fortran

        subroutine gcf(gammcf,a,x,gln)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`gammp`

 
.. c:function:: get_inverse:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine get_inverse(A,LDA,m,C,LDC)


    Returns the inverse of the square matrix A

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_ortho_canonical_coef_inv`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgetrf`
       * :c:func:`dgetri`

 
.. c:function:: get_pseudo_inverse:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine get_pseudo_inverse(A,LDA,m,n,C,LDC)


    Find C = A^-1

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_cart_to_sphe_inv`
       * :c:func:`extrapolate_data`
       * :c:func:`find_rotation`
       * :c:data:`s_inv`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgesvd`

 
.. c:function:: give_explicit_poly_and_gaussian:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        subroutine give_explicit_poly_and_gaussian(P_new,P_center,p,fact_k,iorder,alpha,beta,a,b,A_center,B_center,dim)


    Transforms the product of
             (x-x_A)^a(1) (x-x_B)^b(1) (x-x_A)^a(2) (y-y_B)^b(2) (z-z_A)^a(3) (z-z_B)^b(3) exp(-(r-A)^2 alpha) exp(-(r-B)^2 beta)
    into
           fact_k * [ sum (l_x = 0,i_order(1)) P_new(l_x,1) * (x-P_center(1))^l_x ] exp (- p (x-P_center(1))^2 )
                  * [ sum (l_y = 0,i_order(2)) P_new(l_y,2) * (y-P_center(2))^l_y ] exp (- p (y-P_center(2))^2 )
                  * [ sum (l_z = 0,i_order(3)) P_new(l_z,3) * (z-P_center(3))^l_z ] exp (- p (z-P_center(3))^2 )

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ao_two_e_integral`
       * :c:func:`ao_two_e_integral_erf`
       * :c:func:`ao_two_e_integral_schwartz_accel`
       * :c:func:`ao_two_e_integral_schwartz_accel_erf`
       * :c:func:`give_explicit_poly_and_gaussian_double`
       * :c:func:`overlap_gaussian_xyz`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`gaussian_product`
       * :c:func:`multiply_poly`
       * :c:func:`recentered_poly2`

 
.. c:function:: give_explicit_poly_and_gaussian_double:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        subroutine give_explicit_poly_and_gaussian_double(P_new,P_center,p,fact_k,iorder,alpha,beta,gama,a,b,A_center,B_center,Nucl_center,dim)


    Transforms the product of
             (x-x_A)^a(1) (x-x_B)^b(1) (x-x_A)^a(2) (y-y_B)^b(2) (z-z_A)^a(3) (z-z_B)^b(3)
             exp(-(r-A)^2 alpha) exp(-(r-B)^2 beta) exp(-(r-Nucl_center)^2 gama
    
    into
           fact_k * [ sum (l_x = 0,i_order(1)) P_new(l_x,1) * (x-P_center(1))^l_x ] exp (- p (x-P_center(1))^2 )
                  * [ sum (l_y = 0,i_order(2)) P_new(l_y,2) * (y-P_center(2))^l_y ] exp (- p (y-P_center(2))^2 )
                  * [ sum (l_z = 0,i_order(3)) P_new(l_z,3) * (z-P_center(3))^l_z ] exp (- p (z-P_center(3))^2 )

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`gaussian_product`
       * :c:func:`give_explicit_poly_and_gaussian`

 
.. c:function:: give_explicit_poly_and_gaussian_x:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        subroutine give_explicit_poly_and_gaussian_x(P_new,P_center,p,fact_k,iorder,alpha,beta,a,b,A_center,B_center,dim)


    Transform the product of
             (x-x_A)^a(1) (x-x_B)^b(1) (x-x_A)^a(2) (y-y_B)^b(2) (z-z_A)^a(3) (z-z_B)^b(3) exp(-(r-A)^2 alpha) exp(-(r-B)^2 beta)
    into
           fact_k  (x-x_P)^iorder(1)  (y-y_P)^iorder(2)  (z-z_P)^iorder(3) exp(-p(r-P)^2)

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`overlap_gaussian_x`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`multiply_poly`
       * :c:func:`recentered_poly2`

 
.. c:function:: gser:


    File : :file:`utils/need.irp.f`

    .. code:: fortran

        subroutine gser(gamser,a,x,gln)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`gammp`

 
.. c:function:: heap_dsort:


    File : :file:`utils/sort.irp.f_template_261`

    .. code:: fortran

        subroutine heap_dsort(x,iorder,isize)


    Sort array x(isize) using the heap sort algorithm.
    iorder in input should be (1,2,3,...,isize), and in output
    contains the new order of the elements.

 
.. c:function:: heap_dsort_big:


    File : :file:`utils/sort.irp.f_template_261`

    .. code:: fortran

        subroutine heap_dsort_big(x,iorder,isize)


    Sort array x(isize) using the heap sort algorithm.
    iorder in input should be (1,2,3,...,isize), and in output
    contains the new order of the elements.
    This is a version for very large arrays where the indices need
    to be in integer*8 format

 
.. c:function:: heap_i2sort:


    File : :file:`utils/sort.irp.f_template_261`

    .. code:: fortran

        subroutine heap_i2sort(x,iorder,isize)


    Sort array x(isize) using the heap sort algorithm.
    iorder in input should be (1,2,3,...,isize), and in output
    contains the new order of the elements.

 
.. c:function:: heap_i2sort_big:


    File : :file:`utils/sort.irp.f_template_261`

    .. code:: fortran

        subroutine heap_i2sort_big(x,iorder,isize)


    Sort array x(isize) using the heap sort algorithm.
    iorder in input should be (1,2,3,...,isize), and in output
    contains the new order of the elements.
    This is a version for very large arrays where the indices need
    to be in integer*8 format

 
.. c:function:: heap_i8sort:


    File : :file:`utils/sort.irp.f_template_261`

    .. code:: fortran

        subroutine heap_i8sort(x,iorder,isize)


    Sort array x(isize) using the heap sort algorithm.
    iorder in input should be (1,2,3,...,isize), and in output
    contains the new order of the elements.

 
.. c:function:: heap_i8sort_big:


    File : :file:`utils/sort.irp.f_template_261`

    .. code:: fortran

        subroutine heap_i8sort_big(x,iorder,isize)


    Sort array x(isize) using the heap sort algorithm.
    iorder in input should be (1,2,3,...,isize), and in output
    contains the new order of the elements.
    This is a version for very large arrays where the indices need
    to be in integer*8 format

 
.. c:function:: heap_isort:


    File : :file:`utils/sort.irp.f_template_261`

    .. code:: fortran

        subroutine heap_isort(x,iorder,isize)


    Sort array x(isize) using the heap sort algorithm.
    iorder in input should be (1,2,3,...,isize), and in output
    contains the new order of the elements.

 
.. c:function:: heap_isort_big:


    File : :file:`utils/sort.irp.f_template_261`

    .. code:: fortran

        subroutine heap_isort_big(x,iorder,isize)


    Sort array x(isize) using the heap sort algorithm.
    iorder in input should be (1,2,3,...,isize), and in output
    contains the new order of the elements.
    This is a version for very large arrays where the indices need
    to be in integer*8 format

 
.. c:function:: heap_sort:


    File : :file:`utils/sort.irp.f_template_261`

    .. code:: fortran

        subroutine heap_sort(x,iorder,isize)


    Sort array x(isize) using the heap sort algorithm.
    iorder in input should be (1,2,3,...,isize), and in output
    contains the new order of the elements.

 
.. c:function:: heap_sort_big:


    File : :file:`utils/sort.irp.f_template_261`

    .. code:: fortran

        subroutine heap_sort_big(x,iorder,isize)


    Sort array x(isize) using the heap sort algorithm.
    iorder in input should be (1,2,3,...,isize), and in output
    contains the new order of the elements.
    This is a version for very large arrays where the indices need
    to be in integer*8 format

 
.. c:function:: hermite:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        double precision function hermite(n,x)


    Hermite polynomial

 
.. c:function:: i2set_order:


    File : :file:`utils/sort.irp.f_template_347`

    .. code:: fortran

        subroutine i2set_order(x,iorder,isize)


    array A has already been sorted, and iorder has contains the new order of
    elements of A. This subroutine changes the order of x to match the new order of A.

 
.. c:function:: i2set_order_big:


    File : :file:`utils/sort.irp.f_template_412`

    .. code:: fortran

        subroutine i2set_order_big(x,iorder,isize)


    array A has already been sorted, and iorder has contains the new order of
    elements of A. This subroutine changes the order of x to match the new order of A.
    This is a version for very large arrays where the indices need
    to be in integer*8 format

 
.. c:function:: i2sort:


    File : :file:`utils/sort.irp.f_template_315`

    .. code:: fortran

        subroutine i2sort(x,iorder,isize)


    Sort array x(isize).
    iorder in input should be (1,2,3,...,isize), and in output
    contains the new order of the elements.

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`quick_i2sort`

 
.. c:function:: i8set_order:


    File : :file:`utils/sort.irp.f_template_347`

    .. code:: fortran

        subroutine i8set_order(x,iorder,isize)


    array A has already been sorted, and iorder has contains the new order of
    elements of A. This subroutine changes the order of x to match the new order of A.

 
.. c:function:: i8set_order_big:


    File : :file:`utils/sort.irp.f_template_412`

    .. code:: fortran

        subroutine i8set_order_big(x,iorder,isize)


    array A has already been sorted, and iorder has contains the new order of
    elements of A. This subroutine changes the order of x to match the new order of A.
    This is a version for very large arrays where the indices need
    to be in integer*8 format

 
.. c:function:: i8sort:


    File : :file:`utils/sort.irp.f_template_315`

    .. code:: fortran

        subroutine i8sort(x,iorder,isize)


    Sort array x(isize).
    iorder in input should be (1,2,3,...,isize), and in output
    contains the new order of the elements.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`make_selection_buffer_s2`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`psi_occ_pattern`
       * :c:func:`remove_duplicates_in_selection_buffer`
       * :c:func:`sort_dets_by_det_search_key`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`quick_i8sort`

 
.. c:function:: insertion_dsort:


    File : :file:`utils/sort.irp.f_template_261`

    .. code:: fortran

        subroutine insertion_dsort (x,iorder,isize)


    Sort array x(isize) using the insertion sort algorithm.
    iorder in input should be (1,2,3,...,isize), and in output
    contains the new order of the elements.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`dsort`

 
.. c:function:: insertion_dsort_big:


    File : :file:`utils/sort.irp.f_template_412`

    .. code:: fortran

        subroutine insertion_dsort_big (x,iorder,isize)


    Sort array x(isize) using the insertion sort algorithm.
    iorder in input should be (1,2,3,...,isize), and in output
    contains the new order of the elements.
    This is a version for very large arrays where the indices need
    to be in integer*8 format

 
.. c:function:: insertion_i2sort:


    File : :file:`utils/sort.irp.f_template_261`

    .. code:: fortran

        subroutine insertion_i2sort (x,iorder,isize)


    Sort array x(isize) using the insertion sort algorithm.
    iorder in input should be (1,2,3,...,isize), and in output
    contains the new order of the elements.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`i2radix_sort`

 
.. c:function:: insertion_i2sort_big:


    File : :file:`utils/sort.irp.f_template_412`

    .. code:: fortran

        subroutine insertion_i2sort_big (x,iorder,isize)


    Sort array x(isize) using the insertion sort algorithm.
    iorder in input should be (1,2,3,...,isize), and in output
    contains the new order of the elements.
    This is a version for very large arrays where the indices need
    to be in integer*8 format

 
.. c:function:: insertion_i8sort:


    File : :file:`utils/sort.irp.f_template_261`

    .. code:: fortran

        subroutine insertion_i8sort (x,iorder,isize)


    Sort array x(isize) using the insertion sort algorithm.
    iorder in input should be (1,2,3,...,isize), and in output
    contains the new order of the elements.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`i8radix_sort`

 
.. c:function:: insertion_i8sort_big:


    File : :file:`utils/sort.irp.f_template_412`

    .. code:: fortran

        subroutine insertion_i8sort_big (x,iorder,isize)


    Sort array x(isize) using the insertion sort algorithm.
    iorder in input should be (1,2,3,...,isize), and in output
    contains the new order of the elements.
    This is a version for very large arrays where the indices need
    to be in integer*8 format

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`i8radix_sort_big`

 
.. c:function:: insertion_isort:


    File : :file:`utils/sort.irp.f_template_261`

    .. code:: fortran

        subroutine insertion_isort (x,iorder,isize)


    Sort array x(isize) using the insertion sort algorithm.
    iorder in input should be (1,2,3,...,isize), and in output
    contains the new order of the elements.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`iradix_sort`

 
.. c:function:: insertion_isort_big:


    File : :file:`utils/sort.irp.f_template_412`

    .. code:: fortran

        subroutine insertion_isort_big (x,iorder,isize)


    Sort array x(isize) using the insertion sort algorithm.
    iorder in input should be (1,2,3,...,isize), and in output
    contains the new order of the elements.
    This is a version for very large arrays where the indices need
    to be in integer*8 format

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`iradix_sort_big`

 
.. c:function:: insertion_sort:


    File : :file:`utils/sort.irp.f_template_261`

    .. code:: fortran

        subroutine insertion_sort (x,iorder,isize)


    Sort array x(isize) using the insertion sort algorithm.
    iorder in input should be (1,2,3,...,isize), and in output
    contains the new order of the elements.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`sort`

 
.. c:function:: insertion_sort_big:


    File : :file:`utils/sort.irp.f_template_412`

    .. code:: fortran

        subroutine insertion_sort_big (x,iorder,isize)


    Sort array x(isize) using the insertion sort algorithm.
    iorder in input should be (1,2,3,...,isize), and in output
    contains the new order of the elements.
    This is a version for very large arrays where the indices need
    to be in integer*8 format

 
.. c:function:: iset_order:


    File : :file:`utils/sort.irp.f_template_347`

    .. code:: fortran

        subroutine iset_order(x,iorder,isize)


    array A has already been sorted, and iorder has contains the new order of
    elements of A. This subroutine changes the order of x to match the new order of A.

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`

 
.. c:function:: iset_order_big:


    File : :file:`utils/sort.irp.f_template_412`

    .. code:: fortran

        subroutine iset_order_big(x,iorder,isize)


    array A has already been sorted, and iorder has contains the new order of
    elements of A. This subroutine changes the order of x to match the new order of A.
    This is a version for very large arrays where the indices need
    to be in integer*8 format

 
.. c:function:: isort:


    File : :file:`utils/sort.irp.f_template_315`

    .. code:: fortran

        subroutine isort(x,iorder,isize)


    Sort array x(isize).
    iorder in input should be (1,2,3,...,isize), and in output
    contains the new order of the elements.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`select_singles_and_doubles`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`quick_isort`

 
.. c:function:: lapack_diag:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine lapack_diag(eigvalues,eigvectors,H,nmax,n)


    Diagonalize matrix H
    
    H is untouched between input and ouptut
    
    eigevalues(i) = ith lowest eigenvalue of the H matrix
    
    eigvectors(i,j) = <i|psi_j> where i is the basis function and psi_j is the j th eigenvector
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ci_electronic_energy`
       * :c:func:`davidson_diag_hjj_sjj`
       * :c:func:`mo_as_eigvectors_of_mo_matrix`
       * :c:data:`psi_coef_cas_diagonalized`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dsyev`

 
.. c:function:: lapack_diagd:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine lapack_diagd(eigvalues,eigvectors,H,nmax,n)


    Diagonalize matrix H
    
    H is untouched between input and ouptut
    
    eigevalues(i) = ith lowest eigenvalue of the H matrix
    
    eigvectors(i,j) = <i|psi_j> where i is the basis function and psi_j is the j th eigenvector
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`inertia_tensor_eigenvectors`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dsyevd`

 
.. c:function:: logfact:


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        double precision function logfact(n)


    n!

 
.. c:function:: lowercase:


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        subroutine lowercase(txt,n)


    Transform to lower case

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`end_parallel_job`
       * :c:func:`new_parallel_job`

 
.. c:function:: map_load_from_disk:


    File : :file:`utils/map_functions.irp.f`

    .. code:: fortran

        subroutine map_load_from_disk(filename,map)



    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_integrals_erf_in_map`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`mo_two_e_integrals_erf_in_map`
       * :c:data:`mo_two_e_integrals_in_map`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`c_f_pointer`
       * :c:func:`mmap`

 
.. c:function:: map_save_to_disk:


    File : :file:`utils/map_functions.irp.f`

    .. code:: fortran

        subroutine map_save_to_disk(filename,map)



    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_integrals_erf_in_map`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`mo_two_e_integrals_erf_in_map`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:func:`save_erf_two_e_integrals_ao`
       * :c:func:`save_erf_two_e_integrals_mo`
       * :c:func:`save_erf_two_e_ints_ao_into_ints_ao`
       * :c:func:`save_erf_two_e_ints_mo_into_ints_mo`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`c_f_pointer`
       * :c:func:`map_sort`
       * :c:func:`mmap`
       * :c:func:`msync`

 
.. c:function:: memory_of_double:


    File : :file:`utils/memory.irp.f`

    .. code:: fortran

        double precision function memory_of_double(n)


    Computes the memory required for n double precision elements in gigabytes.

 
.. c:function:: memory_of_int:


    File : :file:`utils/memory.irp.f`

    .. code:: fortran

        double precision function memory_of_int(n)


    Computes the memory required for n double precision elements in gigabytes.

 
.. c:function:: multiply_poly:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        subroutine multiply_poly(b,nb,c,nc,d,nd)


    Multiply two polynomials
    D(t) =! D(t) +( B(t)*C(t))

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`general_primitive_integral`
       * :c:func:`general_primitive_integral_erf`
       * :c:func:`give_explicit_poly_and_gaussian`
       * :c:func:`give_explicit_poly_and_gaussian_x`
       * :c:func:`give_polynomial_mult_center_one_e`
       * :c:func:`give_polynomial_mult_center_one_e_erf`
       * :c:func:`give_polynomial_mult_center_one_e_erf_opt`
       * :c:func:`i_x1_pol_mult_a1`
       * :c:func:`i_x1_pol_mult_a2`
       * :c:func:`i_x1_pol_mult_one_e`
       * :c:func:`i_x1_pol_mult_recurs`
       * :c:func:`i_x2_pol_mult`
       * :c:func:`i_x2_pol_mult_one_e`

 
.. c:function:: normalize:


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        subroutine normalize(u,sze)


    Normalizes vector u

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`copy_h_apply_buffer_to_wf`
       * :c:func:`davidson_diag_hjj_sjj`
       * :c:func:`save_wavefunction_general`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dscal`

 
.. c:function:: ortho_canonical:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine ortho_canonical(overlap,LDA,N,C,LDC,m)


    Compute C_new=C_old.U.s^-1/2 canonical orthogonalization.
    
    overlap : overlap matrix
    
    LDA : leftmost dimension of overlap array
    
    N : Overlap matrix is NxN (array is (LDA,N) )
    
    C : Coefficients of the vectors to orthogonalize. On exit,
        orthogonal vectors
    
    LDC : leftmost dimension of C
    
    m : Coefficients matrix is MxN, ( array is (LDC,N) )
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_ortho_canonical_coef`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`svd`

 
.. c:function:: ortho_lowdin:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine ortho_lowdin(overlap,LDA,N,C,LDC,m)


    Compute C_new=C_old.S^-1/2 orthogonalization.
    
    overlap : overlap matrix
    
    LDA : leftmost dimension of overlap array
    
    N : Overlap matrix is NxN (array is (LDA,N) )
    
    C : Coefficients of the vectors to orthogonalize. On exit,
        orthogonal vectors
    
    LDC : leftmost dimension of C
    
    M : Coefficients matrix is MxN, ( array is (LDC,N) )
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_ortho_lowdin_coef`
       * :c:func:`orthonormalize_mos`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`svd`

 
.. c:function:: ortho_qr:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine ortho_qr(A,LDA,m,n)


    Orthogonalization using Q.R factorization
    
    A : matrix to orthogonalize
    
    LDA : leftmost dimension of A
    
    n : Number of rows of A
    
    m : Number of columns of A
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`davidson_diag_hjj_sjj`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgeqrf`
       * :c:func:`dorgqr`

 
.. c:function:: ortho_qr_unblocked:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine ortho_qr_unblocked(A,LDA,m,n)


    Orthogonalization using Q.R factorization
    
    A : matrix to orthogonalize
    
    LDA : leftmost dimension of A
    
    n : Number of rows of A
    
    m : Number of columns of A
    

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgeqr2`
       * :c:func:`dorg2r`

 
.. c:function:: overlap_gaussian_x:


    File : :file:`utils/one_e_integration.irp.f`

    .. code:: fortran

        double precision function overlap_gaussian_x(A_center,B_center,alpha,beta,power_A,power_B,dim)


    .. math::
    
     \sum_{-infty}^{+infty} (x-A_x)^ax (x-B_x)^bx exp(-alpha(x-A_x)^2) exp(-beta(x-B_X)^2) dx
    

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_explicit_poly_and_gaussian_x`

 
.. c:function:: overlap_x_abs:


    File : :file:`utils/one_e_integration.irp.f`

    .. code:: fortran

        subroutine overlap_x_abs(A_center,B_center,alpha,beta,power_A,power_B,overlap_x,lower_exp_val,dx,nx)


    .. math                      ::
    
     \int_{-infty}^{+infty} (x-A_center)^(power_A) * (x-B_center)^power_B * exp(-alpha(x-A_center)^2) * exp(-beta(x-B_center)^2) dx
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_overlap_abs`

 
.. c:function:: print_memory_usage:


    File : :file:`utils/memory.irp.f`

    .. code:: fortran

        subroutine print_memory_usage()


    Prints the memory usage in the output

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`write_time`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`resident_memory`
       * :c:func:`total_memory`

 
.. c:function:: quick_dsort:


    File : :file:`utils/sort.irp.f_template_261`

    .. code:: fortran

        subroutine quick_dsort(x, iorder, isize)


    Sort array x(isize) using the quicksort algorithm.
    iorder in input should be (1,2,3,...,isize), and in output
    contains the new order of the elements.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nproc`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`dsort`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`rec_d_quicksort`

 
.. c:function:: quick_i2sort:


    File : :file:`utils/sort.irp.f_template_261`

    .. code:: fortran

        subroutine quick_i2sort(x, iorder, isize)


    Sort array x(isize) using the quicksort algorithm.
    iorder in input should be (1,2,3,...,isize), and in output
    contains the new order of the elements.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nproc`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`i2sort`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`rec_i2_quicksort`

 
.. c:function:: quick_i8sort:


    File : :file:`utils/sort.irp.f_template_261`

    .. code:: fortran

        subroutine quick_i8sort(x, iorder, isize)


    Sort array x(isize) using the quicksort algorithm.
    iorder in input should be (1,2,3,...,isize), and in output
    contains the new order of the elements.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nproc`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`i8sort`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`rec_i8_quicksort`

 
.. c:function:: quick_isort:


    File : :file:`utils/sort.irp.f_template_261`

    .. code:: fortran

        subroutine quick_isort(x, iorder, isize)


    Sort array x(isize) using the quicksort algorithm.
    iorder in input should be (1,2,3,...,isize), and in output
    contains the new order of the elements.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nproc`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`isort`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`rec_i_quicksort`

 
.. c:function:: quick_sort:


    File : :file:`utils/sort.irp.f_template_261`

    .. code:: fortran

        subroutine quick_sort(x, iorder, isize)


    Sort array x(isize) using the quicksort algorithm.
    iorder in input should be (1,2,3,...,isize), and in output
    contains the new order of the elements.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nproc`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`sort`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`rec__quicksort`

 
.. c:function:: recentered_poly2:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        subroutine recentered_poly2(P_new,x_A,x_P,a,P_new2,x_B,x_Q,b)


    Recenter two polynomials

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`binom`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`give_explicit_poly_and_gaussian`
       * :c:func:`give_explicit_poly_and_gaussian_x`

 
.. c:function:: resident_memory:


    File : :file:`utils/memory.irp.f`

    .. code:: fortran

        subroutine resident_memory(value)


    Returns the current used memory in gigabytes used by the current process.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`check_mem`
       * :c:func:`davidson_diag_hjj_sjj`
       * :c:func:`print_memory_usage`
       * :c:func:`zmq_pt2`

 
.. c:function:: rint:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        double precision function rint(n,rho)


    .. math::
    
      \int_0^1 dx \exp(-p x^2) x^n
    

 
.. c:function:: rint1:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        double precision function rint1(n,rho)


    Standard version of rint

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`inv_int`
       * :c:data:`fact_inv`

 
.. c:function:: rint_large_n:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        double precision function rint_large_n(n,rho)


    Version of rint for large values of n

 
.. c:function:: rint_sum:


    File : :file:`utils/integration.irp.f`

    .. code:: fortran

        double precision function rint_sum(n_pt_out,rho,d1)


    Needed for the calculation of two-electron integrals.

 
.. c:function:: rinteg:


    File : :file:`utils/need.irp.f`

    .. code:: fortran

        double precision function rinteg(n,u)



 
.. c:function:: rintgauss:


    File : :file:`utils/need.irp.f`

    .. code:: fortran

        double precision function rintgauss(n)



 
.. c:function:: sabpartial:


    File : :file:`utils/need.irp.f`

    .. code:: fortran

        double precision function SABpartial(zA,zB,A,B,nA,nB,gamA,gamB,l)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`binom`

 
.. c:function:: set_order:


    File : :file:`utils/sort.irp.f_template_347`

    .. code:: fortran

        subroutine set_order(x,iorder,isize)


    array A has already been sorted, and iorder has contains the new order of
    elements of A. This subroutine changes the order of x to match the new order of A.

 
.. c:function:: set_order_big:


    File : :file:`utils/sort.irp.f_template_412`

    .. code:: fortran

        subroutine set_order_big(x,iorder,isize)


    array A has already been sorted, and iorder has contains the new order of
    elements of A. This subroutine changes the order of x to match the new order of A.
    This is a version for very large arrays where the indices need
    to be in integer*8 format

 
.. c:function:: sort:


    File : :file:`utils/sort.irp.f_template_293`

    .. code:: fortran

        subroutine sort(x,iorder,isize)


    Sort array x(isize).
    iorder in input should be (1,2,3,...,isize), and in output
    contains the new order of the elements.

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`insertion_sort`
       * :c:func:`quick_sort`

 
.. c:function:: sorted_dnumber:


    File : :file:`utils/sort.irp.f_template_261`

    .. code:: fortran

        subroutine sorted_dnumber(x,isize,n)


    Returns the number of sorted elements

 
.. c:function:: sorted_i2number:


    File : :file:`utils/sort.irp.f_template_261`

    .. code:: fortran

        subroutine sorted_i2number(x,isize,n)


    Returns the number of sorted elements

 
.. c:function:: sorted_i8number:


    File : :file:`utils/sort.irp.f_template_261`

    .. code:: fortran

        subroutine sorted_i8number(x,isize,n)


    Returns the number of sorted elements

 
.. c:function:: sorted_inumber:


    File : :file:`utils/sort.irp.f_template_261`

    .. code:: fortran

        subroutine sorted_inumber(x,isize,n)


    Returns the number of sorted elements

 
.. c:function:: sorted_number:


    File : :file:`utils/sort.irp.f_template_261`

    .. code:: fortran

        subroutine sorted_number(x,isize,n)


    Returns the number of sorted elements

 
.. c:function:: svd:


    File : :file:`utils/linear_algebra.irp.f`

    .. code:: fortran

        subroutine svd(A,LDA,U,LDU,D,Vt,LDVt,m,n)


    Compute A = U.D.Vt
    
    LDx : leftmost dimension of x
    
    Dimsneion of A is m x n
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`mo_as_svd_vectors_of_mo_matrix`
       * :c:func:`mo_as_svd_vectors_of_mo_matrix_eig`
       * :c:func:`ortho_canonical`
       * :c:func:`ortho_lowdin`
       * :c:data:`s_half`
       * :c:data:`s_half_inv`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgesvd`

 
.. c:function:: total_memory:


    File : :file:`utils/memory.irp.f`

    .. code:: fortran

        subroutine total_memory(value)


    Returns the current used memory in gigabytes used by the current process.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`print_memory_usage`

 
.. c:function:: u_dot_u:


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        double precision function u_dot_u(u,sze)


    Compute <u|u>

 
.. c:function:: u_dot_v:


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        double precision function u_dot_v(u,v,sze)


    Compute <u|v>

 
.. c:function:: wall_time:


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        subroutine wall_time(t)


    The equivalent of cpu_time, but for the wall time.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`add_integrals_to_map`
       * :c:func:`add_integrals_to_map_erf`
       * :c:func:`add_integrals_to_map_no_exit_34`
       * :c:func:`add_integrals_to_map_three_indices`
       * :c:data:`ao_pseudo_integrals_local`
       * :c:data:`ao_pseudo_integrals_non_local`
       * :c:data:`ao_two_e_integrals_erf_in_map`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:func:`davidson_converged`
       * :c:func:`davidson_diag_hjj_sjj`
       * :c:data:`output_wall_time_0`
       * :c:func:`pt2_collector`
       * :c:func:`run_pt2_slave`
       * :c:func:`run_slave_main`
       * :c:func:`write_time`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`system_clock`

 
.. c:function:: write_git_log:


    File : :file:`utils/util.irp.f`

    .. code:: fortran

        subroutine write_git_log(iunit)


    Write the last git commit in file iunit.

