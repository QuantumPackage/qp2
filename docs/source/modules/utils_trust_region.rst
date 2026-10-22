.. _module_utils_trust_region: 
 
.. program:: utils_trust_region 
 
.. default-role:: option 
 
utils_trust_region
==================

No README.rst in this module; the providers and
subroutines below are documented from the irpf90 man pages. 
 
 
EZFIO parameters 
---------------- 
 
.. option:: thresh_delta
 
    Threshold to stop the optimization if the radius of the trust region delta < thresh_delta
 
    Default: 1.e-10
 
.. option:: thresh_rho
 
    Threshold for the step acceptance in the trust region algorithm, if (rho .geq. thresh_rho) the step is accepted, else the step is cancelled and a smaller step is tried until (rho .geq. thresh_rho)
 
    Default: 0.1
 
.. option:: thresh_eig
 
    Threshold to consider when an eigenvalue is 0 in the trust region algorithm
 
    Default: 1.e-12
 
.. option:: thresh_model
 
    If if ABS(criterion - criterion_model) < thresh_model, the program exit the trust region algorithm
 
    Default: 1.e-12
 
.. option:: absolute_eig
 
    If True, the algorithm replace the eigenvalues of the hessian by their absolute value to compute the step (in the trust region)
 
    Default: false
 
.. option:: thresh_wtg
 
    Threshold in the trust region algorithm to considere when the dot product of the eigenvector W by the gradient v_grad is equal to 0. Must be smaller than thresh_eig by several order of magnitude to avoid numerical problem. If the research of the optimal lambda cannot reach the condition (||x|| .eq. delta) because (||x|| .lt. delta), the reason might be that thresh_wtg is too big or/and thresh_eig is too small
 
    Default: 1.e-6
 
.. option:: thresh_wtg2
 
    Threshold in the trust region algorithm to considere when the dot product of the eigenvector W by the gradient v_grad is 0 in the case of avoid_saddle .eq. true. There is no particular reason to put a different value that thresh_wtg, but it can be useful one day
 
    Default: 1.e-6
 
.. option:: avoid_saddle
 
    Test to avoid saddle point, active if true
 
    Default: false
 
.. option:: version_avoid_saddle
 
    cf. trust region, not stable
 
    Default: 3
 
.. option:: thresh_rho_2
 
    Threshold for the step acceptance for the research of lambda in the trust region algorithm, if (rho_2 .geq. thresh_rho_2) the step is accepted, else the step is rejected
 
    Default: 0.1
 
.. option:: thresh_cc
 
    Threshold to stop the research of the optimal lambda in the trust region algorithm when (dabs(1d0-||x||^2/delta^2) < thresh_cc)
 
    Default: 1.e-6
 
.. option:: thresh_model_2
 
    if (ABS(criterion - criterion_model) < thresh_model_2), i.e., the difference between the actual criterion and the predicted next criterion, during the research of the optimal lambda in the trust region algorithm it prints a warning
 
    Default: 1.e-12
 
.. option:: version_lambda_search
 
    Research of the optimal lambda in the trust region algorithm to constrain the norm of the step by solving: 1 -> ||x||^2 - delta^2 .eq. 0, 2 ->  1/||x||^2 - 1/delta^2 .eq. 0
 
    Default: 2
 
.. option:: nb_it_max_lambda
 
    Maximal number of iterations for the research of the optimal lambda in the trust region algorithm
 
    Default: 100
 
.. option:: nb_it_max_pre_search
 
    Maximal number of iterations for the pre-research of the optimal lambda in the trust region algorithm
 
    Default: 40
 
 
Subroutines / functions 
----------------------- 
 
.. c:macro:: apply_mo_rotation:


    File : :file:`utils_trust_region/apply_mo_rotation.irp.f`

    .. code:: fortran

        subroutine apply_mo_rotation(R,prev_mos)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_coef`
       * :c:data:`mo_label`
       * :c:data:`mo_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_orb_opt_trust_v2`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`ezfio_set_determinants_mo_label`
       * :c:func:`wall_time`

 
.. c:macro:: d1_norm_inverse_trust_region:


    File : :file:`utils_trust_region/trust_region_optimal_lambda.irp.f`

    .. code:: fortran

        function d1_norm_inverse_trust_region(n,e_val,w,v_grad,lambda,delta)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`thresh_eig`

 
.. c:macro:: d1_norm_inverse_trust_region_omp:


    File : :file:`utils_trust_region/trust_region_optimal_lambda.irp.f`

    .. code:: fortran

        function d1_norm_inverse_trust_region_omp(n,e_val,tmp_wtg,lambda,delta)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`thresh_eig`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`omp_set_max_active_levels`

 
.. c:macro:: d1_norm_trust_region:


    File : :file:`utils_trust_region/trust_region_optimal_lambda.irp.f`

    .. code:: fortran

        function d1_norm_trust_region(n,e_val,w,v_grad,lambda,delta)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`thresh_eig`

 
.. c:macro:: d1_norm_trust_region_omp:


    File : :file:`utils_trust_region/trust_region_optimal_lambda.irp.f`

    .. code:: fortran

        function d1_norm_trust_region_omp(n,e_val,tmp_wtg,lambda,delta)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`thresh_eig`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`omp_set_max_active_levels`

 
.. c:macro:: d2_norm_inverse_trust_region:


    File : :file:`utils_trust_region/trust_region_optimal_lambda.irp.f`

    .. code:: fortran

        function d2_norm_inverse_trust_region(n,e_val,w,v_grad,lambda,delta)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`thresh_eig`

 
.. c:macro:: d2_norm_inverse_trust_region_omp:


    File : :file:`utils_trust_region/trust_region_optimal_lambda.irp.f`

    .. code:: fortran

        function d2_norm_inverse_trust_region_omp(n,e_val,tmp_wtg,lambda,delta)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`thresh_eig`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`omp_set_max_active_levels`

 
.. c:macro:: d2_norm_trust_region:


    File : :file:`utils_trust_region/trust_region_optimal_lambda.irp.f`

    .. code:: fortran

        function d2_norm_trust_region(n,e_val,w,v_grad,lambda,delta)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`thresh_eig`

 
.. c:macro:: d2_norm_trust_region_omp:


    File : :file:`utils_trust_region/trust_region_optimal_lambda.irp.f`

    .. code:: fortran

        function d2_norm_trust_region_omp(n,e_val,tmp_wtg,lambda,delta)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`thresh_eig`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`omp_set_max_active_levels`

 
.. c:macro:: f_norm_trust_region:


    File : :file:`utils_trust_region/trust_region_optimal_lambda.irp.f`

    .. code:: fortran

        function f_norm_trust_region(n,e_val,tmp_wtg,lambda)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`thresh_eig`

 
.. c:macro:: f_norm_trust_region_omp:


    File : :file:`utils_trust_region/trust_region_optimal_lambda.irp.f`

    .. code:: fortran

        function f_norm_trust_region_omp(n,e_val,tmp_wtg,lambda)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`thresh_eig`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`omp_set_max_active_levels`

 
.. c:macro:: factorial:


    File : :file:`utils_trust_region/rotation_matrix_iterative.irp.f`

    .. code:: fortran

        function factorial(n)



 
.. c:macro:: mat_to_vec_index:


    File : :file:`utils_trust_region/mat_to_vec_index.irp.f`

    .. code:: fortran

        subroutine mat_to_vec_index(p,q,i)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`vec_to_mat_v2`

 
.. c:macro:: rotation_matrix:


    File : :file:`utils_trust_region/rotation_matrix.irp.f`

    .. code:: fortran

        subroutine rotation_matrix(A,LDA,R,LDR,n,info,enforce_step_cancellation)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_orb_opt_trust_v2`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`abort`
       * :c:func:`dgemm`
       * :c:func:`dsyev`
       * :c:func:`wall_time`

 
.. c:macro:: rotation_matrix_iterative:


    File : :file:`utils_trust_region/rotation_matrix_iterative.irp.f`

    .. code:: fortran

        subroutine rotation_matrix_iterative(m,X,R)



    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`abort`
       * :c:func:`dgemm`
       * :c:func:`wall_time`

 
.. c:macro:: sub_to_full_rotation_matrix:


    File : :file:`utils_trust_region/sub_to_full_rotation_matrix.irp.f`

    .. code:: fortran

        subroutine sub_to_full_rotation_matrix(m,tmp_list,tmp_R,R)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_orb_opt_trust_v2`

 
.. c:macro:: trust_region_expected_e:


    File : :file:`utils_trust_region/trust_region_expected_e.irp.f`

    .. code:: fortran

        subroutine trust_region_expected_e(n,n2,v_grad,H,x,prev_energy,e_model)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`trust_region_step_w_expected_e`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemv`
       * :c:func:`wall_time`

 
.. c:macro:: trust_region_is_step_cancelled:


    File : :file:`utils_trust_region/algo_trust.irp.f`

    .. code:: fortran

        subroutine trust_region_is_step_cancelled(nb_iter,prev_criterion, criterion, criterion_model,rho,cancel_step)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`thresh_rho`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_orb_opt_trust_v2`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`trust_region_rho`

 
.. c:macro:: trust_region_optimal_lambda:


    File : :file:`utils_trust_region/trust_region_optimal_lambda.irp.f`

    .. code:: fortran

        subroutine trust_region_optimal_lambda(n,e_val,tmp_wtg,delta,lambda)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nb_it_max_lambda`
       * :c:data:`nb_it_max_pre_search`
       * :c:data:`thresh_cc`
       * :c:data:`thresh_model_2`
       * :c:data:`thresh_rho_2`
       * :c:data:`version_lambda_search`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`trust_region_step`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`wall_time`

 
.. c:macro:: trust_region_rho:


    File : :file:`utils_trust_region/trust_region_rho.irp.f`

    .. code:: fortran

        subroutine trust_region_rho(prev_energy, energy,e_model,rho)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`thresh_rho`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`trust_region_is_step_cancelled`

 
.. c:macro:: trust_region_step:


    File : :file:`utils_trust_region/trust_region_step.irp.f`

    .. code:: fortran

        subroutine trust_region_step(n,n2,nb_iter,v_grad,rho,e_val,w,x,delta)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`absolute_eig`
       * :c:data:`avoid_saddle`
       * :c:data:`thresh_eig`
       * :c:data:`thresh_wtg`
       * :c:data:`thresh_wtg2`
       * :c:data:`version_avoid_saddle`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`trust_region_step_w_expected_e`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`trust_region_optimal_lambda`
       * :c:func:`wall_time`

 
.. c:macro:: trust_region_step_w_expected_e:


    File : :file:`utils_trust_region/algo_trust.irp.f`

    .. code:: fortran

        subroutine trust_region_step_w_expected_e(n,n2,H,W,e_val,v_grad,prev_criterion,rho,nb_iter,delta,criterion_model,x,must_exit)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`thresh_delta`
       * :c:data:`thresh_model`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_orb_opt_trust_v2`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`trust_region_expected_e`
       * :c:func:`trust_region_step`

 
.. c:macro:: vec_to_mat_index:


    File : :file:`utils_trust_region/vec_to_mat_index.irp.f`

    .. code:: fortran

        subroutine vec_to_mat_index(i,p,q)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`diag_hessian_list_opt`
       * :c:func:`diag_hessian_opt`
       * :c:func:`first_diag_hessian_list_opt`
       * :c:func:`first_diag_hessian_opt`
       * :c:func:`first_gradient_list_opt`
       * :c:func:`first_gradient_opt`
       * :c:func:`first_hessian_list_opt`
       * :c:func:`first_hessian_opt`
       * :c:func:`gradient_list_opt`
       * :c:func:`gradient_opt`
       * :c:func:`hessian_list_opt`
       * :c:func:`hessian_opt`

 
.. c:macro:: vec_to_mat_v2:


    File : :file:`utils_trust_region/vec_to_mat_v2.irp.f`

    .. code:: fortran

        subroutine vec_to_mat_v2(n,m,v_x,m_x)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_orb_opt_trust_v2`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`mat_to_vec_index`

