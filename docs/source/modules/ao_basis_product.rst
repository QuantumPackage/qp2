.. _module_ao_basis_product: 
 
.. program:: ao_basis_product 
 
.. default-role:: option 
 
================
ao_basis_product
================

Basis of products of |AOs|.

This module builds the set of |AO| products used as a projector basis.
:c:data:`basis_prod_num` is the number of products (at most
:math:`N_\text{AO}(N_\text{AO}+1)/2`) and :c:data:`basis_prod_idx`
gives the indices of the two |AOs| of each product. A product is kept
when the overlap of its two |AOs| exceeds
:c:data:`ao_basis_prod_threshold`, which is set to the value of the
`ao_one_e_integrals_threshold` parameter.

:c:data:`basis_prod_overlap` is the overlap matrix of the product
basis and :c:data:`basis_prod_overlap_inv` its inverse. The overlap of
two products is computed by :c:func:`basis_prod_overlap_func`, which
reduces it to :c:func:`general_overlap_integral`, a general expression
for the overlap of two products of primitives.
 
 
 
Providers 
--------- 
 
.. c:var:: ao_basis_prod_threshold


    File : :file:`ao_basis_product/basis_prod.irp.f`

    .. code:: fortran

        double precision	:: ao_basis_prod_threshold	


    Threshold to select a product of basis functions

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_one_e_integrals_threshold`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`basis_prod_idx`
       * :c:data:`basis_prod_overlap`

 
.. c:var:: basis_prod_idx


    File : :file:`ao_basis_product/basis_prod.irp.f`

    .. code:: fortran

        integer, allocatable	:: basis_prod_idx	(2,basis_prod_num)


    Indices of main basis to create projector basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_basis_prod_threshold`
       * :c:data:`ao_num`
       * :c:data:`ao_overlap_abs`
       * :c:data:`basis_prod_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`basis_prod_overlap`
       * :c:data:`mu_of_r_projector_ao_prod`

 
.. c:var:: basis_prod_num


    File : :file:`ao_basis_product/basis_prod.irp.f`

    .. code:: fortran

        integer	:: basis_prod_num	


    Maximum number of basis functions for the projector

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`basis_prod_idx`
       * :c:data:`basis_prod_overlap`
       * :c:data:`mu_of_r_projector_ao_prod`

 
.. c:var:: basis_prod_overlap


    File : :file:`ao_basis_product/basis_prod.irp.f`

    .. code:: fortran

        double precision, allocatable	:: basis_prod_overlap	(basis_prod_num,basis_prod_num)
        double precision, allocatable	:: basis_prod_overlap_inv	(basis_prod_num,basis_prod_num)


    Overlap matrix of the basis of products

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_basis_prod_threshold`
       * :c:data:`ao_overlap_abs`
       * :c:data:`basis_prod_idx`
       * :c:data:`basis_prod_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mu_of_r_projector_ao_prod`

 
.. c:var:: basis_prod_overlap_inv


    File : :file:`ao_basis_product/basis_prod.irp.f`

    .. code:: fortran

        double precision, allocatable	:: basis_prod_overlap	(basis_prod_num,basis_prod_num)
        double precision, allocatable	:: basis_prod_overlap_inv	(basis_prod_num,basis_prod_num)


    Overlap matrix of the basis of products

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_basis_prod_threshold`
       * :c:data:`ao_overlap_abs`
       * :c:data:`basis_prod_idx`
       * :c:data:`basis_prod_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mu_of_r_projector_ao_prod`

 
 
Subroutines / functions 
----------------------- 
 
.. c:macro:: basis_prod_overlap_func:


    File : :file:`ao_basis_product/basis_prod.irp.f`

    .. code:: fortran

        double precision function basis_prod_overlap_func(a,b)


    Overlap matrix element of the basis of products

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_basis_prod_threshold`
       * :c:data:`ao_overlap_abs`
       * :c:data:`basis_prod_idx`

 
.. c:macro:: general_overlap_integral:


    File : :file:`ao_basis_product/basis_prod.irp.f`

    .. code:: fortran

        double precision function general_overlap_integral(dim,            &
     P_new,P_center,fact_p,p,p_inv,iorder_p,                       &
     Q_new,Q_center,fact_q,q,q_inv,iorder_q)


    Computes the overlap integral <p|q> where p,q are products of primitives

