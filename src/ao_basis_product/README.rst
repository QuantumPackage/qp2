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
