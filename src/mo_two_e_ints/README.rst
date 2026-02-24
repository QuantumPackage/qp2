==================
mo_two_e_ints
==================

This module computes and manages molecular orbital (MO) two-electron integrals.

The module implements several algorithms and optimization strategies for handling
the large number of two-electron integrals that arise in quantum chemical calculations.

Key Features
------------

1. **Convention Awareness**: 
   - AO integrals use the chemists' convention: (ik|jl) = (11|22)
   - MO integrals use the physicist's convention: <ij|kl> = <12|12>

2. **Efficient Storage**: Two-electron integrals are stored in a sparse map structure
   to minimize memory usage, as most integrals are zero due to symmetry and
   selection rules.

3. **Specialized Dense Arrays**:
   - For integrals with 2 identical indices such as <ij|ik> and <ij|jk>,
     dedicated dense arrays are maintained to accelerate computations
     involving single excitations
   - These arrays significantly improve performance for CI calculations
     that frequently compute single excitation amplitudes

4. **Multi-Level Caching System**:
   - Implements a sophisticated cache mechanism for frequently accessed integrals
   - Cache size is configurable via ``mo_integrals_cache_shift`` parameter
   - Automatically manages memory usage while maximizing hit rates for common integrals
   - Reduces redundant computations and improves overall performance

5. **Multiple Computation Strategies**:
   - Direct AO-to-MO transformation using DGEMM operations
   - Sparse transformation methods for memory-efficient calculation
   - Cholesky decomposition for large-scale calculations

6. **Hardware Acceleration**:
   - GPU acceleration support for faster computation
   - Multi-threading with OpenMP
   - Optimized memory access patterns

7. **Flexible Configuration**:
   - Configurable integral thresholds for pruning negligible values
   - Adjustable cache sizes for optimal performance
   - Support for both single and double precision

8. **I/O Options**:
   - Read/write integrals from/to disk for memory-constrained environments
   - Integration with Cholesky decompositions for enhanced efficiency

Implementation Details
----------------------

The module provides several key routines:

- ``get_two_e_integral(i,j,k,l,mo_integrals_map)``: Fetches a specific MO integral
- ``mo_two_e_integral(i,j,k,l)``: Direct access to MO integrals
- ``mo_two_e_integrals_jj``, ``mo_two_e_integrals_jj_exchange``, ``mo_two_e_integrals_jj_anti``:
  Specialized routines for 2-index integrals
- ``big_array_exchange_integrals`` and ``big_array_coulomb_integrals``:
  3-index arrays for efficient exchange and Coulomb integrals computation

Storage Format
--------------

Integrals are stored in a sparse map data structure optimized for:
- Memory efficiency by storing only non-zero integrals
- Fast lookup of integrals by their 4-index notation
- Support for both symmetric and antisymmetric integral combinations
- Integration with multi-level caching for frequently accessed integrals

The caching system automatically manages memory usage while maximizing hit rates for common integrals,
reducing redundant computations and improving overall performance. The cache size can be adjusted via
the ``mo_integrals_cache_shift`` parameter in the EZFIO configuration.

The module automatically determines the most efficient computation strategy based on:
- Available memory (``qp_max_mem``)
- System architecture (CPU/GPU)
- Configuration settings (``do_mo_cholesky``, ``mo_integrals_threshold``)

Cholesky Decomposition
----------------------

When ``do_mo_cholesky`` is enabled, the module uses Cholesky decomposition to:
- Reduce memory footprint for large systems
- Enable faster integral evaluation in CI calculations
- Support efficient handling of the large number of integrals in correlated wavefunction methods

Usage Example
-------------

To use the two-electron integrals in your code:

.. code-block:: fortran

   ! Get a specific 4-index MO integral
   double precision :: integral_value
   integral_value = mo_two_e_integral(i, j, k, l)

   ! Or use the map-based approach for more control
   integral_value = get_two_e_integral(i, j, k, l, mo_integrals_map)

The module is automatically initialized during the configuration process and
provides the integrals as part of the standard quantum package environment.



