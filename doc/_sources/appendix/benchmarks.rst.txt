Benchmarks
==========

The determinant selection, MR-PT2 and diagonalization are parallelized with
distributed parallelism. Benchmarks for the [NH2-CH-NH2]+ molecule in the
aug-cc-pVDZ basis set are presented with up to 50 nodes (1800 cores) on
CALMIP's `Olympe`_ supercomputer, and 200 nodes (9600 cores) on GENCI's
`Irene`_ supercomputer. This represents an active space of 18 electrons
in 111 MOs.

- Nodes of Olympe have two Skylake sockets, 2x18 cores @ 2.3GHz.
- Nodes of Irene  have two Skylake sockets, 2x24 cores @ 2.7GHz.

Convergence of the energy
-------------------------

.. figure:: /_static/cn3_energy.png
   :alt: Convergence of the energy.

   Convergence of the variational energy, with and without the PT2 correction.
   Both energies converge to the (frozen core) FCI energy. 
   The plot is displayed for the ground state and for the 1st excited state.


Variational energy
^^^^^^^^^^^^^^^^^^

================ ================ ================ ===============
 Number of dets    Ground state     Excited state  Excitation (eV)
================ ================ ================ ===============
              7    -149.489 186     -149.207 354     7.67  
            123    -149.536 265     -149.261 860     7.47  
          3 083    -149.685 606     -149.404 450     7.65  
         29 409    -149.826 151     -149.547 275     7.59  
        168 595    -149.900 352     -149.626 058     7.46  
      1 322 537    -149.946 655     -149.675 032     7.39  
      8 495 334    -149.972 032     -149.704 145     7.29  
      9 356 952    -149.973 375     -149.706 822     7.25  
     42 779 636    -149.987 370     -149.721 470     7.24  
    186 978 487    -149.998 582     -149.733 039     7.23  
================ ================ ================ ===============


Variational energy + PT2 correction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

================ ================ ================ ===============
 Number of dets    Ground state     Excited state  Excitation (eV)
================ ================ ================ ===============
             7    -150.161 107     -149.904 883      6.97  
           123    -150.116 958     -149.849 465      7.28  
         3 083    -150.043 5(2)    -149.780 8(2)     7.15  
        29 409    -150.022 2(2)    -149.758 3(2)     7.18  
       168 595    -150.019 9(1)    -149.754 5(1)     7.22  
     1 322 537    -150.017 89(7)   -149.752 55(7)    7.22  
     8 495 334    -150.015 97(4)   -149.750 87(5)    7.21  
     9 356 952    -150.015 89(3)   -149.750 66(3)    7.22  
    42 959 496    -150.016 75(2)   -149.751 88(2)    7.21  
   186 978 487    -150.017 51(2)   -149.752 90(2)    7.20  
================ ================ ================ ===============


Davidson's diagonalization
--------------------------

We present the parallel speedup curve, and the wall-clock time in seconds
required to compute one iteration for two wave functions measured on Olympe
and Irene.

.. figure:: /_static/speedup_davidson.png
   :alt: Parallel speedup of Davidson's diagonalization.

   Parallel speedup of Davidson's diagonalization measured on Olympe and Irene.

Olympe
^^^^^^

======================= ====================== =======================
Number of 36-core Nodes 9 356 952 determinants 42 959 496 determinants
======================= ====================== =======================
               1                775.55               11 198.70
               5                169.88                2 288.58
              10                 93.22                1 213.95
              20                 56.86                  626.41
              30                 43.76                  445.65
              40                 36.18                  350.25
              50                 33.67                  295.25
======================= ====================== =======================


Irene 
^^^^^

======================= ====================== =======================
Number of 48-core Nodes 9 356 952 determinants 42 959 496 determinants
======================= ====================== =======================
               1               572.98                 9 154.30 *
              10                72.55                   922.07
              25                38.88                   412.34
              50                27.95                   241.35
              75                27.54                   183.63
             100                27.86                   165.68
             150                28.14                   134.05
             200                27.77                   134.64
======================= ====================== =======================


PT2 correction
--------------

We present the parallel speedup curve, and the wall-clock time in seconds
required to compute the PT2 correction for two wave functions measured on
Olympe and Irene.


.. figure:: /_static/speedup_pt2.png
   :alt: Parallel speedup of the PT2 computation of the ground state.

   Parallel speedup of the PT2 computation of the ground state measured
   on Olympe and Irene.


Olympe
^^^^^^

======================= ====================== =======================
Number of 36-core Nodes   Ground state (9.3M)   Excited state (9.3M)
======================= ====================== =======================
               1               7 883.74            9 829.19 
               5               1 629.06            2 022.36 
              10                 832.89            1 029.91 
              20                 440.76              537.37 
              30                 303.31              378.69 
              40                 246.12              296.31 
              50                 201.84              241.55 
======================= ====================== =======================


Irene 
^^^^^

======================= ====================== ======================= ====================== =======================
Number of 48-core Nodes   Ground state (9.3M)   Excited state (9.3M)    Ground state (42.9M)   Excited state (42.9M)
======================= ====================== ======================= ====================== =======================
               1               4 935.81                6 152.29               24 586.62             37 440.59  
              10                 525.95                  652.23                2 458.66              3 086.19  
              25                 237.47                  286.06                1 041.69              1 295.43  
              50                 144.39                  174.12                  588.35                724.25  
              75                 109.13                  129.17                  446.74                537.59  
             100                 100.75                  103.43                  367.21                450.32  
             150                  82.04                   91.77                  298.63                358.25  
             200                  75.62                   85.25                  268.96                312.23  
======================= ====================== ======================= ====================== =======================





