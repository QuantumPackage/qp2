.. _module_dressing: 
 
.. program:: dressing 
 
.. default-role:: option 
 
=========
dress_zmq
=========

Module to facilitate the construction of modules using dressed
Hamiltonians, parallelized with |ZeroMQ|.

 
 
 
EZFIO parameters 
---------------- 
 
.. option:: thresh_dressed_ci
 
    Threshold on the convergence of the dressed |CI| energy
 
    Default: 1.e-5
 
.. option:: n_it_max_dressed_ci
 
    Maximum number of dressed |CI| iterations
 
    Default: 10
 
.. option:: dress_relative_error
 
    Stop stochastic dressing when the relative error is smaller than :option:`perturbation PT2_relative_error`
 
    Default: 0.001
