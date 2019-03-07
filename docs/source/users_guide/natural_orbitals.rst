Natural orbitals
================

Summary
-------

To produce state-average natural orbitals, run ::

    qp_run save_natorb file.ezfio

The |MOs| will be replaced, so the two-electron integrals and the wave
function are invalidated as well.



Extracting natural orbitals
---------------------------

Once obtained the near |FCI| wave function, one can obtain many         
quantities related to it. One of these quantities are the natural       
orbitals which have the property of diagonalizing the one-body      
density matrix:                                                       

   .. math::

       \rho_{ij} = \delta_{ij}

where the element of the one-body density matrix :math:`\rho_{ij}` is
defined as:


   .. math::

       \rho_{ij} = \langle \Psi | \left( a^{\dagger}_{j,\alpha} a_{i,\alpha} + a^{\dagger}_{j,\beta} a_{i,\beta} \right) | \Psi \rangle


These orbitals are in general known to be better than the usual |HF|
|MOs| as they are obtained from a correlated wave function. To use these
orbitals for future calculations, one has to replace the current |MOs|
by the natural orbitals. To do so, just run:

.. code::

    qp_run save_natorb file.ezfio


Hands on
--------

.. important::

   As the |MOs| are changed, for the sake of coherence of future        
   calculations, the :ref:`save_natorb` program *automatically removes the     
   current wave function* stored in the |EZFIO| database and replaces   
   it by a single Slater determinant corresponding to a |HF| occupation 
   of the new spin orbitals. Also, all the keywords to read the one-    
   and two-electron integrals on the |MO| basis are set to ``None`` in  
   order to be sure to avoid reading integrals incompatible with the    
   current set of |MOs|.                                                

.. seealso:: 

    The documentation of the :ref:`save_natorb` program. 
 
