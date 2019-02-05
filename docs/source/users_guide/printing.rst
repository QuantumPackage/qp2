Printing the near-|FCI| wave function 
=====================================

Once obtained the near-|FCI| energy, one can also take a closer look at 
the wave function stored in the |EZFIO| database. If the wave function  
contains less than :math:`10^4` determinants, you can directly read it  
with the :ref:`qp_edit` command. Just run                               

.. code::

    qp_edit file.ezfio


.. important::

   The :ref:`qp_edit` mode virtually makes human-friendly the           
   architecture of the |EZFIO| database through the use of a
   the text editor defined by the :envvar:`EDITOR` environment
   variable.


Then, look for the word *hand* when you are in the :ref:`qp_edit`       
mode. If the research is negative, then it means that the wave          
function stored in the |EZFIO| database is too large to be edited       
interactively in :ref:`qp_edit` mode. An alternative is to use the      
:command:`print_wf` command:                                            

.. code::

    qp_run print_wf file.ezfio | tee file.ezfio.fci_natorb.wf

This program will, by default, print out the first :math:`10^4`
determinants whatever the size of the wave function stored in the
|EZFIO| directory. If you want to change the number of printed Slater
determinants, just change the :option:`determinants n_det_print_wf`
keyword using the :ref:`qp_edit` tool.

.. seealso:: 

    The documentation of the :ref:`print_wf` program. 


