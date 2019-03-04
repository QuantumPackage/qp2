=================
Quick-start guide
=================

This tutorial should teach you everything you need to get started with
the the basics of the |qp|. As an example, we will run a frozen core
|CIPSI| calculation on the HCN molecule in the 631-G basis set.


Demo video
==========

This tutorial can be directly watched at: 


`<https://www.youtube.com/watch?v=4nmdCAPkZlc>`_


Hands on
========

.. important::

   Before using the |qp|, it is required to load the environment variables 
   relatives to the |QP| or to be in the |qpsh| mode. 


Please execute in the current shell: 

.. code:: bash 

  ${QP_ROOT}/bin/qpsh 

where :code:`${QP_ROOT}` is the path to the source files of the |QP| installed on your architecture.  

The |QPSH| mode: a bash-like experience for quantum chemistry
-------------------------------------------------------------

The |QP| has been designed pretty much as an *interactive* environment for quantum-chemistry calculations, 
in order to facilitate the user experience. 

Just like in bash, there are many commands in the |QP| (see for instance :ref:`qp_edit` or :ref:`qp_run`) 
which help in handling useful data or running executables (see for instance :ref:`scf` or :ref:`fci`). 

All commands designed within the |qp| **begins** with `qp`, and are two ways of running a **command**: 

* the *executable* associated to the command: 

.. code:: bash

  qp_command 

or the *qp* command which calls the *executable* :code:`qp_command`: 

.. code:: bash

  qp command 

Usually, when using the :command:`qp` command, the name of the |EZFIO| database is omitted.

The advantage or using :code:`qp command` is that you can, just like in bash, have: 

* the :kbd:`Tab` key for the auto-completion for basically any command of the |QP| 

* man pages with -h, --help or qp man 


Just try, for instance: 

.. code:: bash 

  qp 

and then use the auto-completion. You will show appear all possible commands that you can run: 


.. code:: bash 

  convert_output_to_ezfio  -h                       plugins                  unset_file               
  create_ezfio             man                      set_file                 update            

Then, try, still with the auto-completion, 

.. code:: bash

  qp create

You will see appear all the options for the :ref:`qp_create_ezfio` commands. 


Create the EZFIO database
-------------------------

The data relative to calculations are stored in an |EZFIO| database.
|EZFIO| is a hierarchical data format which uses the hierarchy of the
file system to organize the data, as files stored in a directory. The
data in the |EZFIO| directory are stored as plain text files, so it can
be opened with any text editor.
To access the data of the |EZFIO| database, the APIs (Fortran, |Python|,
|OCaml| or Bash) provided by |EZFIO| should be used, or tools using
these APIs such as :ref:`qp_edit` provided with the |qp|.

First, create an `xyz` file containing the coordinates of the molecule.
The file :file:`hcn.xyz` contains::

   3
   HCN molecule
   C    0.0    0.0    0.0
   H    0.0    0.0    1.064
   N    0.0    0.0    -1.156


This xyz file is now used with the :ref:`qp_create_ezfio` command to
create an |EZFIO| database with the 6-31G basis set:

.. code:: bash

  qp create_ezfio -b "6-31G" hcn.xyz -o hcn

The EZFIO database now contains data relative to the nuclear coordinates
and the atomic basis set:

.. code:: bash

  $ ls hcn
  ao_basis           becke_numerical_grid  dft_keywords  mo_one_e_ints      perturbation
  ao_one_e_ints      davidson              dressing      mo_two_e_erf_ints  pseudo
  ao_two_e_erf_ints  density_for_dft       electrons     mo_two_e_ints      scf_utils
  ao_two_e_ints      determinants          ezfio         nuclei             work


Run a Hartree-Fock calculation
------------------------------

The program :ref:`qp_run` is the driver program of the |qp|. To run a
|scf| calculation, just run

.. code:: bash

    qp run scf

The expected energy is ``-92.827856698`` au.

.. seealso:: 

    The documentation of the :ref:`module_hartree_fock` module and that of the
    :ref:`scf` program.

This creates the |MOs| in the |EZFIO| database that will be used to
perform any other post-SCF method. The |qp| does not handle symmetry and
the |MOs| are stored by increasing order of Fock energies.

Choose the target |MO| space
----------------------------

Now, we will modify the |EZFIO| database to make |CIPSI| calculation only in the
full set of valence |MOs|, keeping the core |MOs| frozen. The simple
command :ref:`qp_set_frozen_core` does this automatically:

.. code:: bash

    qp set_frozen_core


The general command to specify core and active orbitals is :ref:`qp_set_mo_class`. 
In the case of HCN molecule in the 631G basis, one has 20 |MOs| in total and the two first orbitals are frozen:

.. code::

    qp set_mo_class --core "[1-2]" --act "[3-20]"



Run the |CIPSI| calculation
----------------------------

We will now use the |CIPSI| algorithm to estimate the |FCI| energy.

.. code::

    qp run fci | tee hcn.fci.out 


The program will start with a single determinant and will iteratively:

* Select the most important determinants from the external space and add them to the
  internal space
* Add all the necessary determinants to allow the eigenvector of |H| to be
  also an eigenstate of |S^2|
* Diagonalize |H| in the enlarged internal space
* Compute (stochastically) the second-order perturbative contribution to the energy 
* Extrapolate the variational energy by fitting
  :math:`E=E_\text{FCI} - \alpha\, E_\text{PT2}`

By default, the program will stop when more than one million determinants have
entered in the internal space, or when the |PT2| energy is below :math:`10^{-4}`.

To have a pictural illustration of the convergence of the |CIPSI| algorithm, just run 

.. code::

    qp_e_conv_fci

This will create the files "hcn.fci.out.conv" containing the data of the convergence of the energy that can be plotted, together with the file "hcn.fci.out.conv.1.eps" which is obtained from the gnuplot plot file "hcn.fci.out.conv.plt". 


The estimated |FCI| energy of HCN is ``-93.0501`` au.

.. seealso:: 

    The documentation of the :ref:`module_fci` module and that of the :ref:`fci` program.


