.. _qp_tunnel:

=========
qp_tunnel
=========

.. TODO

.. program:: qp_tunnel

Establishes a tunnel to allow communications between machines within
different networks, for example multiple MPI slave jobs running on
different clusters.


Usage
-----

.. code:: bash

  qp_tunnel [-g] (ADDRESS|EZFIO_DIR)

``EZFIO_DIR`` is the name of the |EZFIO| directory containing the data,
and ``ADDRESS`` is the address of another tunnel.


.. option:: -h, --help

   Displays the help message


.. option:: -g, --get-input

   Download the EZFIO directory from the remote instance of qp_tunnel.


Example
-------

.. code:: text

    +-------------------+    +------------------+
    |                   |    |                  |
    |   N1_1 N1_2 N1_3  |    |  N2_1 N2_2 N2_3  |
    |    |    |    |    |    |   |    |    |    |
    |    +----+----+    |    |   +----+----+    |
    |         |         |    |        |         |
    | C1      F1        |    |        F2     C2 |
    |         +---------=----=--------+         |
    |                   |    |                  |
    +-------------------+    +------------------+


Imagine you have two clusters, C1 and C2. Each cluster is accessible via SSH
on a front-end named respectively F1 and F2. Groups of nodes N1 and N2 have
been reserved by the batch scheduling system on both clusters.
Each node in N1 is on the same network as the other nodes of N1, but they
can't access the network on which the nodes of N2 are.

1) Start a parallel simulation on the cluster C1, running on nodes N1.
   We assume that there is a shared file system, such that F1 can access
   the EZFIO directory. We also assume that F1 can communicate with the
   nodes of N1.

2) Run a tunnel on the front-end F1 and keep it running:

.. code:: bash

    me@f1 $ qp_tunnel my_directory.ezfio
    Connect to:
    tcp://31.122.230.47:42379
    Ready   

3) On the front-end F2, run another instance connecting to the other one,
   which will fetch the |EZFIO| directory:

.. code:: bash

    me@f2 $ qp_tunnel --get-input tcp://31.122.230.47:42379
    Connect to:
    tcp://31.122.209.139:42379
    Communication [ OK ]
    Getting input... my_directory.ezfio ...done
    Ready

4) Keep the tunnel running, and you can now run a slave simulation within the
   nodes N2.


