=============================
Working with external plugins
=============================


|qp| has very few executables out of the box. Most of the time,
external plugins need to be downloaded and installed in the
:file:`$QP_ROOT/plugins` directory.

Plugins are usually hosted in external repositories. To download a
plugin, the remote repository needs to be downloaded, and the plugins of
the repository can be selected for installation.

To download an external repository of plugins, run the following
command:

.. code-block:: bash

    qp_plugins download http://somewhere/over/the/rainbow/ext_repo


This downloads a copy of the repository of external plugins :file:`ext_repo`
in :file:`$QP_ROOT/plugins`.

The list of available uninstalled plugins can be seen using:

.. code-block:: bash

    qp_plugins list -u


Now, the specific plugin :file:`ext_module` contained in the repository
:file:`ext_repo` can be installed using:

.. code-block:: bash

        qp_plugins install ext_module


The module is now accessible via a symbolic link in :file:`$QP_ROOT/src`,
and can be compiled as any module, running |Ninja|.


To remove the module, run

.. code-block:: bash

        qp_plugins uninstall ext_module


.. seealso:: 

     For a more detailed explanation and an example, see :ref:`qp_plugins`. 

