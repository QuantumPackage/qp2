.. _qp_plugins:

==========
qp_plugins
==========

.. program:: qp_plugins

This command deals with all external plugins of |qp|. Plugin
repositories can be downloaded, and the plugins in these repositories
can be installed/uninstalled or created.

Usage
-----

.. code:: bash

    qp_plugins list [-i] [-u] [-q]
    qp_plugins download <url>
    qp_plugins install <name>...
    qp_plugins uninstall <name>
    qp_plugins create -n <name> [-r <repo>] [<needed_modules>...]

.. option:: list

    List all the available plugins.

.. option:: -i, --installed 

    List all the *installed* plugins.

.. option:: -u, --uninstalled

    List all the *uninstalled* plugins.

.. option:: -q, --repositories

    List all the downloaded repositories.

.. option:: download <url>

    Download an external repository. The URL points to a tar.gz file or a
    git repository, for example:

    * http://example.com/site/example.tar.gz
    * git@gitlab.com:user/example_repository

.. option:: install <plugin_name>

    Install the plugin ``plugin_name``.

.. option:: uninstall <plugin_name>

    Uninstall the plugin ``plugin_name``.

.. option:: update 

    Update the repositories of the plugins. Should be followed by a re-compilation.

.. option:: -n, --name=<plugin_name>

    Create a new plugin named ``plugin_name`` (in local repository by default).

.. option:: -r, --repository=<repo>

    Specify in which repository the new plugin will be created.



Example
-------

Let us download, install and compile some specific external plugins from
`<https://gitlab.com/eginer/qp_plugins_eginer>`_ .

First, download the git repo associated to these plugins. To do so,
first go to the `plugins` directory in the |QP| and execute:

.. code:: bash

    qp_plugins download https://gitlab.com/eginer/qp_plugins_eginer


This will create in the directory `plugins` a local copy of
the git repo located at the URL you indicated. Then, go in
`qp_plugins_eginer/stable/`

.. code:: bash

    cd qp_plugins_eginer/stable/

In the directory `stable`, there are many directories which all
correspond to a specific plugin that have been developed by the person
in charge of the repository. All these plugins might use some global
variables and routines contained in the core modules of the |QP|.

Now let us install the plugin `rsdft_cipsi`: 

.. code:: bash

    qp_plugins install rsdft_cipsi

This will link this directory to the |QP| which means that when the code
will be compiled, this plugin will be compiled to and therefore all the
executables/scripts/input keywords contained in this module will be
available as if there were part of the core of the |QP|.

Then, to compile the new plugin, just recompile the |QP| as usual by
going at the root of the |QP| directory:

.. code:: bash

    cd $QP_ROOT
    ninja 

Finally, if you go back to the plugin directory you just installed, you
should see all the executables/scripts which have been created and which
are now available with the `qp_run` command.

