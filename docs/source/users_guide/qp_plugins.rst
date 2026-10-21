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
    qp_plugins download <url> [-n <name>]
    qp_plugins install <name>...
    qp_plugins uninstall <name>
    qp_plugins remove <name>
    qp_plugins update [-r <repo>]
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

.. option:: -n, --name=<name>

    With ``download``, give the name of the local directory created by
    ``git clone`` instead of using the name of the repository. With
    ``create``, create a new plugin named ``name`` (in the local
    repository by default).

.. option:: install <plugin_name>

    Install the plugin ``plugin_name``.

.. option:: uninstall <plugin_name>

    Uninstall the plugin ``plugin_name``.

.. option:: remove <plugin_name>

    Same as ``uninstall``.

.. option:: update

    Update the repositories of the plugins by running ``git pull`` in
    each of them. Should be followed by a re-compilation.

.. option:: -r, --repository=<repo>

    With ``update``, restrict the update to the repository ``repo``. With
    ``create``, specify in which repository the new plugin will be created.



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

