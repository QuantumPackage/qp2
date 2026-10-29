===================
Developing plugins
===================

Once :file:`quantum_package.rc` has been sourced (or you are in the
|QPSH|), the plugin commands are run with ``qp plugins ...``; the
underlying script is :file:`bin/qp_plugins`.

Creating a repository of plugins
--------------------------------

The purpose of :file:`$QP_ROOT/plugins` is to contain local copies of
external repositories of plugins.

Create a repository, for example :file:`qp_plugins_user`, hosted somewhere
(GitLab, GitHub, etc...), and clone it in the :file:`$QP_ROOT/plugins`
directory. Instead of cloning by hand, you can also run::

    qp plugins download <url>

where ``<url>`` points to a git repository or to a :file:`.tar.gz`
archive. ``qp plugins update`` runs ``git pull`` in every repository of
:file:`$QP_ROOT/plugins`, and ``qp plugins list -q`` prints their names.


Creating a new plugin
---------------------

To create a new plugin named :file:`my_plugin` in this repository, run::

    qp plugins create -n my_plugin -r qp_plugins_user

(the ``-r`` option is optional; without it the plugin is created in
:file:`$QP_ROOT/plugins/local`).

``create`` fills the directory :file:`$QP_ROOT/plugins/qp_plugins_user/my_plugin`
with a :file:`NEED` file, a :file:`README.rst`, and a
:file:`my_plugin.irp.f` program that prints "Hello world", and installs
the plugin right away. "Installing" means that
:file:`$QP_ROOT/src/my_plugin` becomes a symbolic link to the plugin, and
that the modules listed in its :file:`NEED` file that are not already in
:file:`$QP_ROOT/src` are symbolic-linked there too. The plugin can then
be compiled like any other module::

    cd $QP_ROOT/plugins/qp_plugins_user/my_plugin
    ninja

After editing the :file:`NEED` file, reinstall the plugin to update the
symbolic links::

    qp plugins install my_plugin


Useful commands
---------------

::

    qp plugins list         list all the plugins of the repositories
    qp plugins list -i      only the installed ones (symbolic-linked in $QP_ROOT/src)
    qp plugins list -u      only the uninstalled ones
    qp plugins list -q      the repository names
    qp plugins uninstall my_plugin

``uninstall`` removes the symbolic links (and the generated files) but
leaves the plugin in its repository, so you can reinstall it later.
