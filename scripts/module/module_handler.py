#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module utilitary

Usage:
    module_handler.py print_descendant      [<module_name>...]
    module_handler.py clean                 [ --all | <module_name>...]
    module_handler.py tidy                  [ --all | <module_name>...]
    module_handler.py create_git_ignore     [ --all | <module_name>...]

Options:
    print_descendant        Print the genealogy of the needed modules
    clean                   Used for ninja clean
    tidy                    A light version of clean, where only the intermediate
                            files are removed
    create_git_ignore       deprecated
    NEED                    The path of NEED file.
                            by default try to open the file in the current path

"""
import os
import sys
import os.path
from collections import namedtuple
import shutil



try:
    from docopt import docopt
    from qp_path import QP_SRC, QP_ROOT, QP_PLUGINS, QP_EZFIO
except ImportError:
    print("source quantum_package.rc")
    raise


def is_module(path_module_rel):
    return os.path.isfile(os.path.join(QP_SRC, path_module_rel, "NEED"))


def is_plugin(path_module_rel):
    return os.path.isfile(os.path.join(QP_PLUGINS, path_module_rel, "NEED"))


def get_binaries(path_module):
    """
    Return the list of binaries
    """
    import subprocess

    try:
        cmd = 'grep -l -i --regexp="^\\s*program\\s" {0}/*.irp.f'.format(path_module)
        process = subprocess.Popen([cmd],
                                   shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        stdout = stdout.decode()
    except OSError:
        return []
    else:
        if not stdout:
            return []
        elif "No such file or directory" not in stdout:
            l_bin = [i.replace(".irp.f", "", 1) for i in stdout.split()]
            return [os.path.realpath(bin_) for bin_ in l_bin]
        else:
            return []



def get_dict_child(l_root_abs=None):
    """Loop over MODULE in  QP_ROOT/src, open all the NEED
    and create a dict[MODULE] = [sub module needed, ...]
    """
    d_ref = dict()

    if not l_root_abs:
        l_root_abs = [QP_SRC]

    for root_abs in l_root_abs:
        for module_rel in os.listdir(root_abs):

            module_abs = os.path.join(root_abs, module_rel)
            try:
                path_file = os.path.join(module_abs, "NEED")

                with open(path_file, "r") as f:
                    l_children = f.read().split()
            except IOError:
                pass
            else:
                if module_rel not in d_ref:
                    d_ref[module_rel] = l_children
               #else:
               #    print "Module {0} alredy defined"
               #    print "Abort"
               #    sys.exit(1)

    return d_ref


def get_l_module_descendant(d_child, l_module):
    """
    From a list of module return the module and descendant
    """

    l = []
    for module in l_module:
        if module not in l:
            l.append(module)
            try:
                l.extend(get_l_module_descendant(d_child, d_child[module]))
            except KeyError:
                print("Error: ", file=sys.stderr)
                print("`{0}` is not a submodule".format(module), file=sys.stderr)
                print("Check the typo (spelling, case, '/', etc.) ", file=sys.stderr)
#                pass
                sys.exit(1)

    return list(set(l))


class ModuleHandler():
    def __init__(self, l_root_abs=None):
        self.dict_child = get_dict_child(l_root_abs)

    @property
    def l_module(self):
        return list(self.dict_child.keys())

    @property
    def dict_parent(self):
        """
        Get a dic of the first parent
        """
        d_child = self.dict_child

        d = {}

        for module_name in d_child:
            d[module_name] = [i for i in list(d_child.keys())
                              if module_name in d_child[i]]

        return d

    @property
    def dict_descendant(self):
        """
        Get a dic of all the genealogy desc (children and all_children)
        """
        d = {}

        d_child = self.dict_child

        for module_name in d_child:
            try:
                d[module_name] = get_l_module_descendant(d_child,
                                                         d_child[module_name])
            except KeyError:
                print("Check NEED for {0}".format(
                    module_name))
                sys.exit(1)

        return d

    @property
    def dict_root(self):
        """
        Return a dict(module_name) = module_boss
        The top node in a tree.
        """
        d_asc = self.dict_parent
        d_desc = self.dict_descendant

        l_all_module = self.l_module

        dict_root = {}

        for module in l_all_module:
            dict_root[module] = [p for p in l_all_module
                                 if module in [p] + d_desc[p] and not d_asc[p]
                                 ][0]

        return dict_root

    def l_descendant_unique(self, l_module):
        d_desc = self.dict_descendant

        d = {}
        for module in l_module:
            for e in d_desc[module]:
                d[e] = 1

        return list(d.keys())

    def l_reduce_tree(self, l_module):
        """For a list of module in input return only the root"""
        l_d_u = self.l_descendant_unique(l_module)
        l_module_reduce = []
        for module in l_module:
            if module not in l_d_u:
                l_module_reduce.append(module)

        return l_module_reduce


if __name__ == '__main__':

    arguments = docopt(__doc__)

    if arguments['--all']:
        l_module = [f for f in os.listdir(QP_SRC)
                    if os.path.isdir(os.path.join(QP_SRC, f))]
        # Remove all produced ezfio_config files
        for filename in os.listdir( os.path.join(QP_EZFIO, "config") ):
            os.remove( os.path.join(QP_EZFIO, "config", filename) )


    elif not arguments['<module_name>']:
        dir_ = os.getcwd()
        l_module = [os.path.basename(dir_)]
    else:
        l_module = arguments['<module_name>']

    for module in l_module:
        if not is_module(module):
            print("{0} is not a valid module. Abort".format(module))
            print("No NEED in it")
            sys.exit(1)

    m = ModuleHandler()

    if arguments['print_descendant']:

        for module in l_module:
            print(" ".join(sorted(m.l_descendant_unique([module]))))

    if arguments["clean"] or arguments["tidy"]:

        l_dir = ['IRPF90_temp', 'IRPF90_man']
        l_file = ["irpf90_entities", "tags", "irpf90.make", "Makefile",
                  "Makefile.depend", ".ninja_log", ".ninja_deps",
                  "ezfio_interface.irp.f"]

        for module in l_module:
            module_abs = os.path.realpath(os.path.join(QP_SRC, module))
            l_symlink = m.l_descendant_unique([module])
            l_exe = get_binaries(module_abs)

            for f in l_dir:
                try:
                    shutil.rmtree(os.path.join(module_abs, f))
                except:
                    pass

            for symlink in l_symlink:
                try:
                    os.unlink(os.path.join(module_abs, symlink))
                except:
                    pass

            for f in l_file:
                try:
                    os.remove(os.path.join(module_abs, f))
                except:
                    pass

            if arguments["clean"]:
                for f in l_exe:

                    try:
                        os.remove(os.path.join(module_abs, f))
                    except:
                        pass


