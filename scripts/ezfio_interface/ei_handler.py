#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Welcome to the ei_handler.
We will create all the ezfio related stuff from a EZFIO.cfg file.

Usage:
    ei_handler.py [--path_module=<module>]
                  [--irpf90]
                  [--ezfio_config]
                  [--ocaml]
                  [--ezfio_default]
    ei_handler.py list_supported_types
    ei_handler.py ocaml_global

By default all the option are executed.

Options:
    -h --help
    --irpf90        Create the `<module>/ezfio_interface.irpf90`
                         which contains all the providers needed
    --ezfio_config  Create the `<module_lower>_ezfio_interface_config` in
                          `${QP_EZFIO}/config/`
    --ocaml         Create all the stuff needed by *qp_edit*:
                             -`Input_<module_lower>.ml` and
                             - <module_lower>_ezfio_interface_default`
    ocaml_global    Create the *qp_edit*

Format specification :

Required:
    [<provider_name>]   The name of the provider in irp.f90 and in the EZFIO lib
    doc:<str>           The plain text documentation
    type:<str>          A Fancy_type supported by the ocaml.
                            type `ei_handler.py get_supported_type` for a list
    interface:<str>     The interface is list of string sepeared by ","  which can contain :
                          - ezfio (if you only whant the ezfiolib)
                          - provider (if you want the provider)
                          - ocaml (if you want the ocaml gestion)
Optional:
    default: <str>      The default value needed,
                            if 'ocaml' is in interface list.
                           ! No list is allowed for now !
    size: <str>         The size information.
                            (by default is one)
                            Example : 1, =sum(ao_num); (ao_num,3)
                            ATTENTION : The module and the value are separed by a "." not a "_".
                            For example (determinants.n_det)
    ezfio_name: <str>   The name for the EZFIO lib
                             (by default is <provider_name>)
    ezfio_dir: <str>    Will be the folder of EZFIO.
                              (by default is <module_lower>)

Example of EZFIO.cfg:
```
[thresh_SCF]
doc: Threshold on the convergence of the Hartree Fock energy
type: Threshold
default: 1.e-10
interface: provider,ezfio,ocaml
size: 1

[energy]
type: double precision
doc: Calculated HF energy
interface: ezfio
```
"""
from docopt import docopt

import sys
import os
import os.path

import configparser

from collections import defaultdict
from collections import namedtuple

from qp_decorator import cache

from os import listdir
from os.path import isdir, join, exists


from qp_path import QP_ROOT, QP_SRC, QP_OCAML, QP_EZFIO

Type = namedtuple('Type', 'fancy ocaml fortran')
Module = namedtuple('Module', 'path lower')


def is_bool(str_):
    """
    Take a string, if is a bool return the conversion into
        fortran and ocaml.
    """
    if "true" in str_.strip().lower():
        return Type(None, "true", ".True.")
    elif "false" in str_.strip().lower():
        return Type(None, "false", ".False")
    else:
        raise TypeError


@cache
def get_type_dict():
    """
    This function makes the correspondance between the type of value read in
    EZFIO.cfg into the f90 and OCaml type.
    return fancy_type[fancy_type] = namedtuple('Type', 'ocaml fortran')
    For example fancy_type['Ndet'].fortran = integer
                                  .ocaml   = int
    """

    # ~#~#~#~ #
    # I n i t #
    # ~#~#~#~ #

    fancy_type = defaultdict(dict)

    # ~#~#~#~#~#~#~#~ #
    # R a w _ t y p e #
    # ~#~#~#~#~#~#~#~ #

    fancy_type['integer'] = Type(None, "int", "integer")
    fancy_type['integer*8'] = Type(None, "int", "integer*8")

    fancy_type['int'] = Type(None, "int", "integer")
    fancy_type['int64'] = Type(None, "int64", "integer*8")

    fancy_type['float'] = Type(None, "float", "double precision")
    fancy_type['double precision'] = Type(None, "float", "double precision")

    fancy_type['logical'] = Type(None, "bool", "logical")
    fancy_type['bool'] = Type(None, "bool", "logical")

    fancy_type['character*(32)'] = Type(None, "string", "character*(32)")
    fancy_type['character*(64)'] = Type(None, "string", "character*(64)")
    fancy_type['character*(256)'] = Type(None, "string", "character*(256)")

    # ~#~#~#~#~#~#~#~ #
    # q p _ t y p e s #
    # ~#~#~#~#~#~#~#~ #

    # Dict to change ocaml LowLevel type into FortranLowLevel type
    ocaml_to_fortran = {"int": "integer",
                        "int64": "integer*8",
                        "float": "double precision",
                        "logical": "logical",
                        "string": "character*32"}

    # Read and parse qptype generate
    src = join(QP_OCAML, "qptypes_generator.ml")

    with open(src, "r") as f:
        r = f.read()

        # Generate
        l_gen = [i for i in r.splitlines() if i.strip().startswith("*")]

        # Untouch
        b = r.find('let untouched = "')
        e = r.find('let parse_input', b)

        l_un = [i for i in r[b:e].splitlines() if i.strip().startswith("module")]

    # ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~ #
    # q p _ t y p e s _ g e n e r a t e #
    # ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~ #

    # Read the fancy_type, the ocaml. and convert the ocaml to the fortran
    for i in l_gen + l_un:
        str_fancy_type = i.split()[1].strip()
        str_ocaml_type = i.split()[3]

        if str_ocaml_type != 'sig':
            str_fortran_type = ocaml_to_fortran[str_ocaml_type]
        else:
            str_fortran_type = 'character*(32)'
            str_ocaml_type = 'string'

        fancy_type[str_fancy_type] = Type(str_fancy_type,
                                          str_ocaml_type,
                                          str_fortran_type)

    fancy_type["MO_class"] = Type("MO_class", "MO_class", "character*(32)")

    # ~#~#~#~#~#~#~#~ #
    # F i n a l i z e #
    # ~#~#~#~#~#~#~#~ #
    return dict(fancy_type)


type_dict = get_type_dict()


def get_dict_config_file(module_obj):
    """
    Input:
        module_obj.path is the config file 
            (for example FULL_PATH/EZFIO.cfg)
        module_obj.lower is the MODULE name lowered
            (Ex fullci)

    Return a dict d[provider_name] = {type,
                                      doc,
                                      ezfio_name,
                                      ezfio_dir,
                                      size,
                                      interface,
                                      default}
    """
    # ~#~#~#~ #
    # I n i t #
    # ~#~#~#~ #
    d = defaultdict(dict)
    l_info_optional = ["ezfio_dir", "ezfio_name", "size"]

    # ~#~#~#~#~#~#~#~#~#~#~ #
    # L o a d _ C o n f i g #
    # ~#~#~#~#~#~#~#~#~#~#~ #

    config_file = configparser.ConfigParser()
    config_file.read_file(open(module_obj.path))

    # ~#~#~#~#~#~#~#~#~ #
    # F i l l _ d i c t #
    # ~#~#~#~#~#~#~#~#~ #

    def error(o, p, c):
        "o option ; p provider_name ;c module_obj.path"
        print("You need a {0} for {1} in {2}".format(o, p, c))

    for section in config_file.sections():
        # pvd = provider
        pvd = section.lower()

        d[pvd]["module"] = module_obj

        # Create the dictionary which contains the default value 
        d_default = {"ezfio_name": pvd,
                     "ezfio_dir": module_obj.lower,
                     "size": "1"}

        # Check if type is avalaible
        try:
            type_ = config_file.get(section, "type").strip()
        except configparser.NoOptionError:
            error("type", pvd, module_obj.path)
            sys.exit(1)

        if type_ not in type_dict:
            print("{0} not avalaible. Choose in:".format(type_).strip())
            print(", ".join(sorted([i for i in type_dict])))
            sys.exit(1)
        else:
            d[pvd]["type"] = type_dict[type_]

        # Fill the dict with REQUIRED information
        try:
            d[pvd]["doc"] = config_file.get(section, "doc")
        except configparser.NoOptionError:
            error("doc", pvd, module_obj.path)
            sys.exit(1)

        try:
            interface = [i.lower().strip() for i in config_file.get(section, "interface").split(",")]
        except configparser.NoOptionError:
            error("doc", pvd, module_obj.path)
            sys.exit(1)
        else:
            if not any(i in ["ezfio", "provider", "ocaml"] for i in interface):
                print("Bad keyword for interface for {0}".format(pvd))
                sys.exit(1)
            else:
                d[pvd]["interface"] = interface

        # Fill the dict with OPTIONAL information
        for option in l_info_optional:
            try:
                d[pvd][option] = config_file.get(section, option).lower()
            except configparser.NoOptionError:
                if option in d_default:
                    d[pvd][option] = d_default[option]

        # If interface is input we need a default value information

        try:
            default_raw = config_file.get(section, "default")
        except configparser.NoOptionError:
            if "ocaml" in d[pvd]["interface"]:
                error("default", pvd, module_obj.path)
                sys.exit(1)
            else:
                pass
        else:
            try:
                d[pvd]["default"] = is_bool(default_raw)
            except TypeError:
                d[pvd]["default"] = Type(None, default_raw, default_raw)

    return dict(d)


def create_ezfio_provider(dict_ezfio_cfg):
    import re

    """
    From dict d[provider_name] = {type,
                                  doc,
                                  ezfio_name,
                                  ezfio_dir,
                                  interface,
                                  default
                                  size}
    create the a list which contains all the code for the provider
    output = output_dict_info['ezfio_dir'
    return [code, ...]
    """

    from ezfio_generate_provider import EZFIO_Provider, gen_ezfio_provider_disk_access
    dict_code_provider = dict()

    ez_p = EZFIO_Provider()
    for provider_name, dict_info in dict_ezfio_cfg.items():
        if "provider" in dict_info["interface"]:
            ez_p.set_type(dict_info['type'].fortran)
            ez_p.set_name(provider_name)
            ez_p.set_doc(dict_info['doc'])
            ez_p.set_ezfio_dir(dict_info['ezfio_dir'])
            ez_p.set_ezfio_name(dict_info['ezfio_name'])
            ez_p.set_output("6")

            # (nuclei.nucl_num,pseudo.klocmax) => (nucl_num,klocmax)
            ez_p.set_size(re.sub(r'\w+\.', "", dict_info['size']))

            str_ = str(ez_p) + "\n"
            if dict_info['type'].fancy == 'Disk_access':

                allowed_prefix = ['disk_access', 'io']
                assert (any(provider_name.startswith(p) for p in allowed_prefix))

                provider_name_c = provider_name
                for p in allowed_prefix:
                    if provider_name_c.startswith(p):
                        provider_name_c = provider_name_c.replace(p+'_','',1)

                str_ += gen_ezfio_provider_disk_access(provider_name, provider_name_c)

            dict_code_provider[provider_name] = str_

    return dict_code_provider


def save_ezfio_provider(path_head, dict_code_provider):
    """
    Write in path_head/"ezfio_interface.irp.f" the value of dict_code_provider
    """

    path = "{0}/ezfio_interface.irp.f".format(path_head)

    l_output = ["! DO NOT MODIFY BY HAND",
                "! Created by $QP_ROOT/scripts/ezfio_interface/ei_handler.py",
                "! from file {0}/EZFIO.cfg".format(path_head),
                "\n"]

    l_output += [code for code in list(dict_code_provider.values())]

    output = "\n".join(l_output)

    with open(path, "w+") as f:
        f.write(output)


def create_ezfio_stuff(dict_ezfio_cfg, config_or_default="config"):
    """
    From dict_ezfio_cfg[provider_name] = {type, default, ezfio_name,ezfio_dir,doc}
    Return the string ezfio_interface_config
    """

    def size_format_to_ezfio(size_raw):
        """
        If size_raw == "=" is a formula -> do nothing; return
        Else convert the range of a multidimential array
           (12,begin:end) into (12,begin+end+1) for example
        If the values are between parenthses ->  do nothing; return
        """

        size_raw = str(size_raw)
        if size_raw.startswith('='):
            size_convert = size_raw.replace('.', '_')
        else:
            size_raw = provider_info["size"].translate(str.maketrans("","","()"))
            size_raw = size_raw.replace('.', '_')

            a_size_raw = []
            for dim in size_raw.split(","):
                try:
                    (begin, end) = list(map(str.strip, dim.split(":")))
                except ValueError:
                    a_size_raw.append(dim.strip())
                else:
                    if begin[0] == '-':
                        a_size_raw.append("{0}+{1}+1".format(end, begin[1:]))
                    else:
                        a_size_raw.append("{0}-{1}+1".format(end, begin))

            size_raw = ",".join(a_size_raw)

            size_convert = "({0})".format(size_raw)
        return size_convert

    def create_format_string(size):
        """
        Take a size number and
            return the string format for being right align with this offset
        """
        return "{{0:<{0}}}".format(size).format

    # ~#~#~#~#~#~#~#~#~#~#~# #
    #  F o r m a t _ i n f o #
    # ~#~#~#~#~#~#~#~#~#~#~# #

    lenmax_name = max([len(i) for i in dict_ezfio_cfg])
    lenmax_type = max([len(i["type"].fortran)
                       for i in list(dict_ezfio_cfg.values())])

    str_name_format = create_format_string(lenmax_name + 2)
    str_type_format = create_format_string(lenmax_type + 2)

    # ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~# #
    #  C r e a t e _ t h e _ s t r i n g #
    # ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~# #

    # Checking is many ezfio_dir provided
    l_ezfio_dir = [d['ezfio_dir'] for d in list(dict_ezfio_cfg.values())]

    if not l_ezfio_dir.count(l_ezfio_dir[0]) == len(l_ezfio_dir):
        print("You have many ezfio_dir. Not supported yet", file=sys.stderr)
        raise TypeError
    else:
        result = [l_ezfio_dir[0]]

    for provider_name, provider_info in sorted(dict_ezfio_cfg.items()):

        # Get the value from dict
        name_raw = provider_info["ezfio_name"].lower()

        fortran_type_raw = provider_info["type"].fortran

        if "size" in provider_info and not provider_info["size"] == "1":
            size_raw = provider_info["size"]
        else:
            size_raw = None

        # It is the last so we don't need to right align it
        str_size = size_format_to_ezfio(size_raw) if size_raw else ""

        if "default" in provider_info and provider_info["default"].fortran.startswith("="):
            str_default = provider_info["default"].fortran.replace('.', '_')
        else:
            str_default = ""

        # Get the string in to good format (left align and co)
        str_name = str_name_format(name_raw)
        str_fortran_type = str_type_format(fortran_type_raw)

        # Return the string
        if config_or_default == "config":
            s = "  {0} {1} {2} {3}".format(str_name, str_fortran_type, str_size, str_default)
        elif config_or_default == "default":
            try:
                str_value = provider_info["default"].ocaml
            except KeyError:
                continue
            else:
                s = "  {0} {1}".format(str_name, str_value)
        else:
            raise KeyError
        # Append
        result.append(s)
    result.append("\n")

    return "\n".join(result)


def create_ezfio_config(dict_ezfio_cfg):
    return create_ezfio_stuff(dict_ezfio_cfg,
                              config_or_default="config")


def save_ezfio_config(module_lower, str_ezfio_config):
    """
    Write the str_ezfio_config in
    "$QP_ROOT/EZFIO/{0}.ezfio_interface_config".format(module_lower)
    """
    name = "{0}.ezfio_interface_config".format(module_lower)
    path = os.path.join(QP_EZFIO, "config", name)

    with open(path, "w+") as f:
        f.write(str_ezfio_config)


def create_ezfio_default(dict_ezfio_cfg):
    return create_ezfio_stuff(dict_ezfio_cfg,
                              config_or_default="default")


def save_ezfio_default(module_lower, str_ezfio_default):
    """
    Write the str_ezfio_config in
    "$QP_ROOT/data/ezfio_defaults/{0}.ezfio_interface_default".format(module_lower)
    """

    root_ezfio_default = "{0}/data/ezfio_defaults/".format(
        QP_ROOT)
    path = "{0}/{1}.ezfio_interface_default".format(root_ezfio_default,
                                                    module_lower)
    with open(path, "w+") as f:
        f.write(str_ezfio_default)


def create_ocaml_input(dict_ezfio_cfg, module_lower):

    # ~#~#~#~# #
    #  I n i t #
    # ~#~#~#~# #

    from ezfio_generate_ocaml import EZFIO_ocaml

    l_ezfio_name = []
    l_type = []
    l_doc = []

    for k, v in dict_ezfio_cfg.items():
        if "ocaml" in v['interface']:
            l_ezfio_name.append(v['ezfio_name'])
            l_type.append(v["type"])
            l_doc.append(v["doc"])

    if not l_ezfio_name:
        raise ValueError

    e_glob = EZFIO_ocaml(l_ezfio_name=l_ezfio_name,
                         l_type=l_type,
                         l_doc=l_doc)

    # ~#~#~#~#~#~#~#~# #
    #  C r e a t i o n #
    # ~#~#~#~#~#~#~#~# #

    template = ['(* =~=~ *)',
                '(* Init *)',
                '(* =~=~ *)',
                ""]

    template += ["open Qptypes;;",
                 "open Qputils;;",
                 "open Sexplib.Std;;",
                 "",
                 "module {0} : sig".format(module_lower.capitalize())]

    template += [e_glob.create_type()]

    template += ["  val read  : unit -> t option",
                 "  val write : t-> unit",
                 "  val to_string : t -> string",
                 "  val to_rst : t -> Rst_string.t",
                 "  val of_rst : Rst_string.t -> t option",
                 "end = struct"]

    template += [e_glob.create_type()]

    template += ['',
                 '  let get_default = Qpackage.get_ezfio_default "{0}";;'.format(module_lower),
                 '']

    template += ['(* =~=~=~=~=~=~==~=~=~=~=~=~ *)',
                 '(* Generate Special Function *)',
                 '(* =~=~=~==~=~~=~=~=~=~=~=~=~ *)',
                 ""]

    for provider_name, d_val in sorted(dict_ezfio_cfg.items()):

        if 'default' not in d_val:
            continue

        ezfio_dir = d_val["ezfio_dir"]
        ezfio_name = d_val["ezfio_name"]

        e = EZFIO_ocaml(ezfio_dir=ezfio_dir,
                        ezfio_name=ezfio_name,
                        type=d_val["type"])

        template += [e.create_read(),
                     e.create_write(),
                     ""]

    template += ['(* =~=~=~=~=~=~=~=~=~=~=~=~ *)',
                 '(* Generate Global Function *)',
                 '(* =~=~=~=~=~=~=~=~=~=~=~=~ *)',
                 ""]

    template += [e_glob.create_read_global(),
                 e_glob.create_write_global(),
                 e_glob.create_to_string(),
                 e_glob.create_to_rst()]

    template += ["  include Generic_input_of_rst;;",
                 "  let of_rst = of_rst t_of_sexp;;",
                 "",
                 "end"]

    result = "\n".join(template)
    result = result.replace("String.of_string","string_of_string")
    result = result.replace("String.to_string","string_of_string")
    result = result.replace("Int.of_string","int_of_string")
    result = result.replace("Int.to_string","string_of_int")
    result = result.replace("Float.of_string","float_of_string")
    result = result.replace("Float.to_string","string_of_float")
    result = result.replace("Bool.of_string","bool_of_string")
    result = result.replace("Bool.to_string","string_of_bool")
    return result


def save_ocaml_input(module_lower, str_ocaml_input):
    """
    Write the str_ocaml_input in
    qp_path.QP_OCAML/Input_{0}.ml".format(module_lower)
    """

    name = "Input_{0}.ml".format(module_lower)

    path = join(QP_OCAML, name)

    with open(path, "w+") as f:
        f.write(str_ocaml_input)


def get_l_module_with_auto_generate_ocaml_lower():
    """
    Get all modules which have EZFIO.cfg with OCaml data
        (NB `search` in all the lines and `match` only in one)
    """

    # ~#~#~#~#~#~#~#~ #
    # L _ f o l d e r #
    # ~#~#~#~#~#~#~#~ #

    l_folder = [f for f in listdir(QP_SRC) if isdir(join(QP_SRC, f))]

    # ~#~#~#~#~#~#~#~#~#~#~#~#~#~ #
    # L _ m o d u l e _ l o w e r #
    # ~#~#~#~#~#~#~#~#~#~#~#~#~#~ #

    l_module_lower = []
    import re
    p = re.compile(r'interface:.*ocaml')

    for f in l_folder:
        path = "{0}/{1}/EZFIO.cfg".format(QP_SRC, f)
        if exists(path):
            with open(path, 'r') as file_:
                if p.search(file_.read()):
                    l_module_lower.append(f.lower())

    # ~#~#~#~#~#~ #
    # R e t u r n #
    # ~#~#~#~#~#~ #

    return l_module_lower


def create_ocaml_input_global(l_module_with_auto_generate_ocaml_lower):
    """
    Create the Input_auto_generated.ml and qp_edit.ml str
    """

    # ~#~#~#~#~#~#~#~# #
    #  C r e a t i o n #
    # ~#~#~#~#~#~#~#~# #

    from ezfio_generate_ocaml import EZFIO_ocaml

    path = QP_ROOT + "/scripts/ezfio_interface/qp_edit_template"

    with open(path, "r") as f:
        template_raw = f.read()

    e = EZFIO_ocaml(l_module_lower=l_module_with_auto_generate_ocaml_lower)

    template = template_raw.format(keywords=e.create_qp_keywords(),
                                   keywords_to_string=e.create_qp_keywords_to_string(),
                                   section_to_rst=e.create_qp_section_to_rst(),
                                   write=e.create_qp_write(),
                                   tasks=e.create_qp_tasks())

    input_auto = e.create_input_auto_generated()

    return (template, input_auto)


def save_ocaml_input_auto(str_ocaml_input_global):
    """
    Write the str_ocaml_input in
    qp_path.QP_OCAML/Input_auto_generated.ml
    """

    name = "Input_auto_generated.ml"
    path = join(QP_OCAML, name)

    with open(path, "w+") as f:
        f.write(str_ocaml_input_global)


def save_ocaml_qp_edit(str_ocaml_qp_edit):
    """
    Write the str_ocaml_qp_edit in
    qp_path.QP_OCAML/qp_edit.ml
    """

    name = "qp_edit.ml"
    path = join(QP_OCAML, name)

    with open(path, "w+") as f:
        f.write(str_ocaml_qp_edit)


def code_generation(arguments, dict_ezfio_cfg, m):

    module_lower = m.lower
    path_dirname = m.path.replace("/EZFIO.cfg", "")

    # ~#~#~#~#~#~#~#~#~#~ #
    # W h a t _ t o _ d o #
    # ~#~#~#~#~#~#~#~#~#~ #
    if any([arguments[i] for i in ["--irpf90",
                                   "--ezfio_config",
                                   "--ocaml",
                                   "--ezfio_default"]]):
        # User changer somme argument, do what he want
        do_all = False
    else:
        # Do all the stuff
        do_all = True

    # ~#~#~#~#~#~#~ #
    # I R P . f 9 0 #
    # ~#~#~#~#~#~#~ #

    if do_all or arguments["--irpf90"]:
        l_str_code = create_ezfio_provider(dict_ezfio_cfg)
        save_ezfio_provider(path_dirname, l_str_code)

    # ~#~#~#~#~#~#~#~#~#~#~#~ #
    # e z f i o _ c o n f i g #
    # ~#~#~#~#~#~#~#~#~#~#~#~ #

    if do_all or arguments["--ezfio_config"]:
        str_ezfio_config = create_ezfio_config(dict_ezfio_cfg)
        save_ezfio_config(module_lower, str_ezfio_config)

    # ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~ #
    # O c a m l & e z f i o _ d e f a u l t #
    # ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~ #
    if do_all or arguments["--ocaml"]:
        try:
            str_ocaml_input = create_ocaml_input(dict_ezfio_cfg, module_lower)
        except ValueError:
            pass
        else:
            save_ocaml_input(module_lower, str_ocaml_input)

        str_ezfio_default = create_ezfio_default(dict_ezfio_cfg)
        save_ezfio_default(module_lower, str_ezfio_default)

if __name__ == "__main__":
    arguments = docopt(__doc__)
    # ___
    #  |  ._  o _|_
    # _|_ | | |  |_
    #
    if arguments["list_supported_types"]:
        for i in sorted(get_type_dict()):
            print(i)
        sys.exit(0)

    if arguments["ocaml_global"]:

        # ~#~#~#~# #
        #  I n i t #
        # ~#~#~#~# #

        l_module = get_l_module_with_auto_generate_ocaml_lower()

        str_ocaml_qp_edit, str_ocaml_input_auto = create_ocaml_input_global(l_module)
        save_ocaml_input_auto(str_ocaml_input_auto)
        save_ocaml_qp_edit(str_ocaml_qp_edit)
        sys.exit(0)

    # ~#~#~#~#~#~#~#~#~#~#~#~#~#~# #
    #  G e t _ m o d u l e _ d i r #
    # ~#~#~#~#~#~#~#~#~#~#~#~#~#~# #

    if arguments["--path_module"]:
        path_dirname = os.path.abspath(arguments["--path_module"])
    else:
        path_dirname = os.getcwd()

    root_module = os.path.split(path_dirname)[1]

    l_module = [root_module]
    # ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~# #
    #  G e t _ l _ d i c t _ e z f i o _ c f g #
    # ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~# #

    l_module_with_ezfio = []

    for f in l_module:
        path = join(QP_SRC, f, "EZFIO.cfg")
        if exists(path):
            l_module_with_ezfio.append(Module(path, f.lower()))

    l_dict_ezfio_cfg = [(m, get_dict_config_file(m)) for m in l_module_with_ezfio]

    #  _
    # /   _   _|  _     _   _  ._   _  ._ _. _|_ o  _  ._
    # \_ (_) (_| (/_   (_| (/_ | | (/_ | (_|  |_ | (_) | |
    #                   _|

    for (m, dict_ezfio_cfg) in l_dict_ezfio_cfg:
        code_generation(arguments, dict_ezfio_cfg, m)
