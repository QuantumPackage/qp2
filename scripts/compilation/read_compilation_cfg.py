#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import os, sys
import ConfigParser


def get_l_option_section(config):
    """List of options chosen by the user"""
    l = [o for o in ['OPENMP'] if config.getboolean("OPTION", o)]
    l.append(config.get("OPTION", "MODE").strip())
    return l


def get_compilation_option(pwd_cfg, flag_name):
    """
    Return the flag compilation of a compile.cfg located in pwd_cfg
    """
    if not os.path.isfile(pwd_cfg):
        print "Configuration file %s not found"%(pwd_cfg)
        sys.exit(1)

    config = ConfigParser.ConfigParser()
    config.read(pwd_cfg)

    if flag_name == "FC" and config.getboolean("OPTION","CACHE"):
        l = ["cache_compile.py"]
    else:
        l = []

    l_option_section = get_l_option_section(config)

    for section in ["COMMON"] + l_option_section:
        try:
            l.extend(config.get(section, flag_name).split())
        except ConfigParser.NoOptionError:
            pass

    return " ".join(l)

if __name__ == '__main__':

    qpackage_root = os.environ['QP_ROOT']
    pwd_cfg = os.path.join(qpackage_root, "config/ifort_gpi2.cfg")

    print get_compilation_option(pwd_cfg, "FC")
    print get_compilation_option(pwd_cfg, "FCFLAGS")
