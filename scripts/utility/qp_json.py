#!/usr/bin/env python3

"""Strict reader of the QP3 JSON result documents.

Every run writes one complete, valid JSON document to
``<EZFIO>/json/NNNNN.json`` (see ``src/json/README.rst``).  Nothing here
repairs files: a file that does not parse is skipped, never rewritten.
"""

import json
import os
import re

_NAME = re.compile(r'^\d{5}\.json$')


def files(ezfio_filename):
    """Names of the JSON documents of an EZFIO, oldest first."""
    prefix = os.path.join(ezfio_filename, 'json')
    if not os.path.isdir(prefix):
        return []
    return sorted(n for n in os.listdir(prefix) if _NAME.match(n))


def load(filename):
    """Parses one document."""
    with open(filename) as f:
        return json.load(f)


def load_all(ezfio_filename):
    """All valid documents of an EZFIO, keyed by file name."""
    prefix = os.path.join(ezfio_filename, 'json')
    out = {}
    for n in files(ezfio_filename):
        try:
            d = load(os.path.join(prefix, n))
        except (OSError, ValueError):
            continue
        if d.get('schema') == 'qp3/1':
            out[n] = d
    return out


def load_last(ezfio_filename):
    """The newest document of an EZFIO, the state of the calculation."""
    names = files(ezfio_filename)
    if not names:
        raise FileNotFoundError(
            os.path.join(ezfio_filename, 'json', 'NNNNN.json'))
    return load(os.path.join(ezfio_filename, 'json', names[-1]))
