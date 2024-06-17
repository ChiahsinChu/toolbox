# SPDX-License-Identifier: LGPL-3.0-or-later
from .cp2k import Cp2kCalculator
from .lammps import LammpsCalculator

try:
    from .dp import CP2KDPDispatcher, DPDispatcher
except:
    pass

__all__ = ["Cp2kCalculator", "LammpsCalculator", "DPDispatcher", "CP2KDPDispatcher"]


def get_template_json(name):
    from ..utils.utils import save_dict
    from . import template

    d = getattr(template, name)
    save_dict(d, name + ".json")
