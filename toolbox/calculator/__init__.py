# SPDX-License-Identifier: LGPL-3.0-or-later
"""Calculator interfaces for computational chemistry and materials science.

This module provides calculator interfaces for various simulation packages:
- CP2KCalculator: Interface for CP2K quantum chemistry package
- LammpsCalculator: Interface for LAMMPS molecular dynamics package
- DPDispatcher: Job dispatcher for Deep Potential calculations
- CP2KDPDispatcher: Specialized dispatcher for CP2K+DP workflows

Examples
--------
>>> from toolbox.calculator import Cp2kCalculator, LammpsCalculator
>>> cp2k_calc = Cp2kCalculator()
>>> lammps_calc = LammpsCalculator()
"""
from contextlib import suppress

from .cp2k import Cp2kCalculator
from .lammps import LammpsCalculator

with suppress(ImportError):
    from .dp import CP2KDPDispatcher, DPDispatcher

__all__ = ["Cp2kCalculator", "LammpsCalculator", "DPDispatcher", "CP2KDPDispatcher"]


def get_template_json(name):
    """Save calculator template as JSON file.
    
    Parameters
    ----------
    name : str
        Name of the template to save
    """
    from ..utils.utils import save_dict
    from . import template

    d = getattr(template, name)
    save_dict(d, name + ".json")
