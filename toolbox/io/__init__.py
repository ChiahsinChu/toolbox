# SPDX-License-Identifier: LGPL-3.0-or-later
"""Input/output operations for various file formats.

This module provides tools for reading and writing:
- CP2K input and output files
- LAMMPS data and dump files
- Template files for various calculators
- Molecular structure files

Examples
--------
>>> from toolbox.io import lammps, cp2k
>>> atoms = lammps.read_dump("trajectory.lammpstrj")
"""
# from . import cp2k, lammps, template

__all__ = ["cp2k", "lammps", "template"]
