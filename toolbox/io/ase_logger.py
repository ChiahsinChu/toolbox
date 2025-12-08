# SPDX-License-Identifier: LGPL-3.0-or-later
"""ASE MD Logger module.

This module provides an extended MDLogger class that can write
trajectory files during molecular dynamics simulations.
"""

from typing import IO, Any

from ase import Atoms
from ase.md import MDLogger as _MDLogger


class MDLogger(_MDLogger):
    """Extended MD Logger for ASE molecular dynamics.
    
    This class extends ASE's MDLogger to provide additional functionality
    for writing trajectory files during MD simulations.
    """
    
    def __init__(
        self,
        dyn: Any,  # not fully annotated so far to avoid a circular import
        atoms: Atoms,
        logfile: IO | str,
        fname: str,
        write_kwargs: dict = None,
        **kwargs,
    ):
        if write_kwargs is None:
            write_kwargs = {}
        super().__init__(dyn, atoms, logfile, **kwargs)
        self.fname = fname
        self.write_kwargs = write_kwargs

    def __call__(self):
        """Write trajectory data to file.
        
        This method is called during MD simulation to write the current
        atomic configuration to the specified file.
        """
        super().__call__()
        self.atoms.write(self.fname, append=True, **self.write_kwargs)
