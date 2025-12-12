# SPDX-License-Identifier: LGPL-3.0-or-later
"""D3 dispersion correction calculator.

This module provides a calculator for Grimme's D3 dispersion
corrections using the tad_dftd3 library.
"""

from typing import Optional

import tad_dftd3 as d3
import torch
from ase import Atoms

torch.set_default_dtype(torch.float64)


class D3Calculator:
    """Calculator for D3 dispersion corrections.
    
    This calculator computes Grimme's D3 dispersion corrections
    using the tad_dftd3 library.
    """
    
    def __init__(self, params: Optional[dict] = None) -> None:
        """Initialize D3Calculator.
        
        Parameters
        ----------
        params : Optional[Dict], optional
            D3 parameters dictionary, by default None.
            If None, uses r²SCAN-D3(BJ) parameters.
        """
        self.ref = d3.reference.Reference()
        if params is None:
            # r²SCAN-D3(BJ)
            self.params = dict(
                a1=torch.tensor(0.49484001),
                s8=torch.tensor(0.78981345),
                a2=torch.tensor(5.73083694),
            )
        else:
            self.params = params

    def run(self, atoms: Atoms):
        """Calculate D3 dispersion energy.
        
        Parameters
        ----------
        atoms : ase.Atoms
            ASE Atoms object with atomic numbers and positions
            
        Returns
        -------
        torch.Tensor
            D3 dispersion energy
        """
        # nframes * natoms
        numbers = atoms.get_atomic_numbers().reshape(-1, len(atoms))
        numbers = torch.tensor(numbers, dtype=torch.int64)
        # nframes * 3 * natoms
        positions = atoms.get_positions().reshape(-1, len(atoms), 3)
        positions = torch.tensor(positions)

        cn = d3.ncoord.cn_d3(numbers, positions)
        weights = d3.model.weight_references(numbers, cn, self.ref)
        c6 = d3.model.atomic_c6(numbers, weights, self.ref)
        energy = d3.disp.dispersion(numbers, positions, self.params, c6)
        return energy.sum()
