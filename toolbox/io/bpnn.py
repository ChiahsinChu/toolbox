# SPDX-License-Identifier: LGPL-3.0-or-later
"""BPNN data reader module.

This module provides functionality to read BPNN (Behler-Parrinello
Neural Network) data files for machine learning potentials.
"""

import numpy as np

from toolbox.utils.unit import AU_TO_ANG, AU_TO_EV


def read_data(fname="input.data"):
    """Read BPNN data file.
    
    Parameters
    ----------
    fname : str, optional
        Path to BPNN data file, by default "input.data"
        
    Returns
    -------
    tuple
        Tuple of (box, coord, charge, symbol, energy, force) where:
        - box: array of box parameters with shape (n_frames, 9)
        - coord: array of coordinates with shape (n_frames, n_atoms*3)
        - charge: array of charges with shape (n_frames, n_atoms)
        - symbol: list of element symbols for each frame
        - energy: array of energies with shape (n_frames,)
        - force: array of forces with shape (n_frames, n_atoms*3)
    """
    box = []
    coord = []
    charge = []
    symbol = []
    energy = []
    force = []

    flag = False
    count = 0
    _symbol = []
    with open(fname, encoding="UTF-8") as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip()
            if line == "begin":
                flag = True
                count = count + 1
            if line == "end":
                flag = False
                symbol.append(_symbol)
                _symbol = []
            if flag is False:
                continue
            line = line.split()
            if line[0] == "lattice":
                box.append(line[1:])
            if line[0] == "atom":
                coord.append(line[1:4])
                charge.append(line[5])
                force.append(line[7:10])
                _symbol.append(line[4])
            if line[0] == "energy":
                energy.append(line[1])

    box = np.array(box, dtype=np.float64)
    if len(box) > 0:
        box = np.reshape(box, (count, 9)) * AU_TO_ANG
    charge = np.array(charge, dtype=np.float64)
    charge = np.reshape(charge, (count, -1))
    coord = np.array(coord, dtype=np.float64)
    coord = np.reshape(coord, (count, -1)) * AU_TO_ANG
    energy = np.array(energy, dtype=np.float64) * AU_TO_EV
    force = np.array(force, dtype=np.float64)
    force = np.reshape(force, (count, -1)) * AU_TO_EV / AU_TO_ANG

    return box, coord, charge, symbol, energy, force
