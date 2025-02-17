# SPDX-License-Identifier: LGPL-3.0-or-later
import glob
import os
from typing import List

import numpy as np
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator


def set_energy_and_forces(atoms: Atoms, energy: float, forces: np.ndarray):
    """
    Set energy and forces to atoms object.
    Make atoms.get_potential_energy() and atoms.get_forces() return the given energy and forces.

    Parameters
    ----------
    atoms : ase.Atoms
        The atoms object.
    energy : float
        The energy.
    forces : np.ndarray
        The forces.

    """
    calc = SinglePointCalculator(atoms)
    atoms.set_calculator(calc)
    atoms.calc.results["energy"] = energy
    atoms.calc.results["forces"] = forces


def update_atype(dname: str, type_map: List[str]):
    fnames = glob.glob(os.path.join(dname, "**/type.raw"), recursive=True)
    fnames.sort()
    for fname in fnames:
        dname = os.path.dirname(fname)
        _type_map = np.loadtxt(os.path.join(dname, "type_map.raw"), dtype=str)
        _atype = np.loadtxt(fname, dtype=int)
        symbols = _type_map[_atype]
        # generate new type.raw according to the input type_map
        assert np.all(np.isin(symbols, type_map)), "Invalid type_map"
        new_atype = np.array([type_map.index(s) for s in symbols], dtype=int)
        np.savetxt(fname, new_atype, fmt="%d")
        np.savetxt(os.path.join(dname, "type_map.raw"), type_map, fmt="%s")


# if __name__ == "__main__":
#     import glob

#     from ase import io

#     from dpdata import LabeledSystem, MultiSystems
#     from toolbox.utils.unit import AU_TO_EV, AU_TO_EV_EVERY_ANG
#     from toolbox.io.dpdata import set_energy_and_forces

#     fnames_pos = glob.glob("./**/*-pos-*.xyz", recursive=True)
#     fnames_frc = glob.glob("./**/*-frc-*.xyz", recursive=True)
#     assert len(fnames_pos) == len(fnames_frc)
#     fnames_pos.sort()
#     fnames_frc.sort()

#     ms = MultiSystems()

#     for fname_pos, fname_frc in zip(fnames_pos, fnames_frc):
#         traj = io.read(fname_pos, ":")
#         frc_traj = io.read(fname_frc, ":")
#         for atoms, atoms_frc in zip(traj, frc_traj):
#             forces = atoms_frc.get_positions() * AU_TO_EV_EVERY_ANG
#             set_energy_and_forces(atoms, atoms.info["E"] * AU_TO_EV, forces)
#             ms.append(LabeledSystem(atoms, fmt="ase/structure"))

#     ms.to_deepmd_npy("deepmd")
