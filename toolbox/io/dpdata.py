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
