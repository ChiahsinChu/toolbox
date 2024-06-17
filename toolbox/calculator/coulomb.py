# SPDX-License-Identifier: LGPL-3.0-or-later
import numpy as np
from MDAnalysis.lib.distances import distance_array, minimize_vectors
from scipy import constants, special

EPSILON = constants.epsilon_0 / constants.elementary_charge * constants.angstrom
qqrd2e = 1 / (4 * np.pi * EPSILON)
EWALD_F = 1.12837917
EWALD_P = 0.3275911
A1 = 0.254829592
A2 = -0.284496736
A3 = 1.421413741
A4 = -1.453152027
A5 = 1.061405429


def coul(qi, qj, rij):
    """
    qi/qj: e
    rij: A
    return eV
    """
    e = qqrd2e * qi * qj / rij
    return e


def real_coul(qi, qj, xi, xj, gewald, force=False):
    """
    qi/qj: e
    rij: A
    return eV
    """
    rij = np.linalg.norm(xi - xj)
    prefactor = qqrd2e * qi * qj / rij
    if force:
        grij = gewald * rij
        expm2 = np.exp(-grij * grij)
        t = 1.0 / (1.0 + EWALD_P * grij)
        erfc = t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * expm2
        forcecoul = prefactor * (erfc + EWALD_F * grij * expm2)
        fpair = forcecoul / rij**2
        f = (xi - xj) * fpair
        e = prefactor * erfc
        return e, f
    else:
        e = prefactor * special.erfc(rij * gewald)
        return e


class CoulCutCalculator:
    def __init__(self, cutoff=5.0) -> None:
        self.cutoff = cutoff

    def calculate(self, atoms):
        cellpar = atoms.cell.cellpar()
        coords = atoms.get_positions()
        charges = atoms.get_initial_charges()
        dist_mat = distance_array(coords, coords, box=cellpar)
        nat = len(atoms)
        forces = np.zeros((nat, 3))
        energy = 0.0
        for ii in range(nat):
            force_mask = (dist_mat[ii] < self.cutoff) & (np.arange(nat) > ii)
            sel_ids = np.where(force_mask)[0]
            forcecoul = qqrd2e * charges[ii] * charges[sel_ids] / dist_mat[ii][sel_ids]
            fpair = forcecoul / dist_mat[ii][sel_ids] ** 2
            delta_x = minimize_vectors(
                coords[ii].reshape(1, 3) - coords[sel_ids], box=cellpar
            )
            forces[ii] += np.sum(delta_x * fpair.reshape(-1, 1), axis=0)
            forces[sel_ids] -= delta_x * fpair.reshape(-1, 1)
            energy += np.sum(forcecoul)
        return energy, forces


class CoulLongCalculator:
    def __init__(self, gewald, cutoff=5.0) -> None:
        self.cutoff = cutoff
        self.gewald = gewald

    def calculate(self, atoms):
        cellpar = atoms.cell.cellpar()
        coords = atoms.get_positions()
        charges = atoms.get_initial_charges()
        dist_mat = distance_array(coords, coords, box=cellpar)
        nat = len(atoms)
        forces = np.zeros((nat, 3))
        energy = 0.0
        for ii in range(nat):
            force_mask = (dist_mat[ii] < self.cutoff) & (np.arange(nat) > ii)
            sel_ids = np.where(force_mask)[0]
            prefactor = qqrd2e * charges[ii] * charges[sel_ids] / dist_mat[ii][sel_ids]
            grij = self.gewald * dist_mat[ii][sel_ids]
            expm2 = np.exp(-grij * grij)
            t = 1.0 / (1.0 + EWALD_P * grij)
            erfc = t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * expm2
            forcecoul = prefactor * (erfc + EWALD_F * grij * expm2)
            fpair = forcecoul / dist_mat[ii][sel_ids] ** 2
            delta_x = minimize_vectors(
                coords[ii].reshape(1, 3) - coords[sel_ids], box=cellpar
            )
            forces[ii] += np.sum(delta_x * fpair.reshape(-1, 1), axis=0)
            forces[sel_ids] -= delta_x * fpair.reshape(-1, 1)
            energy += np.sum(prefactor * erfc)
        return energy, forces
