# SPDX-License-Identifier: LGPL-3.0-or-later
"""Coulomb interaction calculators for molecular simulations.

This module provides calculators for:
- Coulomb cutoff interactions
- Ewald summation methods
- Particle Mesh Ewald (PME) calculations
- D3 dispersion corrections

These calculators are used for computing electrostatic interactions
in molecular dynamics simulations.
"""

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

try:
    import torch
except ImportError:
    import warnings

    warnings.warn("torch not found, DMFFPMECalculator cannot be used.", stacklevel=2)
try:
    from dp_dmff.dmff.constants import ENERGY_COEFF
    from dp_dmff.dmff.pme import energy_pme, setup_ewald_parameters
    from dp_dmff.dmff.recip import Ck_1, generate_pme_recip
    from dp_dmff.dmff.settings import DEVICE
except ImportError:
    import warnings

    warnings.warn("dp_dmff not found, DMFFPMECalculator cannot be used.", stacklevel=2)


def coul(qi, qj, rij):
    """Calculate Coulomb energy between two point charges.
    
    Parameters
    ----------
    qi : float
        First charge in elementary charge units (e)
    qj : float
        Second charge in elementary charge units (e)
    rij : float
        Distance between charges in Angstroms (A)
        
    Returns
    -------
    float
        Coulomb energy in electron volts (eV)
    """
    e = qqrd2e * qi * qj / rij
    return e


def real_coul(qi, qj, xi, xj, gewald, force=False):
    """Calculate real-space Coulomb energy with Ewald summation.
    
    Parameters
    ----------
    qi : float
        First charge in elementary charge units (e)
    qj : float
        Second charge in elementary charge units (e)
    xi : array_like
        Position of first charge
    xj : array_like
        Position of second charge
    gewald : float
        Ewald screening parameter
    force : bool, optional
        Whether to calculate forces, by default False
        
    Returns
    -------
    float or tuple
        If force=False, returns Coulomb energy in eV
        If force=True, returns tuple of (energy, force)
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
    """Calculator for Coulomb interactions with cutoff.
    
    This calculator computes Coulomb interactions between atoms
    using a simple cutoff scheme.
    """
    
    def __init__(self, cutoff=5.0) -> None:
        """Initialize CoulCutCalculator.
        
        Parameters
        ----------
        cutoff : float, optional
            Cutoff distance in Angstroms, by default 5.0
        """
        self.cutoff = cutoff

    def calculate(self, atoms):
        """Calculate Coulomb energy and forces.
        
        Parameters
        ----------
        atoms : ase.Atoms
            ASE Atoms object with positions and charges
            
        Returns
        -------
        tuple
            Tuple of (energy, forces) where energy is in eV
            and forces is an array in eV/A
        """
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
    """Calculator for long-range Coulomb interactions with Ewald summation.
    
    This calculator computes Coulomb interactions using Ewald summation
    with complementary error function screening.
    """
    
    def __init__(self, gewald, cutoff=5.0) -> None:
        """Initialize CoulLongCalculator.
        
        Parameters
        ----------
        gewald : float
            Ewald screening parameter
        cutoff : float, optional
            Real-space cutoff distance in Angstroms, by default 5.0
        """
        self.cutoff = cutoff
        self.gewald = gewald

    def calculate(self, atoms):
        """Calculate Coulomb energy and forces with Ewald summation.
        
        Parameters
        ----------
        atoms : ase.Atoms
            ASE Atoms object with positions and charges
            
        Returns
        -------
        tuple
            Tuple of (energy, forces) where energy is in eV
            and forces is an array in eV/A
        """
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


class DMFFPMECalculator:
    """Calculator for Coulomb interactions using Particle Mesh Ewald (PME).
    
    This calculator uses the DMFF library to compute Coulomb interactions
    with PME for efficient long-range electrostatics.
    """
    
    def __init__(self, cutoff=6.0) -> None:
        """Initialize DMFFPMECalculator.
        
        Parameters
        ----------
        cutoff : float, optional
            Real-space cutoff distance in Angstroms, by default 6.0
        """
        self.cutoff = cutoff

    def calculate(self, atoms, ethresh=1e-6) -> float:
        """Calculate Coulomb energy using PME.
        
        Parameters
        ----------
        atoms : ase.Atoms
            ASE Atoms object with positions and charges
        ethresh : float, optional
            Ewald threshold, by default 1e-6
            
        Returns
        -------
        float
            Coulomb energy in eV/particle
        """
        positions = atoms.get_positions()
        box = atoms.get_cell()
        charges = atoms.get_initial_charges()
        pairs = calculate_pairs(atoms, self.cutoff)

        box = torch.tensor(box, dtype=torch.double, requires_grad=False, device=DEVICE)
        positions = torch.tensor(
            positions, dtype=torch.double, requires_grad=True, device=DEVICE
        )
        charges = torch.tensor(
            charges, dtype=torch.double, requires_grad=False, device=DEVICE
        )
        pairs = torch.tensor(
            pairs, dtype=torch.int32, requires_grad=False, device=DEVICE
        )

        mscales = torch.tensor([0.0, 0.0, 1.0, 1.0, 1.0, 1.0], device=DEVICE)
        kappa, K = self.setup_ewald(box, ethresh)
        K1, K2, K3 = K
        pme_recip_fn = generate_pme_recip(
            Ck_fn=Ck_1,
            kappa=kappa,
            gamma=False,
            pme_order=6,
            K1=K1,
            K2=K2,
            K3=K3,
            lmax=0,
        )
        charges = torch.reshape(charges, (-1, 1))
        energy = energy_pme(
            positions,
            box,
            pairs,
            charges,
            None,
            None,
            None,
            mscales,
            None,
            None,
            None,
            pme_recip_fn,
            kappa,
            K1,
            K2,
            K3,
            0,
            False,
            True,
        )[0]
        # from kJ/mol to eV/particle
        return (energy * ENERGY_COEFF).item()

    def setup_ewald(self, box, ethresh):
        """Set up Ewald parameters for PME calculation.
        
        Parameters
        ----------
        box : array_like
            Simulation box vectors
        ethresh : float
            Ewald threshold
            
        Returns
        -------
        tuple
            Tuple of (kappa, K) where kappa is the screening parameter
            and K is the grid dimensions
        """
        kappa, K1, K2, K3 = setup_ewald_parameters(
            torch.tensor(self.cutoff), torch.tensor(ethresh), box, 0.01, "openmm"
        )
        K = (K1, K2, K3)
        return kappa, K


def calculate_pairs(atoms, rcut):
    """Calculate atom pairs within cutoff distance.
    
    Parameters
    ----------
    atoms : ase.Atoms
        ASE Atoms object
    rcut : float
        Cutoff distance in Angstroms
        
    Returns
    -------
    list
        List of atom pairs within cutoff, each as [i, j, 0]
    """
    cellpar = atoms.cell.cellpar()
    positions = atoms.get_positions()

    pairs = []
    dist_mat = distance_array(positions, positions, box=cellpar)
    for ii, dist_vec in enumerate(dist_mat):
        mask = dist_vec < rcut
        jjs = np.where(mask)[0]
        mask = jjs > ii
        jjs = jjs[mask]
        for jj in jjs:
            pairs.append([ii, jj, 0])
    return pairs
