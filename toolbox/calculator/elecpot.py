# SPDX-License-Identifier: LGPL-3.0-or-later
"""Electrostatic potential calculator module.

This module provides classes for calculating electrostatic potential
from charge density using various boundary conditions.

Reference:
- https://altafang.com/2020/10/13/calculating-1d-electrostatic-potential-profile-from-md-simulations/
- Code from Justina Moss (Leiden University, Email: j.h.moss@lic.leidenuniv.nl)
- https://github.com/SINGROUP/Potential_solver.
"""

import os
import time

import numpy as np
from ase import Atoms, io
from ase.geometry.cell import cellpar_to_cell
from scipy import integrate

from toolbox.utils.math import gaussian_int
from toolbox.utils.unit import EPSILON


class ElecPotentialCalculator:
    """Calculator for 1D electrostatic potential from charge density.
    
    This calculator computes electrostatic potential from charge density
    using various boundary conditions.
    
    Parameters
    ----------
    charge : array_like
        Charge density in e/Å³
    grid : array_like
        Grid points in Å
        
    Returns
    -------
    Electrostatic potential in eV
        (different from Hartree potential with a negative sign!)
    """

    def __init__(self, charge, grid) -> None:
        """Initialize ElecPotentialCalculator.
        
        Parameters
        ----------
        charge : array_like
            Charge density in e/Å³
        grid : array_like
            Grid points in Å
        """
        if len(grid) != len(charge):
            raise AttributeError("Grid and charge should have the same dimension.")
        self.grid = grid
        self.charge = charge

    def calculate(self, bc="periodic", **kwargs):
        """Calculate electrostatic potential.
        
        Parameters
        ----------
        bc : str, optional
            Boundary condition, by default "periodic"
            Options: "periodic", "open", "dirichlet", "neumann", "dip_cor"
        **kwargs
            Additional keyword arguments for specific boundary conditions
            
        Returns
        -------
        array_like
            Electrostatic potential in V
            
        Raises
        ------
        AttributeError
            If boundary condition is not supported
        """
        if (not hasattr(self, "int1")) or (not hasattr(self, "int2")):
            self._integrate()
        self.bc = bc
        try:
            self.potential = getattr(self, f"_calculate_{bc}")(**kwargs)
            return self.potential
        except AttributeError as e:
            raise AttributeError(f"Unsupported boundary condition {bc}") from e

    def _integrate(self):
        """Perform double integration of charge density."""
        self.int1 = integrate.cumulative_trapezoid(self.charge, self.grid, initial=0)
        self.int2 = (
            -integrate.cumulative_trapezoid(self.int1, self.grid, initial=0) / EPSILON
        )

    def _calculate_periodic(self, l_box):
        """Calculate potential with periodic boundary conditions.
        
        Parameters
        ----------
        l_box : float
            Box length in Å
        """
        phi = self.solve_fft_poisson(l_box)
        return phi

    def _calculate_open(self):
        """Calculate potential with open boundary conditions."""
        return self.int2

    def _calculate_dirichlet(self):
        """Calculate potential with Dirichlet boundary conditions.
        
        Raises
        ------
        NotImplementedError
            This method is not implemented
        """
        raise NotImplementedError

    def _calculate_neumann(self):
        """Calculate potential with Neumann boundary conditions.
        
        Raises
        ------
        NotImplementedError
            This method is not implemented
        """
        raise NotImplementedError

    def _calculate_dip_cor(self, cell):
        """Calculate potential with dipole correction.
        
        Parameters
        ----------
        cell : array_like
            Unit cell vectors or parameters
            
        Returns
        -------
        array_like
            Potential with dipole correction applied
        """
        if np.shape(cell) == (3,):
            cell = np.diag(cell)
        elif np.shape(cell) == (6,):
            cell = cellpar_to_cell(cell)
        elif np.shape(cell) == (3, 3):
            pass
        else:
            raise AttributeError("")

        z_max = cell[2][2]
        phi = self._calculate_periodic(z_max)
        cross_area = np.linalg.norm(np.cross(cell[0], cell[1]))
        surf_pol = np.sum(self.grid * self.charge) / cross_area
        v_cor = surf_pol / EPSILON * (self.grid / z_max - 0.5)
        return phi + v_cor

    def solve_fft_poisson(self, l_box):
        """Solve Poisson equation using FFT.
        
        Parameters
        ----------
        l_box : float
            Box length in Å
            
        Returns
        -------
        array_like
            Electrostatic potential
        """
        n_grid = len(self.grid)
        grid_spacing = l_box / n_grid

        start = time.process_time()
        charge_k_space = np.fft.fft(self.charge)
        end = time.process_time()
        print(f"| | FFT took {end - start} s")

        start = time.process_time()
        pot_k_space = np.zeros(n_grid, dtype=complex)
        k = np.fft.fftfreq(n_grid, grid_spacing)
        k_squared = k**2
        pot_k_space[1:] = charge_k_space[1:] / k_squared[1:]
        pot_k_space = pot_k_space / (4.0 * np.pi * np.pi * EPSILON)
        end = time.process_time()
        print(f"| | k-space arithmetics took {end - start} s")

        start = time.process_time()
        pot_grid = np.fft.ifftn(pot_k_space).real
        end = time.process_time()
        print(f"| | Inverse FFT took {end - start} s")
        return pot_grid


class GaussianElecPotentialCalculator:
    """Calculator for electrostatic potential from Gaussian charge distributions.
    
    This calculator computes electrostatic potential from atoms
    represented as Gaussian charge distributions.
    """
    
    def __init__(
        self,
        grids: np.ndarray,
        l_box: float,
        cross_area: float,
        spread_dict: dict[str, float] | None = None,
        charge_dict: dict[str, float] | None = None,
    ) -> None:
        """Initialize GaussianElecPotentialCalculator.
        
        Parameters
        ----------
        grids : np.ndarray
            Grid points in Å
        l_box : float
            Box length in Å
        cross_area : float
            Cross-sectional area in Å²
        spread_dict : Optional[Dict[str, float]], optional
            Dictionary of Gaussian spreads for each element, by default None
        charge_dict : Optional[Dict[str, float]], optional
            Dictionary of charges for each element, by default None
        """
        # setup grid
        self.l_box = l_box
        self.grids = grids
        self.grid_edges = np.linspace(0.0, self.l_box, len(self.grids) + 1)
        dx = np.diff(self.grid_edges)[0]
        self.grid_volume = cross_area * dx

        self.spread_dict = spread_dict
        self.charge_dict = charge_dict

    def calc_rho(
        self,
        mu: np.ndarray,
        spread: np.ndarray,
        charge: np.ndarray,
    ):
        """Calculate charge density from Gaussian distributions.
        
        Parameters
        ----------
        mu : np.ndarray
            Centers of Gaussian distributions
        spread : np.ndarray
            Spreads of Gaussian distributions
        charge : np.ndarray
            Charges of Gaussian distributions
            
        Returns
        -------
        np.ndarray
            Charge density on the grid
        """
        grid_edges = np.reshape(self.grid_edges, (1, -1))
        spread = np.reshape(spread, [-1, 1])
        mu = np.reshape(mu, [-1, 1])
        charge = np.reshape(charge, [-1, 1])

        # nat * ngrid
        out = gaussian_int(grid_edges[:, 1:], mu, spread) - gaussian_int(
            grid_edges[:, :-1], mu, spread
        )
        out = np.sum(out * charge, axis=0)
        # deal with periodic boundary
        # left image
        l_out = gaussian_int(grid_edges[:, 1:] - self.l_box, mu, spread) - gaussian_int(
            grid_edges[:, :-1] - self.l_box, mu, spread
        )
        l_out = np.sum(l_out * charge, axis=0)
        out += l_out
        # right image
        r_out = gaussian_int(grid_edges[:, 1:] + self.l_box, mu, spread) - gaussian_int(
            grid_edges[:, :-1] + self.l_box, mu, spread
        )
        r_out = np.sum(r_out * charge, axis=0)
        out += r_out

        np.testing.assert_almost_equal(out.sum(), charge.sum())

        rho = out / self.grid_volume
        return rho

    def run(
        self,
        mu: np.ndarray | None = None,
        spread: np.ndarray | None = None,
        charge: np.ndarray | None = None,
        atoms: Atoms | None = None,
    ):
        """Calculate electrostatic potential from Gaussian charge distributions.
        
        Parameters
        ----------
        mu : Optional[np.ndarray], optional
            Centers of Gaussian distributions, by default None
        spread : Optional[np.ndarray], optional
            Spreads of Gaussian distributions, by default None
        charge : Optional[np.ndarray], optional
            Charges of Gaussian distributions, by default None
        atoms : Optional[Atoms], optional
            ASE Atoms object, by default None
            
        Returns
        -------
        np.ndarray
            Electrostatic potential
        """
        if mu is None:
            mu = atoms.get_positions()[:, 2]
        if spread is None:
            assert self.spread_dict is not None, "spread_dict is not set"
            spread = np.array([self.spread_dict[s] for s in atoms.symbols])
        if charge is None:
            if self.charge_dict is not None:
                charge = np.array([self.charge_dict[s] for s in atoms.symbols])
            else:
                charge = atoms.get_initial_charges()

        rho = self.calc_rho(mu, spread, charge)
        calculator = ElecPotentialCalculator(rho, self.grids)
        phi = calculator.calculate(l_box=self.l_box)
        return phi


class WannierHartreePotentialCalculator(GaussianElecPotentialCalculator):
    """Calculator for Hartree potential from Wannier functions.
    
    This calculator extends GaussianElecPotentialCalculator to
    compute Hartree potential from Wannier functions.
    """
    
    def __init__(
        self,
        grids: np.ndarray,
        l_box: float,
        cross_area: float,
        spread_dict: dict[str, float] | None = None,
        charge_dict: dict[str, float] | None = None,
    ) -> None:
        """Initialize WannierHartreePotentialCalculator.
        
        Parameters
        ----------
        grids : np.ndarray
            Grid points in Å
        l_box : float
            Box length in Å
        cross_area : float
            Cross-sectional area in Å²
        spread_dict : Optional[Dict[str, float]], optional
            Dictionary of Gaussian spreads for each element, by default None
        charge_dict : Optional[Dict[str, float]], optional
            Dictionary of charges for each element, by default None
        """
        super().__init__(grids, l_box, cross_area, spread_dict, charge_dict)

    def run(
        self,
        mu: np.ndarray | None = None,
        spread: np.ndarray | None = None,
        charge: np.ndarray | None = None,
        atoms: Atoms | None = None,
        dname: str | None = None,
        fname_coord: str | None = "coord.xyz",
        fname_wannier: str | None = "wannier.xyz",
    ):
        """Calculate Hartree potential from Wannier functions.
        
        Parameters
        ----------
        mu : Optional[np.ndarray], optional
            Centers of Gaussian distributions, by default None
        spread : Optional[np.ndarray], optional
            Spreads of Gaussian distributions, by default None
        charge : Optional[np.ndarray], optional
            Charges of Gaussian distributions, by default None
        atoms : Optional[Atoms], optional
            ASE Atoms object, by default None
        dname : Optional[str], optional
            Directory name containing coordinate files, by default None
        fname_coord : Optional[str], optional
            Coordinate file name, by default "coord.xyz"
        fname_wannier : Optional[str], optional
            Wannier coordinate file name, by default "wannier.xyz"
            
        Returns
        -------
        np.ndarray
            Negative of the electrostatic potential (Hartree potential)
        """
        if dname is not None:
            _atoms = io.read(os.path.join(dname, fname_coord))
            wannier_atoms = io.read(os.path.join(dname, fname_wannier))
            atoms = _atoms + wannier_atoms
            atoms.set_pbc(True)
            atoms.wrap()

        phi = super().run(mu, spread, charge, atoms)
        return -phi
