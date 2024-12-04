# SPDX-License-Identifier: LGPL-3.0-or-later
"""
Reference:
- https://altafang.com/2020/10/13/calculating-1d-electrostatic-potential-profile-from-md-simulations/
- Code from Justina Moss (Leiden University, Email: j.h.moss@lic.leidenuniv.nl)
- https://github.com/SINGROUP/Potential_solver
"""

import os
import time
from typing import Dict, Optional

import numpy as np
from ase import Atoms, io
from ase.geometry.cell import cellpar_to_cell
from scipy import integrate

from toolbox.utils.math import gaussian_int
from toolbox.utils.unit import EPSILON


class ElecPotentialCalculator:
    """
    Parameters
    ----------
    charge:
        charge density in e/A^3
    grid:
        grid in A

    Return
    ------
    Electrostatic potential in eV
    (different from Hartree potential with a negative sign!)
    """

    def __init__(self, charge, grid) -> None:
        if len(grid) != len(charge):
            raise AttributeError("Grid and charge should have the same dimension.")
        self.grid = grid
        self.charge = charge

    def calculate(self, bc="periodic", **kwargs):
        """
        Return:
            electrostatic potential in V
        """

        if (not hasattr(self, "int1")) or (not hasattr(self, "int2")):
            self._integrate()
        self.bc = bc
        try:
            self.potential = getattr(self, "_calculate_%s" % bc)(**kwargs)
            return self.potential
        except:
            raise AttributeError("Unsupported boundary condition %s" % bc)

    def _integrate(self):
        self.int1 = integrate.cumulative_trapezoid(self.charge, self.grid, initial=0)
        self.int2 = (
            -integrate.cumulative_trapezoid(self.int1, self.grid, initial=0) / EPSILON
        )

    def _calculate_periodic(self, l_box):
        phi = self.solve_fft_poisson(l_box)
        return phi

    def _calculate_open(self):
        return self.int2

    def _calculate_dirichlet(self):
        raise NotImplementedError

    def _calculate_neumann(self):
        raise NotImplementedError

    def _calculate_dip_cor(self, cell):
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
        n_grid = len(self.grid)
        grid_spacing = l_box / n_grid

        start = time.process_time()
        charge_k_space = np.fft.fft(self.charge)
        end = time.process_time()
        print("| | FFT took {} s".format(end - start))

        start = time.process_time()
        pot_k_space = np.zeros(n_grid, dtype=complex)
        k = np.fft.fftfreq(n_grid, grid_spacing)
        k_squared = k**2
        pot_k_space[1:] = charge_k_space[1:] / k_squared[1:]
        pot_k_space = pot_k_space / (4.0 * np.pi * np.pi * EPSILON)
        end = time.process_time()
        print("| | k-space arithmetics took {} s".format(end - start))

        start = time.process_time()
        pot_grid = np.fft.ifftn(pot_k_space).real
        end = time.process_time()
        print("| | Inverse FFT took {} s".format(end - start))
        return pot_grid


class GaussianElecPotentialCalculator:
    def __init__(
        self,
        grids: np.ndarray,
        l_box: float,
        cross_area: float,
        spread_dict: Optional[Dict[str, float]] = None,
        charge_dict: Optional[Dict[str, float]] = None,
    ) -> None:
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
        mu: Optional[np.ndarray] = None,
        spread: Optional[np.ndarray] = None,
        charge: Optional[np.ndarray] = None,
        atoms: Optional[Atoms] = None,
    ):
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
    def __init__(
        self,
        grids: np.ndarray,
        l_box: float,
        cross_area: float,
        spread_dict: Optional[Dict[str, float]] = None,
        charge_dict: Optional[Dict[str, float]] = None,
    ) -> None:
        super().__init__(grids, l_box, cross_area, spread_dict, charge_dict)

    def run(
        self,
        mu: Optional[np.ndarray] = None,
        spread: Optional[np.ndarray] = None,
        charge: Optional[np.ndarray] = None,
        atoms: Optional[Atoms] = None,
        dname: Optional[str] = None,
        fname_coord: Optional[str] = "coord.xyz",
        fname_wannier: Optional[str] = "wannier.xyz",
    ):
        if dname is not None:
            _atoms = io.read(os.path.join(dname, fname_coord))
            wannier_atoms = io.read(os.path.join(dname, fname_wannier))
            atoms = _atoms + wannier_atoms
            atoms.set_pbc(True)
            atoms.wrap()

        phi = super().run(mu, spread, charge, atoms)
        return -phi
