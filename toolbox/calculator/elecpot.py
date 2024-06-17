# SPDX-License-Identifier: LGPL-3.0-or-later
"""
Reference:
- https://altafang.com/2020/10/13/calculating-1d-electrostatic-potential-profile-from-md-simulations/
- Code from Justina Moss (Leiden University, Email: j.h.moss@lic.leidenuniv.nl)
- https://github.com/SINGROUP/Potential_solver
"""

import time

import numpy as np
from ase.geometry.cell import cellpar_to_cell
from scipy import constants, integrate

EPSILON = constants.epsilon_0 / constants.elementary_charge * constants.angstrom


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
        self.int1 = integrate.cumtrapz(self.charge, self.grid, initial=0)
        self.int2 = -integrate.cumtrapz(self.int1, self.grid, initial=0) / EPSILON

    def _calculate_periodic(self, l_box):
        phi = self.solve_fft_poisson(l_box)
        return phi

    def _calculate_open(self):
        return self.int2

    def _calculate_dirichlet(self):
        pass

    def _calculate_neumann(self):
        pass

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
