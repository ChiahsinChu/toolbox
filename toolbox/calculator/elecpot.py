"""
Reference:
- https://altafang.com/2020/10/13/calculating-1d-electrostatic-potential-profile-from-md-simulations/
- Code from Justina Moss (Leiden University, Email: j.h.moss@lic.leidenuniv.nl)
"""
import numpy as np
from scipy import integrate
from ase.geometry.cell import cellpar_to_cell

from ..utils.unit import *

_EPSILON = VAC_PERMITTIVITY / UNIT_CHARGE * ANG_TO_M


class ElecPotentialCalculator:
    """
    Parameters
    ----------
    charge:
        charge density in e/A^3
    grid:
        grid in A
    """
    def __init__(self, charge, grid) -> None:
        if len(grid) != len(charge):
            raise AttributeError(
                "Grid and charge should have the same dimension.")
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
        self.int2 = -integrate.cumtrapz(self.int1, self.grid,
                                        initial=0) / _EPSILON

    def _calculate_periodic(self):
        phi = self.int2 - ((self.int2[-1] - self.int2[0]) /
                           (self.grid[-1] - self.grid[0]) *
                           (self.grid - self.grid[0]) + self.int2[0])
        return phi

    def _calculate_open(self):
        return self.int2

    def _calculate_dirichlet(self):
        pass

    def _calculate_dip_cor(self, cell):
        if np.shape(cell) == (3, ):
            cell = np.diag(cell)
        elif np.shape(cell) == (6, ):
            cell = cellpar_to_cell(cell)
        elif np.shape(cell) == (3, 3):
            pass
        else:
            raise AttributeError("")

        phi = self._calculate_periodic()
        z_max = cell[2][2]
        cross_area = np.linalg.norm(np.cross(cell[0], cell[1]))
        surf_pol = np.sum(self.grid * self.charge) / cross_area
        v_cor = surf_pol / _EPSILON * (self.grid / z_max - 0.5)
        return phi + v_cor