# SPDX-License-Identifier: LGPL-3.0-or-later
from typing import Optional

from scipy import constants


class ParallelPlateCapacitor:
    """
    Parallel plate capacitor

    Parameters
    ----------
    epsilon_r : float
        Relative permittivity
    d : float
        Distance between plates [angstrom]
    capacitance : float
        Capacitance [muF/cm^2]

    Examples
    --------
    obj = ParallelPlateCapacitor(epsilon_r=10.0, d=4.0)
    print(obj)
    obj = ParallelPlateCapacitor(capacitance=20.0, d=4.0)
    print(obj)
    """

    def __init__(
        self,
        epsilon_r: Optional[float] = None,
        d: Optional[float] = None,
        capacitance: Optional[float] = None,
    ):
        # conversion factor from F/m^2 to muF/cm^2
        self._coeff = constants.centi**2 / constants.micro

        if epsilon_r is not None:
            self._epsilon_r = epsilon_r

        if d is not None:
            # d [m]
            self._d = d * constants.angstrom

        if capacitance is not None:
            self._capacitance = capacitance

    @property
    def capacitance(self):
        try:
            return self._capacitance
        except AttributeError:
            return self._epsilon_r * constants.epsilon_0 / self._d * self._coeff

    @property
    def epsilon_r(self):
        try:
            return self._epsilon_r
        except AttributeError:
            return self._capacitance * self._d / (constants.epsilon_0 * self._coeff)

    @property
    def d_angstrom(self):
        return self._d / constants.angstrom

    def __repr__(self):
        return f"ParallelPlateCapacitor(epsilon_r={self.epsilon_r}, d={self.d_angstrom} angstrom, capacitance={self.capacitance} muF/cm^2)"
