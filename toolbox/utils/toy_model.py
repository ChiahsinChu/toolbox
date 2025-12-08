# SPDX-License-Identifier: LGPL-3.0-or-later

from scipy import constants


class ParallelPlateCapacitor:
    """Class for parallel plate capacitor calculations.
    
    This class models a parallel plate capacitor system
    with configurable permittivity, plate distance, and capacitance.
    
    Parameters
    ----------
    epsilon_r : float
        Relative permittivity
    d : float
        Distance between plates in Angstroms
    capacitance : float
        Capacitance in μF/cm²
        
    Examples
    --------
    >>> obj = ParallelPlateCapacitor(epsilon_r=10.0, d=4.0)
    >>> print(obj)
    >>> obj = ParallelPlateCapacitor(capacitance=20.0, d=4.0)
    >>> print(obj)
    """

    def __init__(
        self,
        epsilon_r: float | None = None,
        d: float | None = None,
        capacitance: float | None = None,
    ):
        """Initialize ParallelPlateCapacitor.
        
        Parameters
        ----------
        epsilon_r : Optional[float], optional
            Relative permittivity, by default None
        d : Optional[float], optional
            Distance between plates in Angstroms, by default None
        capacitance : Optional[float], optional
            Capacitance in μF/cm², by default None
        """
        # conversion factor from F/m² to μF/cm²
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
        """Get capacitance in μF/cm²."""
        try:
            return self._capacitance
        except AttributeError:
            return self._epsilon_r * constants.epsilon_0 / self._d * self._coeff

    @property
    def epsilon_r(self):
        """Get relative permittivity."""
        try:
            return self._epsilon_r
        except AttributeError:
            return self._capacitance * self._d / (constants.epsilon_0 * self._coeff)

    @property
    def d_angstrom(self):
        """Get plate distance in Angstroms."""
        return self._d / constants.angstrom

    def __repr__(self):
        return f"ParallelPlateCapacitor(epsilon_r={self.epsilon_r}, d={self.d_angstrom} angstrom, capacitance={self.capacitance} muF/cm^2)"
