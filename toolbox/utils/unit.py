# SPDX-License-Identifier: LGPL-3.0-or-later
"""
Ref:
    https://docs.scipy.org/doc/scipy/reference/constants.html
"""

from scipy import constants

# AU
AU_TO_ANG = 5.29177208590000e-01
AU_TO_EV = 2.72113838565563e01
AU_TO_EV_EVERY_ANG = AU_TO_EV / AU_TO_ANG
DEBYE_TO_EA = 1.0 / 4.8

# epsilon [e/(V Angstrom)]
EPSILON = constants.epsilon_0 / constants.elementary_charge * constants.angstrom
