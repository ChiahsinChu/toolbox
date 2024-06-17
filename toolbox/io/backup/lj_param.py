# SPDX-License-Identifier: LGPL-3.0-or-later
"""
L-J parameters from https://pubs.acs.org/doi/full/10.1021/ct500918t
"""

from toolbox.utils.utils import save_dict

param_dict = {
    "Li": {"epsilon": 0.00274091, "sigma": 1.258},
    "Na": {"epsilon": 0.02639002, "sigma": 1.454},
    "K": {"epsilon": 0.12693448, "sigma": 1.683},
    "Rb": {"epsilon": 0.20665151, "sigma": 1.792},
    "Cs": {"epsilon": 0.34673208, "sigma": 1.953},
    "F": {"epsilon": 0.22878796, "sigma": 1.819},
    "Cl": {"epsilon": 0.64367011, "sigma": 2.308},
    "Br": {"epsilon": 0.74435812, "sigma": 2.470},
    "I": {"epsilon": 0.86877007, "sigma": 2.770},
}

new_param_dict = {}
# from kcal/mol to eV
factor = 0.043361254529175

for k, v in param_dict.items():
    new_param_dict[k] = {
        "epsilon": v["epsilon"] * factor,
        "sigma": ((v["sigma"] * 2) ** 12 / 4.0) ** (1.0 / 12.0),
    }
# print(new_param_dict)
save_dict(new_param_dict, "lj.json")
