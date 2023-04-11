"""
source: https://gitlab.com/jmargraf/kqeq/-/blob/master/kqeq/data.py
"""
uff_xi_vdw = {
    'H': 2.886,
    'He': 2.362,
    'Li': 2.451,
    'Be': 2.745,
    'B': 4.083,
    'C': 3.851,
    'N': 3.66,
    'O': 3.5,
    'F': 3.364,
    'Ne': 3.243,
    'Na': 2.983,
    'Mg': 3.021,
    'Al': 4.499,
    'Si': 4.295,
    'P': 4.147,
    'S': 4.035,
    'Cl': 3.947,
    'Ar': 3.868,
    'K': 3.812,
    'Ca': 3.399,
    'Sc': 3.295,
    'Ti': 3.175,
    'V': 3.144,
    'Cr': 3.023,
    'Mn': 2.961,
    'Fe': 2.912,
    'Co': 2.872,
    'Ni': 2.834,
    'Cu': 3.495,
    'Zn': 2.763,
    'Ga': 4.383,
    'Ge': 4.28,
    'As': 4.23,
    'Se': 4.205,
    'Br': 4.189,
    'Kr': 4.141,
    'Rb': 4.114,
    'Sr': 3.641,
    'Y': 3.345,
    'Zr': 3.124,
    'Nb': 3.165,
    'Mo': 3.052,
    'Tc': 2.998,
    'Ru': 2.963,
    'Rh': 2.929,
    'Pd': 2.899,
    'Ag': 3.148,
    'Cd': 2.848,
    'In': 4.463,
    'Sn': 4.392,
    'Sb': 4.42,
    'Te': 4.47,
    'I': 4.5,
    'Xe': 4.404,
    'Cs': 4.517,
    'Ba': 3.703,
    'La': 3.522,
    'Ce': 3.556,
    'Pr': 3.606,
    'Nd': 3.575,
    'Pm': 3.547,
    'Sm': 3.52,
    'Eu': 3.493,
    'Gd': 3.368,
    'Tb': 3.451,
    'Dy': 3.428,
    'Ho': 3.409,
    'Er': 3.391,
    'Tm': 3.374,
    'Yb': 3.355,
    'Lu': 3.64,
    'Hf': 3.141,
    'Ta': 3.17,
    'W': 3.069,
    'Re': 2.954,
    'Os': 3.12,
    'Ir': 2.84,
    'Pt': 2.754,
    'Au': 3.293,
    'Hg': 2.705,
    'Tl': 4.347,
    'Pb': 4.297,
    'Bi': 4.37,
    'Po': 4.709,
    'At': 4.75,
    'Rn': 4.765,
    'Fr': 4.9,
    'Ra': 3.677,
    'Ac': 3.478,
    'Th': 3.396,
    'Pa': 3.424,
    'U': 3.395,
    'Np': 3.424,
    'Pu': 3.424,
    'Am': 3.381,
    'Cm': 3.326,
    'Bk': 3.339,
    'Cf': 3.313,
    'Es': 3.299,
    'Fm': 3.286,
    'Md': 3.274,
    'No': 3.248,
    'Lr': 3.236
}

uff_Di_vdw = {
    'H': 0.044,
    'He': 0.056,
    'Li': 0.025,
    'Be': 0.085,
    'B': 0.18,
    'C': 0.105,
    'N': 0.069,
    'O': 0.06,
    'F': 0.05,
    'Ne': 0.042,
    'Na': 0.03,
    'Mg': 0.111,
    'Al': 0.505,
    'Si': 0.402,
    'P': 0.305,
    'S': 0.274,
    'Cl': 0.227,
    'Ar': 0.185,
    'K': 0.035,
    'Ca': 0.238,
    'Sc': 0.019,
    'Ti': 0.017,
    'V': 0.016,
    'Cr': 0.015,
    'Mn': 0.013,
    'Fe': 0.013,
    'Co': 0.014,
    'Ni': 0.015,
    'Cu': 0.005,
    'Zn': 0.124,
    'Ga': 0.415,
    'Ge': 0.379,
    'As': 0.309,
    'Se': 0.291,
    'Br': 0.251,
    'Kr': 0.22,
    'Rb': 0.04,
    'Sr': 0.235,
    'Y': 0.072,
    'Zr': 0.069,
    'Nb': 0.059,
    'Mo': 0.056,
    'Tc': 0.048,
    'Ru': 0.056,
    'Rh': 0.053,
    'Pd': 0.048,
    'Ag': 0.036,
    'Cd': 0.228,
    'In': 0.599,
    'Sn': 0.567,
    'Sb': 0.449,
    'Te': 0.398,
    'I': 0.339,
    'Xe': 0.332,
    'Cs': 0.045,
    'Ba': 0.364,
    'La': 0.017,
    'Ce': 0.013,
    'Pr': 0.01,
    'Nd': 0.01,
    'Pm': 0.009,
    'Sm': 0.008,
    'Eu': 0.008,
    'Gd': 0.009,
    'Tb': 0.007,
    'Dy': 0.007,
    'Ho': 0.007,
    'Er': 0.007,
    'Tm': 0.006,
    'Yb': 0.228,
    'Lu': 0.041,
    'Hf': 0.072,
    'Ta': 0.081,
    'W': 0.067,
    'Re': 0.066,
    'Os': 0.037,
    'Ir': 0.073,
    'Pt': 0.08,
    'Au': 0.039,
    'Hg': 0.385,
    'Tl': 0.68,
    'Pb': 0.663,
    'Bi': 0.518,
    'Po': 0.325,
    'At': 0.284,
    'Rn': 0.248,
    'Fr': 0.05,
    'Ra': 0.404,
    'Ac': 0.033,
    'Th': 0.026,
    'Pa': 0.022,
    'U': 0.022,
    'Np': 0.019,
    'Pu': 0.016,
    'Am': 0.014,
    'Cm': 0.013,
    'Bk': 0.013,
    'Cf': 0.013,
    'Es': 0.012,
    'Fm': 0.012,
    'Md': 0.011,
    'No': 0.011,
    'Lr': 0.011
}

uff_zeta_vdw = {
    'H': 12.0,
    'He': 15.24,
    'Li': 12.0,
    'Be': 12.0,
    'B': 12.052,
    'C': 12.73,
    'N': 13.407,
    'O': 14.085,
    'F': 14.762,
    'Ne': 15.44,
    'Na': 12.0,
    'Mg': 12.0,
    'Al': 11.278,
    'Si': 12.175,
    'P': 13.072,
    'S': 13.969,
    'Cl': 14.866,
    'Ar': 15.763,
    'K': 12.0,
    'Ca': 12.0,
    'Sc': 12.0,
    'Ti': 12.0,
    'V': 12.0,
    'Cr': 12.0,
    'Mn': 12.0,
    'Fe': 12.0,
    'Co': 12.0,
    'Ni': 12.0,
    'Cu': 12.0,
    'Zn': 12.0,
    'Ga': 11.0,
    'Ge': 12.0,
    'As': 13.0,
    'Se': 14.0,
    'Br': 15.0,
    'Kr': 16.0,
    'Rb': 12.0,
    'Sr': 12.0,
    'Y': 12.0,
    'Zr': 12.0,
    'Nb': 12.0,
    'Mo': 12.0,
    'Tc': 12.0,
    'Ru': 12.0,
    'Rh': 12.0,
    'Pd': 12.0,
    'Ag': 12.0,
    'Cd': 12.0,
    'In': 11.0,
    'Sn': 12.0,
    'Sb': 13.0,
    'Te': 14.0,
    'I': 15.0,
    'Xe': 12.0,
    'Cs': 12.0,
    'Ba': 12.0,
    'La': 12.0,
    'Ce': 12.0,
    'Pr': 12.0,
    'Nd': 12.0,
    'Pm': 12.0,
    'Sm': 12.0,
    'Eu': 12.0,
    'Gd': 12.0,
    'Tb': 12.0,
    'Dy': 12.0,
    'Ho': 12.0,
    'Er': 12.0,
    'Tm': 12.0,
    'Yb': 12.0,
    'Lu': 12.0,
    'Hf': 12.0,
    'Ta': 12.0,
    'W': 12.0,
    'Re': 12.0,
    'Os': 12.0,
    'Ir': 12.0,
    'Pt': 12.0,
    'Au': 12.0,
    'Hg': 12.0,
    'Tl': 11.0,
    'Pb': 12.0,
    'Bi': 13.0,
    'Po': 14.0,
    'At': 15.0,
    'Rn': 16.0,
    'Fr': 12.0,
    'Ra': 12.0,
    'Ac': 12.0,
    'Th': 12.0,
    'Pa': 12.0,
    'U': 12.0,
    'Np': 12.0,
    'Pu': 12.0,
    'Am': 12.0,
    'Cm': 12.0,
    'Bk': 12.0,
    'Cf': 12.0,
    'Es': 12.0,
    'Fm': 12.0,
    'Md': 12.0,
    'No': 12.0,
    'Lr': 12.0
}

uff_Xi_qeq = {
    'H': 4.528,
    'He': 9.66,
    'Li': 3.006,
    'Be': 4.877,
    'B': 5.11,
    'C': 5.343,
    'N': 6.899,
    'O': 8.741,
    'F': 10.874,
    'Ne': 11.04,
    'Na': 2.843,
    'Mg': 3.951,
    'Al': 4.06,
    'Si': 4.168,
    'P': 5.463,
    'S': 6.928,
    'Cl': 8.564,
    'Ar': 9.465,
    'K': 2.421,
    'Ca': 3.231,
    'Sc': 3.395,
    'Ti': 3.47,
    'V': 3.65,
    'Cr': 3.415,
    'Mn': 3.325,
    'Fe': 3.76,
    'Co': 4.105,
    'Ni': 4.465,
    'Cu': 4.2,
    'Zn': 5.106,
    'Ga': 3.641,
    'Ge': 4.051,
    'As': 5.188,
    'Se': 6.428,
    'Br': 7.79,
    'Kr': 8.505,
    'Rb': 2.331,
    'Sr': 3.024,
    'Y': 3.83,
    'Zr': 3.4,
    'Nb': 3.55,
    'Mo': 3.465,
    'Tc': 3.29,
    'Ru': 3.575,
    'Rh': 3.975,
    'Pd': 4.32,
    'Ag': 4.436,
    'Cd': 5.034,
    'In': 3.506,
    'Sn': 3.987,
    'Sb': 4.899,
    'Te': 5.816,
    'I': 6.822,
    'Xe': 7.595,
    'Cs': 2.183,
    'Ba': 2.814,
    'La': 2.8355,
    'Ce': 2.774,
    'Pr': 2.858,
    'Nd': 2.8685,
    'Pm': 2.881,
    'Sm': 2.9115,
    'Eu': 2.8785,
    'Gd': 3.1665,
    'Tb': 3.018,
    'Dy': 3.0555,
    'Ho': 3.127,
    'Er': 3.1865,
    'Tm': 3.2514,
    'Yb': 3.2889,
    'Lu': 2.9629,
    'Hf': 3.7,
    'Ta': 5.1,
    'W': 4.63,
    'Re': 3.96,
    'Os': 5.14,
    'Ir': 5.0,
    'Pt': 4.79,
    'Au': 4.894,
    'Hg': 6.27,
    'Tl': 3.2,
    'Pb': 3.9,
    'Bi': 4.69,
    'Po': 4.21,
    'At': 4.75,
    'Rn': 5.37,
    'Fr': 2.0,
    'Ra': 2.843,
    'Ac': 2.835,
    'Th': 3.175,
    'Pa': 2.985,
    'U': 3.341,
    'Np': 3.549,
    'Pu': 3.243,
    'Am': 2.9895,
    'Cm': 2.8315,
    'Bk': 3.1935,
    'Cf': 3.197,
    'Es': 3.333,
    'Fm': 3.4,
    'Md': 3.47,
    'No': 3.475,
    'Lr': 3.5
}

uff_hard_qeq = {
    'H': 6.9452,
    'He': 14.92,
    'Li': 2.386,
    'Be': 4.443,
    'B': 4.75,
    'C': 5.063,
    'N': 5.88,
    'O': 6.682,
    'F': 7.474,
    'Ne': 10.55,
    'Na': 2.296,
    'Mg': 3.693,
    'Al': 3.59,
    'Si': 3.487,
    'P': 4.0,
    'S': 4.486,
    'Cl': 4.946,
    'Ar': 6.355,
    'K': 1.92,
    'Ca': 2.88,
    'Sc': 3.08,
    'Ti': 3.38,
    'V': 3.41,
    'Cr': 3.865,
    'Mn': 4.105,
    'Fe': 4.14,
    'Co': 4.175,
    'Ni': 4.205,
    'Cu': 4.22,
    'Zn': 4.285,
    'Ga': 3.16,
    'Ge': 3.438,
    'As': 3.809,
    'Se': 4.131,
    'Br': 4.425,
    'Kr': 5.715,
    'Rb': 1.846,
    'Sr': 2.44,
    'Y': 2.81,
    'Zr': 3.55,
    'Nb': 3.38,
    'Mo': 3.755,
    'Tc': 3.99,
    'Ru': 4.015,
    'Rh': 4.005,
    'Pd': 4.0,
    'Ag': 3.134,
    'Cd': 3.957,
    'In': 2.896,
    'Sn': 3.124,
    'Sb': 3.342,
    'Te': 3.526,
    'I': 3.762,
    'Xe': 4.975,
    'Cs': 1.711,
    'Ba': 2.396,
    'La': 2.7415,
    'Ce': 2.692,
    'Pr': 2.564,
    'Nd': 2.6205,
    'Pm': 2.673,
    'Sm': 2.7195,
    'Eu': 2.7875,
    'Gd': 2.9745,
    'Tb': 2.834,
    'Dy': 2.8715,
    'Ho': 2.891,
    'Er': 2.9145,
    'Tm': 2.9329,
    'Yb': 2.965,
    'Lu': 2.4629,
    'Hf': 3.4,
    'Ta': 2.85,
    'W': 3.31,
    'Re': 3.92,
    'Os': 3.63,
    'Ir': 4.0,
    'Pt': 4.43,
    'Au': 2.586,
    'Hg': 4.16,
    'Tl': 2.9,
    'Pb': 3.53,
    'Bi': 3.74,
    'Po': 4.21,
    'At': 4.75,
    'Rn': 5.37,
    'Fr': 2.0,
    'Ra': 2.434,
    'Ac': 2.835,
    'Th': 2.905,
    'Pa': 2.905,
    'U': 2.853,
    'Np': 2.717,
    'Pu': 2.819,
    'Am': 3.0035,
    'Cm': 3.1895,
    'Bk': 3.0355,
    'Cf': 3.101,
    'Es': 3.089,
    'Fm': 3.1,
    'Md': 3.11,
    'No': 3.175,
    'Lr': 3.2
}

uff_radius_qeq = {
    'H': 0.371,
    'He': 1.3,
    'Li': 1.557,
    'Be': 1.24,
    'B': 0.822,
    'C': 0.759,
    'N': 0.715,
    'O': 0.669,
    'F': 0.706,
    'Ne': 1.768,
    'Na': 2.085,
    'Mg': 1.5,
    'Al': 1.201,
    'Si': 1.176,
    'P': 1.102,
    'S': 1.047,
    'Cl': 0.994,
    'Ar': 2.108,
    'K': 2.586,
    'Ca': 2.0,
    'Sc': 1.75,
    'Ti': 1.607,
    'V': 1.47,
    'Cr': 1.402,
    'Mn': 1.533,
    'Fe': 1.393,
    'Co': 1.406,
    'Ni': 1.398,
    'Cu': 1.434,
    'Zn': 1.4,
    'Ga': 1.211,
    'Ge': 1.189,
    'As': 1.204,
    'Se': 1.224,
    'Br': 1.141,
    'Kr': 2.27,
    'Rb': 2.77,
    'Sr': 2.415,
    'Y': 1.998,
    'Zr': 1.758,
    'Nb': 1.603,
    'Mo': 1.53,
    'Tc': 1.5,
    'Ru': 1.5,
    'Rh': 1.509,
    'Pd': 1.544,
    'Ag': 1.622,
    'Cd': 1.6,
    'In': 1.404,
    'Sn': 1.354,
    'Sb': 1.404,
    'Te': 1.38,
    'I': 1.333,
    'Xe': 2.459,
    'Cs': 2.984,
    'Ba': 2.442,
    'La': 2.071,
    'Ce': 1.925,
    'Pr': 2.007,
    'Nd': 2.007,
    'Pm': 2.0,
    'Sm': 1.978,
    'Eu': 2.227,
    'Gd': 1.968,
    'Tb': 1.954,
    'Dy': 1.934,
    'Ho': 1.925,
    'Er': 1.915,
    'Tm': 2.0,
    'Yb': 2.158,
    'Lu': 1.896,
    'Hf': 1.759,
    'Ta': 1.605,
    'W': 1.538,
    'Re': 1.6,
    'Os': 1.7,
    'Ir': 1.866,
    'Pt': 1.557,
    'Au': 1.618,
    'Hg': 1.6,
    'Tl': 1.53,
    'Pb': 1.444,
    'Bi': 1.514,
    'Po': 1.48,
    'At': 1.47,
    'Rn': 2.2,
    'Fr': 2.3,
    'Ra': 2.2,
    'Ac': 2.108,
    'Th': 2.018,
    'Pa': 1.8,
    'U': 1.713,
    'Np': 1.8,
    'Pu': 1.84,
    'Am': 1.942,
    'Cm': 1.9,
    'Bk': 1.9,
    'Cf': 1.9,
    'Es': 1.9,
    'Fm': 1.9,
    'Md': 1.9,
    'No': 1.9,
    'Lr': 1.9
}