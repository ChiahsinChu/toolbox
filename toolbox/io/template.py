# SPDX-License-Identifier: LGPL-3.0-or-later
"""
CP2K input files:
- energy (metal with smearing)
- geo_opt
- cell_opt
- bomd
- qmmm
"""

import copy

from ..utils.utils import update_dict

# cp2k input templates
cp2k_default_input = {
    "energy": {
        "FORCE_EVAL": {
            "METHOD": "QS",
            "DFT": {
                "BASIS_SET_FILE_NAME": [
                    "BASIS_MOLOPT",
                    "BASIS_ADMM",
                    "BASIS_ADMM_MOLOPT",
                    "BASIS_MOLOPT-HSE06",
                ],
                "POTENTIAL_FILE_NAME": "GTH_POTENTIALS",
                "MGRID": {"CUTOFF": 400},
                "QS": {"EPS_DEFAULT": 1.0e-13},
                "SCF": {
                    "SCF_GUESS": "RESTART",
                    "EPS_SCF": 1.0e-6,
                    "MAX_SCF": 500,
                    "ADDED_MOS": 500,
                    "CHOLESKY": "INVERSE",
                    "SMEAR": {
                        "_": "ON",
                        "METHOD": "FERMI_DIRAC",
                        "ELECTRONIC_TEMPERATURE": 300,
                    },
                    "DIAGONALIZATION": {"ALGORITHM": "STANDARD"},
                    "MIXING": {
                        "METHOD": "BROYDEN_MIXING",
                        "ALPHA": 0.3,
                        "BETA": 1.5,
                        "NBUFFER": 8,
                    },
                },
                "XC": {
                    "XC_FUNCTIONAL": {"_": "PBE"},
                    "vdW_POTENTIAL": {
                        "DISPERSION_FUNCTIONAL": "PAIR_POTENTIAL",
                        "PAIR_POTENTIAL": {
                            "TYPE": "DFTD3",
                            "PARAMETER_FILE_NAME": "dftd3.dat",
                            "REFERENCE_FUNCTIONAL": "PBE",
                        },
                    },
                },
            },
            "SUBSYS": {
                "TOPOLOGY": {
                    "COORD_FILE_FORMAT": "XYZ",
                    "COORD_FILE_NAME": "coord.xyz",
                },
                "KIND": [
                    {
                        "_": "O",
                        "POTENTIAL": "GTH-PBE-q6",
                        "BASIS_SET": "DZVP-MOLOPT-SR-GTH",
                    },
                    {
                        "_": "H",
                        "POTENTIAL": "GTH-PBE-q1",
                        "BASIS_SET": "DZVP-MOLOPT-SR-GTH",
                    },
                    {
                        "_": "Pt",
                        "POTENTIAL": "GTH-PBE-q10",
                        "BASIS_SET": "DZVP-A5-Q10-323-MOL-T1-DERIVED_SET-1",
                    },
                    {
                        "_": "Ag",
                        "POTENTIAL": "GTH-PBE-q11",
                        "BASIS_SET": "DZVP-MOLOPT-SR-GTH",
                    },
                    {
                        "_": "Na",
                        "POTENTIAL": "GTH-PBE-q9",
                        "BASIS_SET": "DZVP-MOLOPT-SR-GTH",
                    },
                    {
                        "_": "K",
                        "POTENTIAL": "GTH-PBE-q9",
                        "BASIS_SET": "DZVP-MOLOPT-SR-GTH",
                    },
                    {
                        "_": "Li",
                        "POTENTIAL": "GTH-PBE-q3",
                        "BASIS_SET": "DZVP-MOLOPT-SR-GTH",
                    },
                    {
                        "_": "C",
                        "POTENTIAL": "GTH-PBE-q4",
                        "BASIS_SET": "DZVP-MOLOPT-SR-GTH",
                    },
                    {
                        "_": "N",
                        "POTENTIAL": "GTH-PBE-q5",
                        "BASIS_SET": "DZVP-MOLOPT-SR-GTH",
                    },
                    {
                        "_": "Cl",
                        "POTENTIAL": "GTH-PBE-q7",
                        "BASIS_SET": "DZVP-MOLOPT-SR-GTH",
                    },
                    {
                        "_": "F",
                        "POTENTIAL": "GTH-PBE-q7",
                        "BASIS_SET": "DZVP-MOLOPT-SR-GTH",
                    },
                    {
                        "_": "Mg",
                        "POTENTIAL": "GTH-PBE-q10",
                        "BASIS_SET": "DZVP-MOLOPT-SR-GTH",
                    },
                    {
                        "_": "Al",
                        "POTENTIAL": "GTH-PBE-q3",
                        "BASIS_SET": "DZVP-MOLOPT-SR-GTH",
                    },
                ],
            },
            "PRINT": {"FORCES": {"_": "ON"}},
        },
        "GLOBAL": {"PROJECT": "cp2k"},
    }
}

cp2k_default_input.update(
    {
        "geo_opt": {
            "FORCE_EVAL": cp2k_default_input["energy"]["FORCE_EVAL"],
            "GLOBAL": {"PROJECT": "cp2k", "RUN_TYPE": "GEO_OPT"},
            "MOTION": {
                "GEO_OPT": {
                    "TYPE": "MINIMIZATION",
                    "OPTIMIZER": "BFGS",
                    "MAX_ITER": 200,
                    "MAX_FORCE": 4.5e-5,
                },
                "PRINT": {
                    "TRAJECTORY": {"EACH": {"GEO_OPT": 1}},
                    "VELOCITIES": {"EACH": {"GEO_OPT": 1}},
                    "FORCES": {"_": "ON"},
                    "RESTART_HISTORY": {},
                    "RESTART": {"BACKUP_COPIES": 3},
                },
            },
        }
    }
)

# >>>>>>>>>>>>>>>>>>>> CELL_OPT >>>>>>>>>>>>>>>>>>>>
cp2k_default_input["cell_opt"] = copy.deepcopy(cp2k_default_input["geo_opt"])
update_d = {
    "GLOBAL": {"RUN_TYPE": "CELL_OPT"},
    "MOTION": {
        "PRINT": {
            "TRAJECTORY": {"EACH": {"CELL_OPT": 1}},
            "VELOCITIES": {"EACH": {"CELL_OPT": 1}},
        },
        "CELL_OPT": {
            "OPTIMIZER": "BFGS",
            "KEEP_ANGLES": ".TRUE.",
            "MAX_ITER": 200,
            "MAX_FORCE": 4.5e-5,
        },
    },
}
update_dict(cp2k_default_input["cell_opt"], update_d)
# <<<<<<<<<<<<<<<<<<<< CELL_OPT <<<<<<<<<<<<<<<<<<<<

# >>>>>>>>>>>>>>>>>>>> BOMD >>>>>>>>>>>>>>>>>>>>
cp2k_default_input["bomd"] = copy.deepcopy(cp2k_default_input["energy"])

# default thermostat is Nose-Hoover
update_d = {
    "GLOBAL": {"RUN_TYPE": "MD"},
    "MOTION": {
        "MD": {
            "ENSEMBLE": "NVT",
            "TIMESTEP": 0.5,
            "STEPS": 2000000,
            "TEMPERATURE": 330,
            "THERMOSTAT": {
                "NOSE": {},
            },
        },
        "PRINT": {
            "TRAJECTORY": {"EACH": {"MD": 1}},
            "VELOCITIES": {"EACH": {"MD": 1}},
            "FORCES": {"_": "ON"},
            "RESTART_HISTORY": {"EACH": {"MD": 1000}},
            "RESTART": {"BACKUP_COPIES": 3},
        },
    },
}

update_dict(cp2k_default_input["bomd"], update_d)
# <<<<<<<<<<<<<<<<<<<< BOMD <<<<<<<<<<<<<<<<<<<<

# >>>>>>>>>>>>>>>>>>>> SG-CPMD >>>>>>>>>>>>>>>>>>>>

cp2k_default_input["sgcpmd"] = copy.deepcopy(cp2k_default_input["bomd"])
cp2k_default_input["sgcpmd"]["MOTION"]["MD"].pop("THERMOSTAT")
# turn off smearing and diag in sgcpmd (incompatible with OT)
cp2k_default_input["sgcpmd"]["FORCE_EVAL"]["DFT"]["SCF"].pop("SMEAR")
cp2k_default_input["sgcpmd"]["FORCE_EVAL"]["DFT"]["SCF"].pop("DIAGONALIZATION")

update_d = {
    "FORCE_EVAL": {
        "DFT": {
            "QS": {
                "EPS_DEFAULT": 1.0e-13,
                "EXTRAPOLATION": "ASPC",
                "EXTRAPOLATION_ORDER": 0,
            },
            "SCF": {
                "SCF_GUESS": "RESTART",
                "EPS_SCF": 3.0e-7,
                "MAX_SCF": 50,
                "MAX_SCF_HISTORY": 5,
                "CHOLESKY": "INVERSE_DBCSR",
                "OUTER_SCF": {"EPS_SCF": 1.0e-6, "MAX_SCF": 20},
                "OT": {
                    "MINIMIZER": "DIIS",
                    "PRECOND_SOLVER": "INVERSE_UPDATE",
                    "PRECONDITIONER": "FULL_SINGLE_INVERSE",
                    "STEPSIZE": 0.01,
                },
                "PRINT": {"RESTART": {"EACH": {"MD": 20}}},
            },
        }
    },
    "MOTION": {
        "MD": {
            "ENSEMBLE": "LANGEVIN",
            "LANGEVIN": {
                "GAMMA": 0.001,
                "NOISY_GAMMA": 0,
            },
            "THERMAL_REGION": {
                "DO_LANGEVIN_DEFAULT": ".TRUE.",
                "FORCE_RESCALING": ".TRUE.",
                "PRINT": {
                    "TEMPERATURE": {},
                    "LANGEVIN_REGIONS": {
                        "FILENAME": "__STD_OUT__",
                    },
                },
                "DEFINE_REGION": [
                    {
                        "TEMPERATURE": 330.0,
                        "NOISY_GAMMA_REGION": 2.2e-4,
                        "LIST": "",
                    },
                    {
                        "TEMPERATURE": 330.0,
                        "NOISY_GAMMA_REGION": 5.0e-5,
                        "LIST": "",
                    },
                ],
            },
            "TEMP_KIND": ".TRUE.",
        }
    },
}

update_dict(cp2k_default_input["sgcpmd"], update_d)
# <<<<<<<<<<<<<<<<<<<< SG-CPMD <<<<<<<<<<<<<<<<<<<<

cp2k_restart_pbc = {
    "FORCE_EVAL": {
        "DFT": {
            "SCF": {
                "MAX_SCF": 3,
                "EPS_DIIS": 1e-15,
                "DIAGONALIZATION": {"_": ".FALSE."},
                "MIXING": {"ALPHA": 0.0, "METHOD": "DIRECT_P_MIXING"},
            }
        }
    }
}

# >>>>>>>>>>>>>>>>>>>> QMMM >>>>>>>>>>>>>>>>>>>>
cp2k_default_input.update(
    {
        "qmmm": {
            "FORCE_EVAL": {
                "METHOD": "QMMM",
                "DFT": {
                    "BASIS_SET_FILE_NAME": [
                        "BASIS_MOLOPT",
                        "BASIS_ADMM",
                        "BASIS_ADMM_MOLOPT",
                        "BASIS_MOLOPT-HSE06",
                    ],
                    "POTENTIAL_FILE_NAME": "GTH_POTENTIALS",
                    "MGRID": {"CUTOFF": 400, "COMMENSURATE": ".TRUE."},
                    "QS": {"EPS_DEFAULT": 1.0e-13},
                    "SCF": {
                        "SCF_GUESS": "RESTART",
                        "EPS_SCF": 1.0e-6,
                        "MAX_SCF": 500,
                        "ADDED_MOS": 500,
                        "CHOLESKY": "INVERSE",
                        "SMEAR": {
                            "_": "ON",
                            "METHOD": "FERMI_DIRAC",
                            "ELECTRONIC_TEMPERATURE": 300,
                        },
                        "DIAGONALIZATION": {"ALGORITHM": "STANDARD"},
                        "MIXING": {
                            "METHOD": "BROYDEN_MIXING",
                            "ALPHA": 0.3,
                            "BETA": 1.5,
                            "NBUFFER": 8,
                        },
                    },
                    "XC": {
                        "XC_FUNCTIONAL": {"_": "PBE"},
                        "vdW_POTENTIAL": {
                            "DISPERSION_FUNCTIONAL": "PAIR_POTENTIAL",
                            "PAIR_POTENTIAL": {
                                "TYPE": "DFTD3",
                                "PARAMETER_FILE_NAME": "dftd3.dat",
                                "REFERENCE_FUNCTIONAL": "PBE",
                            },
                        },
                    },
                },
                "MM": {
                    "FORCEFIELD": {
                        "IGNORE_MISSING_CRITICAL_PARAMS": ".TRUE.",
                        "DO_NONBONDED": ".FALSE.",
                        "CHARGE": [
                            {"ATOM": "O", "CHARGE": -0.8476},
                            {"ATOM": "H", "CHARGE": 0.4238},
                            {"ATOM": "Pt", "CHARGE": 0.0},
                        ],
                        "SPLINE": {"EMAX_SPLINE": 1e-1},
                    },
                    "POISSON": {
                        "EWALD": {"EWALD_TYPE": "SPME", "ALPHA": 0.4, "GMAX": 80}
                    },
                },
                "QMMM": {
                    "CENTER": "NEVER",
                    "ECOUPL": "GAUSS",
                    "USE_GEEP_LIB": 9,
                    "QM_KIND": [],
                },
                "SUBSYS": {
                    "TOPOLOGY": {
                        "COORD_FILE_FORMAT": "XYZ",
                        "COORD_FILE_NAME": "coord.xyz",
                    },
                    "KIND": [
                        {
                            "_": "O",
                            "POTENTIAL": "GTH-PBE-q6",
                            "BASIS_SET": "DZVP-MOLOPT-SR-GTH",
                        },
                        {
                            "_": "H",
                            "POTENTIAL": "GTH-PBE-q1",
                            "BASIS_SET": "DZVP-MOLOPT-SR-GTH",
                        },
                        {
                            "_": "Pt",
                            "POTENTIAL": "GTH-PBE-q10",
                            "BASIS_SET": "DZVP-A5-Q10-323-MOL-T1-DERIVED_SET-1",
                        },
                    ],
                },
                "PRINT": {"FORCES": {"_": "ON"}},
            },
            "GLOBAL": {"PROJECT": "cp2k"},
        }
    }
)
# <<<<<<<<<<<<<<<<<<<< QMMM <<<<<<<<<<<<<<<<<<<<

# >>>>>>>>>>>>>>>>>>>> DPMD >>>>>>>>>>>>>>>>>>>>
cp2k_default_input["dpmd"] = {
    "FORCE_EVAL": {
        "METHOD": "FIST",
        "MM": {
            "FORCEFIELD": {
                "IGNORE_MISSING_CRITICAL_PARAMS": ".TRUE.",
                "CHARGE": [],
                "NONBONDED": {
                    "DEEPMD": [],
                },
            },
            "POISSON": {"EWALD": {"EWALD_TYPE": "NONE"}},
        },
        "SUBSYS": {
            "TOPOLOGY": {"COORD_FILE_FORMAT": "XYZ", "COORD_FILE_NAME": "coord.xyz"},
        },
    },
    "GLOBAL": {"PROJECT": "cp2k", "RUN_TYPE": "MD"},
    "MOTION": {
        "MD": {
            "ENSEMBLE": "NVT",
            "TIMESTEP": 0.5,
            "STEPS": 2000000,
            "TEMPERATURE": 330,
            "THERMOSTAT": {
                "NOSE": {},
            },
        },
        "PRINT": {
            "TRAJECTORY": {"EACH": {"MD": 100}},
            "VELOCITIES": {"EACH": {"MD": 100}},
            "FORCES": {"_": "ON"},
            "RESTART_HISTORY": {"EACH": {"MD": 50000}},
            "RESTART": {"BACKUP_COPIES": 3},
        },
    },
}
# <<<<<<<<<<<<<<<<<<<< DPMD <<<<<<<<<<<<<<<<<<<<


####################################### setup in the following has not been checked #######################################
default_ffmd = {
    "FORCE_EVAL": {
        "METHOD": "FIST",
        "MM": {
            "FORCEFIELD": {
                "IGNORE_MISSING_CRITICAL_PARAMS": ".TRUE.",
                "CHARGE": [],
                "SPLINE": {"EMAX_SPLINE": 50},
            },
            "POISSON": {"EWALD": {"EWALD_TYPE": "ewald", "ALPHA": 0.5, "GMAX": 21}},
        },
        "SUBSYS": {
            "TOPOLOGY": {"COORD_FILE_FORMAT": "XYZ", "COORD_FILE_NAME": "coord.xyz"},
        },
    },
    "GLOBAL": {"PROJECT": "cp2k", "RUN_TYPE": "MD"},
    "MOTION": {
        "MD": {
            "ENSEMBLE": "NVT",
            "TIMESTEP": 0.5,
            "STEPS": 2000000,
            "TEMPERATURE": 330,
            "THERMOSTAT": {
                "REGION": "MOLECULE",
                "NOSE": {},
            },
        },
        "PRINT": {
            "TRAJECTORY": {"EACH": {"MD": 100}},
            "VELOCITIES": {"EACH": {"MD": 100}},
            "FORCES": {"_": "ON"},
            "RESTART_HISTORY": {"EACH": {"MD": 50000}},
            "RESTART": {"BACKUP_COPIES": 3},
        },
    },
}

# TODO: SPC/E water template
spce_wat = {
    "FORCE_EVAL": {
        "METHOD": "FIST",
        "MM": {
            "FORCEFIELD": {
                "IGNORE_MISSING_CRITICAL_PARAMS": ".TRUE.",
                "BEND": [{"ATOMS": "H O H", "K": 1.0, "THETA0": 1.91}],
                "BOND": [{"ATOMS": "O H", "K": 1.0, "R0": 1.89}],
                "CHARGE": [
                    {"ATOM": "O", "CHARGE": -0.8476},
                    {"ATOM": "H", "CHARGE": 0.4238},
                ],
                "NONBONDED": {"LENNARD-JONES": []},
                "SPLINE": {"EMAX_SPLINE": 50},
            },
            "POISSON": {"EWALD": {"EWALD_TYPE": "ewald", "ALPHA": 0.5, "GMAX": 21}},
        },
        "SUBSYS": {
            "TOPOLOGY": {
                "COORD_FILE_FORMAT": "XYZ",
                "COORD_FILE_NAME": "coord.xyz",
                "DUMP_PSF": {},
                "GENERATE": {"CREATE_MOLECULES": "T", "REORDER": "T"},
            },
        },
    },
    "GLOBAL": {"PROJECT": "cp2k", "RUN_TYPE": "MD"},
    "MOTION": {
        "MD": {
            "ENSEMBLE": "NVT",
            "TIMESTEP": 0.5,
            "STEPS": 2000000,
            "TEMPERATURE": 330,
            "THERMOSTAT": {
                "REGION": "MOLECULE",
                "NOSE": {},
            },
        },
        "PRINT": {
            "TRAJECTORY": {"EACH": {"MD": 100}},
            "VELOCITIES": {"EACH": {"MD": 100}},
            "FORCES": {"_": "ON"},
            "RESTART_HISTORY": {"EACH": {"MD": 50000}},
            "RESTART": {"BACKUP_COPIES": 3},
        },
    },
}
"""
**Lennard-Jones potential**

In CP2K/lammps V(r) = 4.0 * EPSILON * [(SIGMA/r)^12-(SIGMA/r)^6]
In lj_param dict, epsilon in kcal/mol and simga in angstrom

Parameters from:
- https://docs.lammps.org/Howto_spc.html# (SPC/E water)
- https://pubs.acs.org/doi/10.1021/jp801931d (metal)
- https://pubs.acs.org/doi/full/10.1021/ct500918t (ion)

Other sources:
- https://pubs.acs.org/doi/10.1021/acs.jctc.9b00941 (Na, K, Cl)
- https://pubs.acs.org/doi/10.1021/ct600252r
- https://ambermd.org/AmberModels_ions.php
"""
lj_params = {
    "O": {"epsilon": 0.1553, "sigma": 3.166},
    "H": {"epsilon": 0.0, "sigma": 0.4},
    "Ag": {"epsilon": 4.56, "sigma": 2.6326057121047},
    "Al": {"epsilon": 4.02, "sigma": 2.60587875056049},
    "Au": {"epsilon": 5.29, "sigma": 2.62904211723214},
    "Cu": {"epsilon": 4.72, "sigma": 2.33059104665513},
    "Ni": {"epsilon": 5.65, "sigma": 2.27357352869415},
    "Pb": {"epsilon": 2.93, "sigma": 3.17605393017031},
    "Pd": {"epsilon": 6.15, "sigma": 2.51144348643762},
    "Pt": {"epsilon": 7.80, "sigma": 2.53460685310927},
    "Li": {"epsilon": 0.00274091, "sigma": 2.2415011748410936},
    "Na": {"epsilon": 0.02639002, "sigma": 2.5907334723521065},
    "K": {"epsilon": 0.12693448, "sigma": 2.998765085260382},
    "Rb": {"epsilon": 0.20665151, "sigma": 3.192981005814976},
    "Cs": {"epsilon": 0.34673208, "sigma": 3.4798503930561653},
    "F": {"epsilon": 0.22878796, "sigma": 3.2410895365945542},
    "Cl": {"epsilon": 0.64367011, "sigma": 4.112388482935805},
    "Br": {"epsilon": 0.74435812, "sigma": 4.401039667613277},
    "I": {"epsilon": 0.86877007, "sigma": 4.9355788984974795},
}
# lj_param = {
#     "O": {
#         "epsilon": 0.0067,
#         "sigma": 3.166
#     },
#     "H": {
#         "epsilon": 0.0,
#         "sigma": 3.305
#     },
#     "Ag": {
#         "epsilon": 0.19608,
#         "sigma": 2.6326057121047
#     },
#     "Al": {
#         "epsilon": 0.17286,
#         "sigma": 2.60587875056049
#     },
#     "Au": {
#         "epsilon": 0.22747,
#         "sigma": 2.62904211723214
#     },
#     "Cu": {
#         "epsilon": 0.20296,
#         "sigma": 2.33059104665513
#     },
#     "Ni": {
#         "epsilon": 0.24295,
#         "sigma": 2.27357352869415
#     },
#     "Pb": {
#         "epsilon": 0.12599,
#         "sigma": 3.17605393017031
#     },
#     "Pd": {
#         "epsilon": 0.26445,
#         "sigma": 2.51144348643762
#     },
#     "Pt": {
#         "epsilon": 0.26445,
#         "sigma": 2.53460685310927
#     },
#     "Li": {
#         "epsilon": 0.00011884929615156104,
#         "sigma": 2.2415011748410936
#     },
#     "Na": {
#         "epsilon": 0.0011443043742500188,
#         "sigma": 2.5907334723521065
#     },
#     "K": {
#         "epsilon": 0.0055040382958084725,
#         "sigma": 2.998765085260382
#     },
#     "Rb": {
#         "epsilon": 0.008960668723948352,
#         "sigma": 3.192981005814976
#     },
#     "Cs": {
#         "epsilon": 0.015034737974310267,
#         "sigma": 3.4798503930561653
#     },
#     "F": {
#         "epsilon": 0.009920532966770708,
#         "sigma": 3.2410895365945542
#     },
#     "Cl": {
#         "epsilon": 0.027910343472532066,
#         "sigma": 4.112388482935805
#     },
#     "Br": {
#         "epsilon": 0.03227630190217819,
#         "sigma": 4.401039667613277
#     },
#     "I": {
#         "epsilon": 0.03767096013259918,
#         "sigma": 4.9355788984974795
#     }
# }
