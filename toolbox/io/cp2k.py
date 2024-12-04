# SPDX-License-Identifier: LGPL-3.0-or-later
import copy
import datetime
import re

from ase import Atoms
from ase.io.cube import read_cube_data
from cp2kdata import Cp2kPdos as _Cp2kPdos
from cp2kdata.block_parser.cells import parse_all_cells
from cp2kdata.block_parser.coordinates import parse_init_atomic_coordinates
from cp2kdata.block_parser.energies import parse_energies_list
from cp2kdata.block_parser.forces import parse_atomic_forces_list
from scipy import constants

from .. import CONFIGS
from ..utils import *
from ..utils.unit import *
from ..utils.utils import iterdict, save_dict_json, update_dict
from .template import cp2k_default_input

# from cp2kdata.pdos.pdos import gaussian_filter1d


_EPSILON = constants.epsilon_0 / constants.elementary_charge * constants.angstrom


class Cp2kInput:
    """
    Class for CP2K input file generation (on the basis of templates)

    Attributes
    ----------
    atoms: ASE Atoms object
        TBC
    input_type: str
        TBC
    pp_dir: str
        directory for basis set, peusudopotential, etc.
    wfn_restart: str
        wfn file for restart, see ref:
    qm_charge: float
        charge in QS
    multiplicity: int
        ref:
    uks: boolen
        ref:
    cutoff: int
        ref:
    rel_cutoff: int
        ref:

    Examples
    --------
    >>> from ase import io
    >>> from toolbox.io.cp2k import Cp2kInput
    >>> atoms = io.read("POSCAR")
    >>> input = Cp2kInput(atoms,
    >>>                   pp_dir="/data/basis",
    >>>                   hartree=True,
    >>>                   eden=True)
    >>> input.write()
    """

    def __init__(self, atoms, input_type="energy", **kwargs) -> None:
        self.atoms = atoms
        self.input_dict = copy.deepcopy(cp2k_default_input[input_type])
        # print(kwargs)
        # read user setup in config file
        try:
            update_d = CONFIGS["io"]["cp2k"]["input"]
        except:
            update_d = {}
        update_dict(kwargs, update_d)
        self.set_params(kwargs)

    def set_params(self, kwargs):
        for kw, value in kwargs.items():
            update_d = getattr(self, "set_%s" % kw)(value)
            update_dict(self.input_dict, update_d)

    def write(self, output_dir=".", fp_params={}, save_dict=False):
        """
        generate coord.xyz and input.inp for CP2K calculation at output_dir

        Parameters
        ----------
        output_dir : str
            directory to store coord.xyz and input.inp
        fp_params : dict
            dict for updated parameters
        """
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        cell = self.atoms.get_cell()
        cell_a = np.array2string(
            cell[0], formatter={"float_kind": lambda x: "%.4f" % x}
        )
        cell_a = cell_a[1:-1]
        cell_b = np.array2string(
            cell[1], formatter={"float_kind": lambda x: "%.4f" % x}
        )
        cell_b = cell_b[1:-1]
        cell_c = np.array2string(
            cell[2], formatter={"float_kind": lambda x: "%.4f" % x}
        )
        cell_c = cell_c[1:-1]

        user_config = fp_params
        update_dict(self.input_dict, user_config)

        if self.input_dict["FORCE_EVAL"].get("QMMM", None) is not None:
            cell_config = {
                "FORCE_EVAL": {
                    "SUBSYS": {"CELL": {"A": cell_a, "B": cell_b, "C": cell_c}},
                    "QMMM": {
                        "CELL": {
                            "A": cell_a,
                            "B": cell_b,
                            "C": cell_c,
                            "PERIODIC": "XYZ",
                        }
                    },
                }
            }
        else:
            cell_config = {
                "FORCE_EVAL": {
                    "SUBSYS": {"CELL": {"A": cell_a, "B": cell_b, "C": cell_c}}
                }
            }

        update_dict(self.input_dict, cell_config)
        # output list
        input_str = iterdict(self.input_dict, out_list=["\n"], loop_idx=0)
        # del input_str[0]
        # del input_str[-1]
        # print(input_str)
        str = "\n".join(input_str)
        str = str.strip("\n")

        io.write(os.path.join(output_dir, "coord.xyz"), self.atoms)
        with open(os.path.join(output_dir, "input.inp"), "w", encoding="utf-8") as f:
            f.write(str)

        if save_dict:
            save_dict_json(self.input_dict, os.path.join(output_dir, "input.json"))

    def set_project(self, project_name: str):
        update_d = {"GLOBAL": {"PROJECT": project_name}}
        return update_d

    def set_pp_dir(self, pp_dir):
        pp_dir = os.path.abspath(pp_dir)
        update_d = {
            "FORCE_EVAL": {
                "DFT": {
                    "BASIS_SET_FILE_NAME": [
                        os.path.join(pp_dir, "BASIS_MOLOPT"),
                        os.path.join(pp_dir, "BASIS_ADMM"),
                        os.path.join(pp_dir, "BASIS_ADMM_MOLOPT"),
                        os.path.join(pp_dir, "BASIS_MOLOPT-HSE06"),
                    ],
                    "POTENTIAL_FILE_NAME": os.path.join(pp_dir, "GTH_POTENTIALS"),
                    "XC": {
                        "vdW_POTENTIAL": {
                            "PAIR_POTENTIAL": {
                                "PARAMETER_FILE_NAME": os.path.join(pp_dir, "dftd3.dat")
                            }
                        }
                    },
                }
            }
        }
        return update_d

    def set_wfn_restart(self, wfn_file):
        update_d = {}
        if wfn_file is not None:
            update_d = {
                "FORCE_EVAL": {
                    "DFT": {"WFN_RESTART_FILE_NAME": os.path.abspath(wfn_file)}
                }
            }
        return update_d

    def set_qm_charge(self, charge):
        update_d = {"FORCE_EVAL": {"DFT": {"CHARGE": charge}}}
        return update_d

    def set_multiplicity(self, multiplicity):
        update_d = {"FORCE_EVAL": {"DFT": {"MULTIPLICITY": multiplicity}}}
        return update_d

    def set_uks(self, flag):
        if flag:
            update_d = {"FORCE_EVAL": {"DFT": {"UKS": ".TRUE."}}}
            return update_d
        else:
            return {}

    def set_cutoff(self, cutoff):
        update_d = {"FORCE_EVAL": {"DFT": {"MGRID": {"CUTOFF": cutoff}}}}
        return update_d

    def set_rel_cutoff(self, rel_cutoff):
        update_d = {"FORCE_EVAL": {"DFT": {"MGRID": {"REL_CUTOFF": rel_cutoff}}}}
        return update_d

    def set_kp(self, kp_mp):
        update_d = {
            "FORCE_EVAL": {
                "DFT": {
                    "KPOINTS": {
                        "SCHEME MONKHORST-PACK": "%d %d %d"
                        % (kp_mp[0], kp_mp[1], kp_mp[2]),
                        "SYMMETRY": ".TRUE.",
                        "EPS_GEO": 1.0e-8,
                        "FULL_GRID": ".TRUE.",
                        "PARALLEL_GROUP_SIZE": 0,
                    }
                }
            }
        }
        return update_d

    def set_max_scf(self, max_scf: int):
        update_d = {"FORCE_EVAL": {"DFT": {"SCF": {"MAX_SCF": max_scf}}}}
        return update_d

    def set_eps_scf(self, eps_scf: float):
        update_d = {"FORCE_EVAL": {"DFT": {"SCF": {"EPS_SCF": eps_scf}}}}
        return update_d

    def set_dip_cor(self, flag):
        if flag:
            update_d = {"FORCE_EVAL": {"DFT": {"SURFACE_DIPOLE_CORRECTION": ".TRUE."}}}
            return update_d
        else:
            return {}

    def set_eden(self, flag):
        if flag:
            update_d = {
                "FORCE_EVAL": {
                    "DFT": {
                        "PRINT": {
                            "E_DENSITY_CUBE": {"ADD_LAST": "NUMERIC", "STRIDE": "8 8 1"}
                        }
                    }
                }
            }
            return update_d
        else:
            return {}

    def set_mo(self, flag):
        if flag:
            update_d = {
                "FORCE_EVAL": {"DFT": {"PRINT": {"MO_CUBES": {"ADD_LAST": "NUMERIC"}}}}
            }
            return update_d
        else:
            return {}

    def set_pdos(self, flag):
        if flag:
            update_d = {
                "FORCE_EVAL": {
                    "DFT": {
                        "PRINT": {
                            "PDOS": {
                                "COMPONENTS": ".TRUE.",
                                "ADD_LAST": "NUMERIC",
                                "NLUMO": -1,
                                "COMMON_ITERATION_LEVELS": 0,
                            }
                        }
                    }
                }
            }
            return update_d
        else:
            return {}

    def set_hartree(self, flag):
        if flag:
            update_d = {
                "FORCE_EVAL": {
                    "DFT": {
                        "PRINT": {
                            "V_HARTREE_CUBE": {"ADD_LAST": "NUMERIC", "STRIDE": "8 8 1"}
                        }
                    }
                }
            }
            return update_d
        else:
            return {}

    def set_efield(self, flag):
        if flag:
            update_d = {
                "FORCE_EVAL": {
                    "DFT": {
                        "PRINT": {
                            "EFIELD_CUBE": {"ADD_LAST": "NUMERIC", "STRIDE": "8 8 1"}
                        }
                    }
                }
            }
            return update_d
        else:
            return {}

    def set_totden(self, flag):
        if flag:
            update_d = {
                "FORCE_EVAL": {
                    "DFT": {
                        "PRINT": {
                            "TOT_DENSITY_CUBE": {
                                "ADD_LAST": "NUMERIC",
                                "STRIDE": "1 1 1",
                            }
                        }
                    }
                }
            }
            return update_d
        else:
            return {}

    def set_extended_fft_lengths(self, flag):
        if flag:
            update_d = {"GLOBAL": {"EXTENDED_FFT_LENGTHS": ".TRUE."}}
            return update_d
        else:
            return {}

    def set_smear(self, flag):
        if flag == False:
            update_d = {
                "FORCE_EVAL": {
                    "DFT": {
                        "SCF": {
                            "ADDED_MOS": 0,
                            "CHOLESKY": "RESTORE",
                            "SMEAR": {"_": ".FALSE."},
                            "DIAGONALIZATION": {"_": ".FALSE."},
                        }
                    }
                }
            }
            return update_d

    def set_mlwf(self, flag):
        if flag:
            update_d = {
                "FORCE_EVAL": {
                    "DFT": {
                        "LOCALIZE": {
                            "METHOD": "CRAZY",
                            "EPS_LOCALIZATION": 1e-08,
                            "PRINT": {"WANNIER_CENTERS": {"IONS+CENTERS": ".TRUE."}},
                        }
                    }
                }
            }
            return update_d
        else:
            return {}

    def set_kind(self, kind_dict: dict):
        """
        Parameters
        ----------
        kind_dict : dict
            dict to update kind section, for example:
            {
                "S": {
                    "ELEMENT": "O",
                    "BASIS_SET": "DZVP-MOLOPT-SR-GTH",
                    "POTENTIAL": "GTH-PBE-q6"
                },
                "Li": {
                    "ELEMENT": "H",
                    "BASIS_SET": "DZVP-MOLOPT-SR-GTH",
                    "POTENTIAL": "GTH-PBE-q1"
                }
            }
        """

        update_d = {"FORCE_EVAL": {"SUBSYS": {}}}
        update_dict(self.input_dict, update_d)

        old_kind_list = self.input_dict["FORCE_EVAL"]["SUBSYS"].get("KIND", [])
        if len(old_kind_list) > 0:
            for k, v in kind_dict.items():
                tmp_dict = copy.deepcopy(v)
                tmp_dict.update({"_": k})
                flag = False
                for ii, item in enumerate(old_kind_list):
                    if k == item["_"]:
                        # print(v)
                        old_kind_list[ii] = tmp_dict
                        flag = True
                        break
                if flag == False:
                    old_kind_list.append(tmp_dict)
        return {}

    def set_restart(self, flag):
        if flag:
            update_d = {
                "EXT_RESTART": {
                    "RESTART_FILE_NAME": "%s-1.restart"
                    % self.input_dict["GLOBAL"]["PROJECT"]
                }
            }
            return update_d
        else:
            return {}

    def set_md_step(self, md_step):
        update_d = {"MOTION": {"MD": {"STEPS": md_step}}}
        return update_d

    def set_md_temp(self, md_temp):
        update_d = {"MOTION": {"MD": {"TEMPERATURE": md_temp}}}
        return update_d

    def set_md_timestep(self, md_timestep):
        update_d = {"MOTION": {"MD": {"TIMESTEP": md_timestep}}}
        return update_d

    def update_dp(self, dp_model: str):
        from deepmd.infer import DeepPot

        dp = DeepPot(dp_model)
        type_map = dp.tmap
        for ii, atype in enumerate(type_map):
            self.input_dict["FORCE_EVAL"]["MM"]["FORCEFIELD"]["CHARGE"].append(
                {
                    "ATOM": atype,
                    "CHARGE": 0.0,
                }
            )
            self.input_dict["FORCE_EVAL"]["MM"]["FORCEFIELD"]["NONBONDED"][
                "DEEPMD"
            ].append(
                {
                    "ATOMS": "%s %s" % (atype, atype),
                    "POT_FILE_NAME": dp_model,
                    "ATOM_DEEPMD_TYPE": ii,
                }
            )


class Cp2kOutput:
    def __init__(self, fname="output.out", ignore_warning=False) -> None:
        self.output_file = fname
        with open(fname, "r") as f:
            self.content = f.readlines()
        self.string = "".join(self.content)

        self.check_scf = not ignore_warning
        if (self.check_scf) and (self.scf_loop == -1):
            raise Warning("SCF run NOT converged")

        self.natoms = len(self.atoms)

    @property
    def worktime(self):
        pattern = r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}\.\d{3}"
        out = re.findall(pattern, self.string)
        t = re.split(r"\.|:|-|\s", out[0])
        start_time = datetime.datetime(
            int(t[0]), int(t[1]), int(t[2]), int(t[3]), int(t[4]), int(t[5]), int(t[6])
        )
        t = re.split(r"\.|:|-|\s", out[-1])
        end_time = datetime.datetime(
            int(t[0]), int(t[1]), int(t[2]), int(t[3]), int(t[4]), int(t[5]), int(t[6])
        )
        delta_t = end_time - start_time
        run_time = delta_t.total_seconds()
        return run_time

    def grep_text_search(self, pattern):
        search_pattern = re.compile(pattern)
        scf_pattern = re.compile(r"SCF run converged in")

        flag = False
        scf_flag = not self.check_scf
        for line in self.content:
            line = line.strip("\n")
            if scf_pattern.search(line) is not None:
                scf_flag = True
            if scf_flag == False:
                continue
            if search_pattern.search(line) is not None:
                flag = True
                break
        if flag:
            return line
        else:
            return ""

    def grep_texts(self, start_pattern, end_pattern):
        start_pattern = re.compile(start_pattern)
        end_pattern = re.compile(end_pattern)
        scf_pattern = re.compile(r"SCF run converged in")

        flag = False
        scf_flag = not self.check_scf
        data_lines = []
        nframe = 0
        for line in self.content:
            line = line.strip("\n")
            if scf_pattern.search(line) is not None:
                scf_flag = True
            if scf_flag == False:
                continue
            if start_pattern.match(line):
                flag = True
            if end_pattern.match(line):
                assert flag is True, (flag, "No data is found in this file.")
                flag = False
                nframe += 1
            if flag is True:
                data_lines.append(line)
        return nframe, np.reshape(data_lines, (nframe, -1))

    def grep_texts_by_nlines(self, start_pattern, nlines):
        start_pattern = re.compile(start_pattern)
        scf_pattern = re.compile(r"SCF run converged in")

        data_lines = []
        nframe = 0
        scf_flag = not self.check_scf
        for ii, line in enumerate(self.content):
            line = line.strip("\n")
            if scf_pattern.search(line) is not None:
                scf_flag = True
            if scf_flag == False:
                continue
            if start_pattern.search(line) is not None:
                data_lines.append(self.content[ii : ii + nlines])
                nframe += 1
                continue
        if nframe == 0:
            raise AttributeError("No data is found in this file.")
        return nframe, np.reshape(data_lines, (nframe, -1))

    @property
    def coord(self):
        """
        get atomic coordinate from cp2k output

        Return:
            coord numpy array (n_atom, 3)
        """
        out = parse_init_atomic_coordinates(self.string)
        self.chemical_symbols = out[2]
        return out[0]

    @property
    def atoms(self):
        positions = self.coord
        atoms = Atoms(symbols=self.chemical_symbols, positions=positions)
        check_scf = self.check_scf
        self.check_scf = False
        out = parse_all_cells(self.string)
        self.check_scf = check_scf

        atoms.set_cell(out[0])
        atoms.set_pbc(True)
        return atoms

    @property
    def force(self):
        """
        get atomic force from cp2k output

        Return:
            force numpy array (n_atom, 3)
        """
        out = parse_atomic_forces_list(self.string)
        return out[0] * AU_TO_EV_EVERY_ANG

    @property
    def energy(self):
        out = parse_energies_list(self.string)
        return out[0] * AU_TO_EV

    @property
    def scf_loop(self):
        pattern = r"\s+SCF\srun\sconverged\sin\s+\d+"
        out = re.findall(pattern, self.string)
        if len(out) == 0:
            return -1
        else:
            return int(out[0].split(" ")[-1])

    @property
    def fermi(self):
        pattern = r"\s+Fermi\senergy:\s+.\d\.\d+"
        out = re.findall(pattern, self.string)
        return float(out[0].split(" ")[-1]) * AU_TO_EV

    @property
    def m_charge(self):
        start_pattern = "Mulliken Population Analysis"
        nframe, data_lines = self.grep_texts_by_nlines(start_pattern, self.natoms + 3)
        data_list = []
        for line in data_lines[-1, 3:]:
            line_list = line.split()
            data_list.append(float(line_list[-1]))
        return np.reshape(data_list, -1)

    @property
    def h_charge(self):
        start_pattern = "Hirshfeld Charges"
        nframe, data_lines = self.grep_texts_by_nlines(start_pattern, self.natoms + 3)
        data_list = []
        for line in data_lines[-1, 3:]:
            line_list = line.split()
            data_list.append(float(line_list[-1]))
        return np.reshape(data_list, -1)

    @property
    def dipole_moment(self):
        pattern = "Dipole moment"
        nframe, data_lines = self.grep_texts_by_nlines(pattern, 2)

        data_list = []
        for line in data_lines[-1, 1:]:
            line_list = line.split()
            data_list.append(list(map(float, line_list[1:-1:2])))
        return np.reshape(data_list, 3)

    @property
    def surf_dipole_moment(self):
        """
        Total dipole moment perpendicular to
        the slab [electrons-Angstroem]:              -1.5878220042
        """
        pattern = "Total dipole moment perpendicular to"
        nframe, data_lines = self.grep_texts_by_nlines(pattern, 2)

        line = data_lines[-1, 1]
        line_list = line.split()
        return float(line_list[-1])

    @property
    def potdrop(self):
        cross_area = np.linalg.norm(np.cross(self.atoms.cell[0], self.atoms.cell[1]))
        DeltaV = self.surf_dipole_moment / cross_area / _EPSILON
        return DeltaV

    @property
    def energy_dict(self):
        energy_dict = {}

        start_pattern = r"  Total charge density g-space grids:"
        end_pattern = r"  Total energy:"
        nframe, data_lines = self.grep_texts(start_pattern, end_pattern)

        for kw in data_lines[-1, 2:-1]:
            kw = kw.split(":")
            k = kw[0].strip(" ")
            v = float(kw[1]) * AU_TO_EV
            energy_dict[k] = v

        energy_dict.pop("Fermi energy", None)

        tot_e = 0.0
        for v in energy_dict.values():
            tot_e += v
        energy_dict["total"] = tot_e

        return energy_dict


class MultiFrameCp2kOutput(Cp2kOutput):
    def __init__(self, fname="output.out", ignore_warning=False) -> None:
        super().__init__(fname, ignore_warning)

    @property
    def atoms(self):
        positions = self.coord
        atoms = Atoms(symbols=self.chemical_symbols, positions=positions)
        a = float(self.grep_text_search(r"Vector a").split()[-1])
        b = float(self.grep_text_search(r"Vector b").split()[-1])
        c = float(self.grep_text_search(r"Vector c").split()[-1])
        alpha = float(self.grep_text_search(r"Angle | alpha").split()[-1])
        beta = float(self.grep_text_search(r"Angle | beta").split()[-1])
        gamma = float(self.grep_text_search(r"Angle | gamma").split()[-1])
        atoms.set_cell([a, b, c, alpha, beta, gamma])
        atoms.set_pbc(True)
        return atoms

    @property
    def force(self):
        """
        get atomic force from cp2k output

        Return:
            force numpy array (n_atom, 3)
        """
        start_pattern = r" ATOMIC FORCES in"
        end_pattern = r" SUM OF ATOMIC FORCES"
        nframe, data_lines = self.grep_texts(start_pattern, end_pattern)
        data_lines = np.reshape(data_lines, (nframe, -1))

        data_list = []
        for line in data_lines[:, 3:].reshape(-1):
            line_list = line.split()
            data_list.append(
                [
                    float(line_list[3]) * AU_TO_EV_EVERY_ANG,
                    float(line_list[4]) * AU_TO_EV_EVERY_ANG,
                    float(line_list[5]) * AU_TO_EV_EVERY_ANG,
                ]
            )
        return np.reshape(data_list, (nframe, -1, 3))

    @property
    def energy(self):
        data = self.grep_text_search("Total energy: ")
        # data = self.grep_text_search("Total FORCE_EVAL")
        data = data.replace("\n", " ")
        data = data.split(" ")
        return float(data[-1]) * AU_TO_EV

    @property
    def scf_loop(self):
        pattern = r"\s+SCF\srun\sconverged\sin\s+\d+"
        out = re.findall(pattern, self.string)
        if len(out) == 0:
            return -1
        else:
            return int(out[0].split(" ")[-1])

    @property
    def fermi(self):
        pattern = r"\s+Fermi\sEnergy\s\[eV\]\s:\s+.\d\.\d+"
        out = re.findall(pattern, self.string)
        return float(out[0].split(" ")[-1]) * AU_TO_EV

    @property
    def m_charge(self):
        start_pattern = "Mulliken Population Analysis"
        nframe, data_lines = self.grep_texts_by_nlines(start_pattern, self.natoms + 3)
        data_lines = np.reshape(data_lines, (nframe, -1))

        data_list = []
        for line in data_lines[:, 3:].reshape(-1):
            line_list = line.split()
            data_list.append(float(line_list[-1]))
        return np.reshape(data_list, (nframe, -1))

    @property
    def h_charge(self):
        start_pattern = "Hirshfeld Charges"
        nframe, data_lines = self.grep_texts_by_nlines(start_pattern, self.natoms + 3)
        data_lines = np.reshape(data_lines, (nframe, -1))

        data_list = []
        for line in data_lines[:, 3:].reshape(-1):
            line_list = line.split()
            data_list.append(float(line_list[-1]))
        return np.reshape(data_list, (nframe, -1))

    @property
    def dipole_moment(self):
        pattern = "Dipole moment"
        nframe, data_lines = self.grep_texts_by_nlines(pattern, 2)
        data_lines = np.reshape(data_lines, (nframe, -1))

        data_list = []
        for line in data_lines[:, 1:].reshape(-1):
            line_list = line.split()
            data_list.append(list(map(float, line_list[1:-1:2])))
        return np.reshape(data_list, (nframe, 3))

    @property
    def surf_dipole_moment(self):
        """
        Total dipole moment perpendicular to
        the slab [electrons-Angstroem]:              -1.5878220042
        """
        pattern = "Total dipole moment perpendicular to"
        nframe, data_lines = self.grep_texts_by_nlines(pattern, 2)
        data_lines = np.reshape(data_lines, (nframe, -1))

        data_list = []
        for line in data_lines[:, 1:].reshape(-1):
            line_list = line.split()
            data_list.append(float(line_list[-1]))
        return np.reshape(data_list, (nframe))

    @property
    def potdrop(self):
        cross_area = np.linalg.norm(np.cross(self.atoms.cell[0], self.atoms.cell[1]))
        DeltaV = self.surf_dipole_moment / cross_area / _EPSILON
        return DeltaV

    @property
    def energy_dict(self):
        energy_dict = {}

        start_pattern = r"  Total charge density g-space grids:"
        end_pattern = r"  Total energy:"
        nframe, data_lines = self.grep_texts(start_pattern, end_pattern)
        data_lines = np.reshape(data_lines, (nframe, -1))

        tot_e = 0.0
        for kw in data_lines[-1, 2:-1].reshape(-1):
            kw = kw.split(":")
            k = kw[0].strip(" ")
            v = float(kw[1]) * AU_TO_EV
            energy_dict[k] = [v]
            tot_e += v
        # energy_dict["Total energy"].append(tot_e)
        # for kws in data_lines[1:, 2:-1].reshape(-1):
        #     tot_e = 0.
        #     for kw in kws:
        #         kw = kw.split(":")
        #         k = kw[0].strip(' ')
        #         v = float(kw[1]) * AU_TO_EV
        #         energy_dict[k].append(v)
        #         tot_e += v
        #     # energy_dict["Total energy"].append(tot_e)

        energy_dict.pop("Fermi energy", None)

        return energy_dict

    # @property
    # def xc_energy(self):
    #     data = self.grep_text_search("Exchange-correlation energy")
    #     data = data.replace('\n', ' ')
    #     data = data.split(' ')
    #     return float(data[-1]) * AU_TO_EV

    # @property
    # def core_hmt_energy(self):
    #     data = self.grep_text_search("Core Hamiltonian energy")
    #     data = data.replace('\n', ' ')
    #     data = data.split(' ')
    #     return float(data[-1]) * AU_TO_EV

    # @property
    # def overlap_energy(self):
    #     data = self.grep_text_search(
    #         "Overlap energy of the core charge distribution")
    #     data = data.replace('\n', ' ')
    #     data = data.split(' ')
    #     return float(data[-1]) * AU_TO_EV

    # @property
    # def self_energy(self):
    #     data = self.grep_text_search("Self energy of the core charge distribution")
    #     data = data.replace('\n', ' ')
    #     data = data.split(' ')
    #     return float(data[-1]) * AU_TO_EV

    # @property
    # def hartree_energy(self):
    #     data = self.grep_text_search("Hartree energy")
    #     data = data.replace('\n', ' ')
    #     data = data.split(' ')
    #     return float(data[-1]) * AU_TO_EV

    # @property
    # def vdw_energy(self):
    #     data = self.grep_text_search("Dispersion energy")
    #     data = data.replace('\n', ' ')
    #     data = data.split(' ')
    #     return float(data[-1]) * AU_TO_EV

    # @property
    # def entropy_e_energy(self):
    #     data = self.grep_text_search("Electronic entropic energy")
    #     data = data.replace('\n', ' ')
    #     data = data.split(' ')
    #     return float(data[-1]) * AU_TO_EV


class Cp2kCube:
    def __init__(self, fname) -> None:
        self.cube_data, self.atoms = read_cube_data(fname)
        self.n_grid = np.array(self.cube_data.shape)
        # 3 * 3
        self.cube_vectors = self.atoms.get_cell() / self.n_grid.reshape(3, 1)
        # cube volume [A^3]
        self.cube_volume = self.atoms.get_volume() / np.prod(self.n_grid)

        # cell_params = self.atoms.cell.cellpar()
        # try:
        #     assert not (False in (90. == cell_params[-3:]))
        #     self.cell_params = cell_params[:3]
        # except:
        #     raise ValueError("Cell is not orthogonal")

    # @property
    # def cube_grids(self):
    #     cube_grids = []
    #     for ii in range(3):
    #         cube_grids.append(
    #             np.arange(
    #                 0, self.cell_params[ii],
    #                 self.cell_params[ii] / self.cube_data.shape[ii]
    #             )[:self.cube_data.shape[ii]]
    #         )
    #     return cube_grids

    @property
    def mesh(self):
        """
        mesh for cube data
        """
        # generate mesh from self.cube_vectors
        x = np.arange(self.n_grid[0])
        y = np.arange(self.n_grid[1])
        z = np.arange(self.n_grid[2])
        coeff = np.meshgrid(x, y, z, indexing="ij")
        coeff_x = coeff[0].reshape(self.n_grid[0], self.n_grid[1], self.n_grid[2], 1)
        coeff_y = coeff[1].reshape(self.n_grid[0], self.n_grid[1], self.n_grid[2], 1)
        coeff_z = coeff[2].reshape(self.n_grid[0], self.n_grid[1], self.n_grid[2], 1)
        unit_x = self.cube_vectors[0].reshape(1, 1, 1, 3)
        unit_y = self.cube_vectors[1].reshape(1, 1, 1, 3)
        unit_z = self.cube_vectors[2].reshape(1, 1, 1, 3)
        # dimx * dimy * dimz * 3
        mesh = coeff_x * unit_x + coeff_y * unit_y + coeff_z * unit_z
        return mesh

    @property
    def dipole(self):
        """
        Dipole moment [e A]
        """
        # cube_data: charge density [e / bohr^3]
        charge = self.cube_data * (self.cube_volume / AU_TO_ANG**3)
        charge = np.reshape(charge, [self.n_grid[0], self.n_grid[1], self.n_grid[2], 1])
        dipole = np.sum(self.mesh * charge, axis=(0, 1, 2))
        # in cp2k, the electron density is positive
        return -dipole

    def get_ave_cube(self, axis=2, gaussian_sigma=0.0):
        if (
            hasattr(self, "axis")
            and self.axis == axis
            and hasattr(self, "ave_cube_data")
        ):
            pass
        else:
            self.axis = axis
            self.ave_grid = self.cube_vectors[self.axis][self.axis] * np.arange(
                self.n_grid[self.axis]
            )
            ave_axis = tuple(np.delete(np.arange(3), self.axis).tolist())
            self.ave_cube_data = np.mean(self.cube_data, axis=ave_axis)

        if gaussian_sigma > 0.0:
            self.ave_cube_data_convolve = gaussian_convolve(
                self.ave_grid, self.ave_cube_data, gaussian_sigma
            )
        else:
            self.ave_cube_data_convolve = copy.deepcopy(self.ave_cube_data)

        return (self.ave_grid, self.ave_cube_data, self.ave_cube_data_convolve)


class Cp2kHartreeCube(Cp2kCube):
    def __init__(
        self,
        fname,
        vac_region: list = None,
    ) -> None:
        super().__init__(fname)
        if vac_region:
            self.set_vac_region(vac_region)

    def set_vac_region(self, vac_region):
        assert len(vac_region) == 2
        self.vac_region = vac_region

    def get_ave_cube(self, axis=2, gaussian_sigma=0):
        if (
            hasattr(self, "axis")
            and self.axis == axis
            and hasattr(self, "ave_cube_data")
        ):
            pass
        else:
            # (self.ave_grid, self.ave_cube_data, self.ave_cube_data_convolve)
            output = super().get_ave_cube(axis, gaussian_sigma)
            self.ave_cube_data_convolve *= AU_TO_EV
            self.ave_cube_data *= AU_TO_EV
        return (self.ave_grid, self.ave_cube_data, self.ave_cube_data_convolve)

    @property
    def potdrop(self):
        start_id = np.argmin(np.abs(self.ave_grid - self.vac_region[0]))
        end_id = np.argmin(np.abs(self.ave_grid - self.vac_region[1]))
        if start_id > end_id:
            _data = np.append(
                self.ave_cube_data[start_id:], self.ave_cube_data[:end_id]
            )
        else:
            _data = np.array(self.ave_cube_data)[
                self.vac_region[0] : self.vac_region[1]
            ]
        dev_data = np.diff(_data, axis=0)
        p_jump = np.argmax(np.abs(dev_data))
        return dev_data[p_jump]

    @property
    def dipole(self):
        """
        Dipole moment of cell [e A]
        """
        d = -self.potdrop * self.cross_area * _EPSILON
        return d

    def set_cross_area(self, cross_area):
        self.cross_area = cross_area


class Cp2kPdos(_Cp2kPdos):
    def __init__(self, file_name, parse_file_name=True) -> None:
        super().__init__(file_name, parse_file_name)

    def get_raw_dos(self, dos_type="total", steplen=0.01):
        energies = self.energies
        fermi = self.fermi

        weights = self._get_raw_dos(dos_type)
        bins = int((energies[-1] - energies[0]) / steplen)
        dos, ener = np.histogram(
            energies, bins, weights=weights, range=(energies[0], energies[-1])
        )
        dos = dos / steplen
        ener = ener[:-1] - fermi + 0.5 * steplen
        self.dos = dos
        self.ener = ener
        return dos, ener

    def _get_raw_dos(self, dos_type):
        try:
            return getattr(self, "_get_raw_dos_%s" % dos_type)()
        except:
            raise NameError("dos type does not exist!")

    def _get_raw_dos_total(self):
        return np.loadtxt(self.file)[:, 3:].sum(axis=1)

    def _get_raw_dos_system(self):
        tmp_len = len(np.loadtxt(self.file, usecols=2))
        return np.ones(tmp_len)

    def _get_raw_dos_s(self):
        return np.loadtxt(self.file, usecols=3)

    def _get_raw_dos_p(self):
        return np.loadtxt(self.file, usecols=np.arange(4, 7)).sum(axis=1)

    def _get_raw_dos_d(self):
        return np.loadtxt(self.file, usecols=np.arange(7, 12)).sum(axis=1)

    def _get_raw_dos_f(self):
        return np.loadtxt(self.file, usecols=np.arange(12, 19)).sum(axis=1)

    def get_dos(self, sigma=0.2, dos_type="total", steplen=0.01):
        # smooth the dos data
        dos, ener = self.get_raw_dos(dos_type=dos_type, steplen=steplen)
        smth_dos = gaussian_filter1d(dos, sigma)
        self.smth_dos = smth_dos
        return smth_dos, ener

    @property
    def homo(self):
        homo_idx = np.where(self.occupation == 0)[0][0] - 1
        return self.energies[homo_idx] - self.fermi

    @property
    def lumo(self):
        return self.energies[self.occupation == 0][0] - self.fermi

    @property
    def vbm(self):
        raw_dos = self._get_raw_dos_total()
        mask = (self.occupation > 1e-5) & (raw_dos > 1e-3)
        try:
            return self.energies[mask].max() - self.fermi
        except:
            print("Warning: No VBM is found!")
            return self.energies.min() - self.fermi

    @property
    def cbm(self):
        raw_dos = self._get_raw_dos_total()
        mask = (self.occupation < 1e-5) & (raw_dos > 1e-3)
        try:
            return self.energies[mask].min() - self.fermi
        except:
            print("Warning: No CBM is found!")
            return self.energies.max() - self.fermi


"""
smoothing func: Gaussian kernel
"""


def gaussian_kernel(bins, sigma):
    if sigma == 0:
        output = np.zeros_like(bins)
        one_id = np.where(a == 0.0)[0][0]
        output[one_id] = 1
        return output
    elif sigma > 0:
        A = 1 / (sigma * np.sqrt(2 * np.pi))
        output = np.exp(-bins * bins / (2 * sigma**2))
        output *= A
        return output
    else:
        raise AttributeError("Sigma should be non-negative value.")


def gaussian_convolve(xs, ys, sigma):
    nbins = len(xs) - 1

    output = []
    for x in xs:
        bins = xs - x
        tmp_out = gaussian_kernel(bins, sigma)
        bin_width = bins[1:] - bins[:-1]
        output.append(
            np.sum(bin_width * ((tmp_out * ys)[1:] + (tmp_out * ys)[:-1]) / 2)
        )
    return np.array(output)


def generate_ids_list(ids):
    ids = np.sort(ids)
    diff = np.diff(ids)
    # ids_std = ""
    ids_list = []
    update_flag = True
    for count, ii in enumerate(diff):
        if update_flag:
            start_id = ids[count]
            update_flag = False
        end_id = ids[count]
        # print(start_id, end_id)
        if ii == 1:
            continue
        else:
            if start_id == end_id:
                # ids_std += "%d" % start_id
                ids_list.append(start_id)
            elif end_id > start_id:
                # ids_std += "%d..%d " % (start_id, end_id)
                ids_list.append("%d..%d" % (start_id, end_id))
            else:
                raise ValueError("end_id is smaller than start_id")
            update_flag = True

    if update_flag:
        # ids_std += "%d " % ids[-1]
        ids_list.append(ids[-1])
    else:
        # ids_std += "%d..%d " % (start_id, ids[-1])
        ids_list.append("%d..%d" % (start_id, ids[-1]))

    # return ids_std.strip()
    return ids_list


def read_wannier_spread(fname: str = "wannier_spread.out"):
    """
    Read wannier spread file generated by cp2k

    Parameters
    ----------
    fname : str
        wannier_spread.out file name

    Returns
    -------
    wannier_spread : numpy array
        wannier spread data
    """
    # skip 1, 2, and the last lines
    # save the others as array
    with open(fname, "r", encoding="UTF-8") as f:
        lines = f.readlines()[2:-1]

    wannier = []
    for line in lines:
        # line to array
        line = line.strip().split()
        wannier.append([float(line[1]), float(line[2])])
    return np.array(wannier)
