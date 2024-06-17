# SPDX-License-Identifier: LGPL-3.0-or-later
from typing import List

import dpdata
import numpy as np
from ase import io
from ase.geometry.cell import cell_to_cellpar
from dpdata.data_type import Axis, DataType
from MDAnalysis.lib.distances import distance_array, minimize_vectors

from .cp2k import Cp2kOutput

# import re


# from dpdata.unit import econvs
# from deepmd.infer.ewald_recp import EwaldRecp


dpdata.LabeledSystem.register_data_type(
    DataType("ext_efield", np.ndarray, (Axis.NFRAMES, 3), required=False),
    DataType("efield", np.ndarray, (Axis.NFRAMES, Axis.NATOMS, 3), required=False),
    DataType("aparam", np.ndarray, (Axis.NFRAMES, Axis.NATOMS, -1), required=False),
    DataType("atomic_dipole", np.ndarray, (Axis.NFRAMES, -1), required=False),
    DataType("charge", np.ndarray, (Axis.NFRAMES, Axis.NATOMS), required=False),
)

"""
todo:
- [ ] add_efield
- [ ] add_aparam
"""


class CP2KDPDataSystem:
    def __init__(
        self,
        fname: str,
    ) -> None:
        self.cp2k_out = Cp2kOutput(fname)
        self.dp_sys = dpdata.LabeledSystem(fname, fmt="cp2k/output")
        self.nframes = self.dp_sys.get_nframes()
        self.natoms = self.dp_sys.get_natoms()

    def add(self, kw, **kwargs):
        getattr(self, "add_%s" % kw)(**kwargs)

    def add_atomic_dipole(self, type_map, sel_type, wannier_fname):
        wannier_atoms = io.read(wannier_fname)
        sel_ids = self.get_sel_ids(type_map, sel_type)

        coords = self.dp_sys.data["coords"].reshape(-1, 3)
        cellpar = cell_to_cellpar(self.dp_sys.data["cells"].reshape(3, 3))
        extended_coords = coords.copy()

        ref_coords = coords[sel_ids].reshape(-1, 3)
        e_coords = wannier_atoms.get_positions()
        dist_mat = distance_array(ref_coords, e_coords, box=cellpar)
        atomic_dipole = []
        for ii, dist_vec in enumerate(dist_mat):
            mask = dist_vec < 1.0
            cn = np.sum(mask)
            assert cn == 4
            # print(cn)
            wc_coord_rel = e_coords[mask] - ref_coords[ii]
            wc_coord_rel = minimize_vectors(wc_coord_rel, box=cellpar)
            _atomic_dipole = wc_coord_rel.mean(axis=0)
            atomic_dipole.append(_atomic_dipole)
            wc_coord = _atomic_dipole + ref_coords[ii]
            extended_coords = np.concatenate(
                (extended_coords, wc_coord.reshape(1, 3)), axis=0
            )
        atomic_dipole = np.reshape(atomic_dipole, (-1, 3))
        assert atomic_dipole.shape[0] == len(sel_ids)
        self.dp_sys.data["atomic_dipole"] = atomic_dipole.reshape(self.nframes, -1)
        self.extended_coords = extended_coords

    def add_ext_efield(self, ext_efield):
        self.dp_sys.data["ext_efield"] = ext_efield.reshape(self.nframes, 3)

    def add_charge(self, type="mulliken"):
        if type == "mulliken":
            self.dp_sys.data["charge"] = self.cp2k_out.m_charge.reshape(
                self.nframes, self.natoms
            )
        elif type == "hirshfeld":
            self.dp_sys.data["charge"] = self.cp2k_out.h_charge.reshape(
                self.nframes, self.natoms
            )
        else:
            raise NotImplementedError

    def get_sel_ids(
        self,
        type_map: List[str],
        sel_type: List[int],
    ):
        symbols = np.array(self.dp_sys.data["atom_names"])[
            self.dp_sys.data["atom_types"]
        ]
        sel_ids = []
        for ii in sel_type:
            atype = type_map[ii]
            sel_ids.append(np.where(symbols == atype)[0])
        sel_ids = np.concatenate(sel_ids)
        return sel_ids
