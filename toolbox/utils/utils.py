# SPDX-License-Identifier: LGPL-3.0-or-later
import collections
import csv
import json
import os
import pickle
from typing import Optional

import h5py
import numpy as np
import yaml
from ase import Atoms
from ase.geometry import wrap_positions
from scipy import constants


def iterdict(input_dict, out_list, loop_idx):
    """
    recursively generate a list of strings for further
    print out CP2K input file

    Args:
        input_dict: dictionary for CP2K input parameters
        out_list: list of strings for printing
        loop_idx: record of loop levels in recursion
    Return:
        out_list
    """
    if len(out_list) == 0:
        out_list.append("\n")
    start_idx = len(out_list) - loop_idx - 2
    for k, v in input_dict.items():
        k = str(k)  # cast key into string
        # if value is dictionary
        if isinstance(v, dict):
            out_list.insert(-1 - loop_idx, "  " * loop_idx + "&" + k)
            out_list.insert(-1 - loop_idx, "  " * loop_idx + "&END " + k)
            iterdict(v, out_list, loop_idx + 1)
        # if value is list
        elif isinstance(v, list):
            if isinstance(v[0], dict):
                for _v in v:
                    out_list.insert(-1 - loop_idx, "  " * loop_idx + "&" + k)
                    out_list.insert(-1 - loop_idx, "  " * loop_idx + "&END " + k)
                    iterdict(_v, out_list, loop_idx + 1)
                # print(loop_idx)
                # print(input_dict)
                # print(out_list)
            else:
                for _v in v:
                    _v = str(_v)
                    out_list.insert(-1 - loop_idx, "  " * loop_idx + k + " " + _v)
        # if value is other type, e.g., int/float/str
        else:
            v = str(v)
            if k == "_":
                out_list[start_idx] = out_list[start_idx] + " " + v
            else:
                out_list.insert(-1 - loop_idx, "  " * loop_idx + k + " " + v)
                # out_list.insert(-1-loop_idx, v)
    return out_list


def update_dict(old_d, update_d):
    """
    source: dpgen.generator.lib.cp2k

    a method to recursive update dict
    :old_d: old dictionary
    :update_d: some update value written in dictionary form
    """
    for k, v in update_d.items():
        if (
            k in old_d
            and isinstance(old_d[k], dict)
            and isinstance(update_d[k], collections.abc.Mapping)
        ):
            update_dict(old_d[k], update_d[k])
        # elif (k in old_d and isinstance(old_d[k], list)
        #       and isinstance(update_d[k], list)):
        #     old_d[k].extend(update_d[k])
        else:
            old_d[k] = update_d[k]


def dict_to_cp2k_input(input_dict):
    input_str = iterdict(input_dict, out_list=["\n"], loop_idx=0)
    s = "\n".join(input_str)
    s = s.strip("\n")
    return s


def symlink(src, _dst):
    dst = os.path.abspath(_dst)
    os.symlink(src, dst)


def save_dict(d: dict, fname: str, fmt=None):
    if fmt is None:
        fmt = os.path.splitext(fname)[1][1:]
    try:
        globals()["save_dict_%s" % fmt](d, fname)
    except:
        raise AttributeError("Unknown format %s" % fmt)


def save_dict_json(d: dict, fname: str):
    with open(fname, "w", encoding="UTF-8") as f:
        json.dump(d, f, indent=4)


def save_dict_csv(d: dict, fname: str):
    with open(fname, "w", newline="", encoding="UTF-8") as f:
        writer = csv.DictWriter(f, fieldnames=d.keys())
        writer.writeheader()
        writer.writerow(d)


def save_dict_pkl(d: dict, fname: str):
    with open(fname, "wb") as f:
        pickle.dump(d, f)


def save_dict_hdf5(d: dict, fname: str):
    with h5py.File(fname, "a") as f:
        n = len(f)
        dts = f.create_group("%02d" % n)
        for k, v in d.items():
            dts.create_dataset(k, data=v)


def save_dict_yaml(d: dict, fname: str):
    with open(fname, "w", encoding="UTF-8") as f:
        yaml.safe_dump(d, f)


def load_dict(fname: str, fmt: str = None):
    if fmt is None:
        fmt = os.path.splitext(fname)[1][1:]
    try:
        return globals()["load_dict_%s" % fmt](fname)
    except KeyError:
        raise KeyError("Unknown format %s" % fmt)


def load_dict_json(fname: str):
    with open(fname, "r", encoding="UTF-8") as f:
        d = json.load(f)
    return d


def load_dict_csv(fname: str):
    with open(fname, "r", encoding="UTF-8") as f:
        data = csv.reader(f)
        d = {rows[0]: rows[1] for rows in data}
    return d


def load_dict_pkl(fname: str):
    with open(fname, "rb") as f:
        d = pickle.load(f)
    return d


def load_dict_yaml(fname: str):
    with open(fname, "r", encoding="UTF-8") as f:
        return yaml.safe_load(f)


def load_dict_hdf5(fname: str):
    raise NotImplementedError


def get_efields(DeltaV, l: list, eps: list):
    r_field = 1.0 / np.array(eps)
    _delta_v = np.sum(np.array(l) * r_field)
    v_coeff = DeltaV / _delta_v
    return r_field * v_coeff


def safe_makedirs(dname):
    if not os.path.exists(dname):
        os.makedirs(dname)


def safe_symlink(src, dst, **kwargs):
    try:
        os.symlink(src, dst, **kwargs)
    except:
        pass


def calc_density(n, v, mol_mass: float):
    """
    calculate density (g/cm^3) from the number of particles

    Parameters
    ----------
    n : int or array
        number of particles
    v : float or array
        volume
    mol_mass : float
        mole mass in g/mol
    """
    rho = (n / constants.Avogadro * mol_mass) / (
        v * (constants.angstrom / constants.centi) ** 3
    )
    return rho


def calc_water_density(n, v):
    return calc_density(n, v, 18.015)


def calc_number(rho, v, mol_mass: float):
    n = (
        rho
        * (v * (constants.angstrom / constants.centi) ** 3)
        * constants.Avogadro
        / mol_mass
    )
    return int(n)


def calc_water_number(rho, v):
    return calc_number(rho, v, 18.015)


def calc_coord_number(
    atoms,
    c_ids,
    neigh_ids,
    cutoff: Optional[float] = None,
    voronoi: bool = False,
):
    from MDAnalysis.lib.distances import distance_array

    p = atoms.get_positions()
    p_c = p[c_ids]
    p_n = p[neigh_ids]
    results = np.empty((len(c_ids), len(neigh_ids)), dtype=np.float64)
    distance_array(p_c, p_n, box=atoms.cell.cellpar(), result=results)
    if cutoff is None and voronoi:
        # use voronoi method
        out = np.unique(np.argmin(results, axis=0), return_counts=True)
        cns = np.zeros(len(p_c), dtype=np.int32)
        cns[out[0]] = out[1]
    else:
        # use cutoff method
        cns = np.count_nonzero(results <= cutoff, axis=1)
    return cns


def calc_water_coord_number(atoms):
    atype = np.array(atoms.get_chemical_symbols())
    c_ids = np.where(atype == "O")[0]
    neigh_ids = np.where(atype == "H")[0]
    return calc_coord_number(atoms, c_ids, neigh_ids, 1.3)


def check_water(atoms):
    cns = calc_water_coord_number(atoms)
    flags = cns == 2
    return False not in flags


def wrap_water(atoms):
    "make water atoms together, pure water only"
    from MDAnalysis.lib.distances import distance_array, minimize_vectors

    atoms = atoms.copy()
    oxygen_mask = atoms.symbols == "O"
    hydrogen_mask = atoms.symbols == "H"
    other_mask = np.logical_not(np.logical_or(oxygen_mask, hydrogen_mask))

    coords = wrap_positions(atoms.get_positions(), atoms.get_cell())
    cellpar = atoms.cell.cellpar()
    oxygen_coords = coords[oxygen_mask]
    hydrogen_coords = coords[hydrogen_mask]
    assert len(oxygen_coords) * 2 == len(hydrogen_coords)
    dist_mat = distance_array(oxygen_coords, hydrogen_coords, box=cellpar)
    new_atoms = Atoms(cell=atoms.cell, pbc=atoms.pbc)
    for ii, ds in enumerate(dist_mat):
        mask = ds < 1.3
        cn = np.sum(mask)
        assert cn == 2
        # print(cn)
        coords_rel = hydrogen_coords[mask] - oxygen_coords[ii]
        coords_rel = minimize_vectors(coords_rel, box=cellpar)
        _coords = np.concatenate(
            (
                oxygen_coords[ii].reshape(1, 3),
                oxygen_coords[ii].reshape(1, 3) + coords_rel,
            ),
            axis=0,
        )
        new_atoms.extend(Atoms("OHH", positions=_coords))
    new_atoms.extend(atoms[other_mask])
    return new_atoms


def calc_lj_params(ks):
    """
    ks
        list of element symbols
    """
    from ..io.template import lj_params

    n_elements = len(ks)
    _sigma = []
    _epsilon = []
    for i in range(n_elements):
        try:
            sigma_i = lj_params[ks[i]]["sigma"]
            epsilon_i = lj_params[ks[i]]["epsilon"]
            _sigma.append(sigma_i)
        except:
            raise KeyError("sigma for %s not found" % ks[i])
        for j in range(i, n_elements):
            try:
                sigma_j = lj_params[ks[j]]["sigma"]
                epsilon_j = lj_params[ks[j]]["epsilon"]
            except KeyError:
                raise KeyError("sigma for %s not found" % ks[j])
            print(ks[i], ks[j], (sigma_i + sigma_j) / 2, np.sqrt(epsilon_i * epsilon_j))
