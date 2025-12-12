# SPDX-License-Identifier: LGPL-3.0-or-later
"""
Utility functions for computational chemistry and materials science.

This module provides a collection of utility functions for various tasks including:
- Dictionary manipulation and file I/O operations
- CP2K input file generation
- Density and concentration calculations
- Water molecule analysis and manipulation
- Lennard-Jones parameter calculations
- Coordinate number calculations

The functions are designed to work with common scientific computing libraries
such as NumPy, ASE, and MDAnalysis.
"""
import collections
import csv
import json
import os
import pickle
from typing import Any, List, Optional, Union

import h5py
import numpy as np
import yaml
from ase import Atoms
from ase.geometry import wrap_positions
from scipy import constants

try:
    from MDAnalysis.lib.distances import distance_array, minimize_vectors
except ImportError:
    import warnings

    warnings.warn("calc_coord_number and wrap_water cannot be used without MDAnalysis", stacklevel=2)
    
import contextlib


def iterdict(
    input_dict: dict[str, Any], out_list: list[str], loop_idx: int
) -> list[str]:
    """
    Recursively generate a list of strings for CP2K input file formatting.

    This function processes a nested dictionary structure and converts it into
    a formatted list of strings suitable for generating CP2K input files.

    Parameters
    ----------
    input_dict : Dict[str, Any]
        Dictionary containing CP2K input parameters
    out_list : List[str]
        List of strings for printing (modified in-place)
    loop_idx : int
        Record of loop levels in recursion

    Returns
    -------
    List[str]
        Formatted list of strings for CP2K input file
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
    return out_list


def update_dict(old_d: dict[str, Any], update_d: dict[str, Any]) -> None:
    """
    Recursively update a dictionary with values from another dictionary.

    Source: dpgen.generator.lib.cp2k

    Parameters
    ----------
    old_d : Dict[str, Any]
        Original dictionary to be updated
    update_d : Dict[str, Any]
        Dictionary containing update values
    """
    for k in update_d:
        if (
            k in old_d
            and isinstance(old_d[k], dict)
            and isinstance(update_d[k], collections.abc.Mapping)
        ):
            update_dict(old_d[k], update_d[k])
        else:
            old_d[k] = update_d[k]


def dict_to_cp2k_input(input_dict: dict[str, Any]) -> str:
    """
    Convert a dictionary to CP2K input file format.

    Parameters
    ----------
    input_dict : Dict[str, Any]
        Dictionary containing CP2K input parameters

    Returns
    -------
    str
        Formatted CP2K input string
    """
    input_str = iterdict(input_dict, out_list=["\n"], loop_idx=0)
    s = "\n".join(input_str)
    s = s.strip("\n")
    return s


def symlink(src: str, _dst: str) -> None:
    """
    Create a symbolic link.

    Parameters
    ----------
    src : str
        Source file path
    _dst : str
        Destination path
    """
    dst = os.path.abspath(_dst)
    os.symlink(src, dst)


def save_dict(d: dict[str, Any], fname: str, fmt: Optional[str] = None) -> None:
    """
    Save a dictionary to a file in various formats.

    Args:
        d: Dictionary to save
        fname: Output filename
        fmt: File format (json, csv, pkl, hdf5, yaml). If None, inferred from extension.

    Raises
    ------
        KeyError: If the format is not supported
    """
    if fmt is None:
        fmt = os.path.splitext(fname)[1][1:]
    try:
        globals()[f"save_dict_{fmt}"](d, fname)
    except KeyError as exc:
        raise KeyError(f"Unknown format {fmt}") from exc


def save_dict_json(d: dict[str, Any], fname: str) -> None:
    """
    Save dictionary to JSON file.

    Parameters
    ----------
    d : Dict[str, Any]
        Dictionary to save
    fname : str
        Output filename
    """
    with open(fname, "w", encoding="UTF-8") as f:
        json.dump(d, f, indent=4)


def save_dict_csv(d: dict[str, Any], fname: str) -> None:
    """
    Save dictionary to CSV file.

    Parameters
    ----------
    d : Dict[str, Any]
        Dictionary to save
    fname : str
        Output filename
    """
    with open(fname, "w", newline="", encoding="UTF-8") as f:
        writer = csv.DictWriter(f, fieldnames=d.keys())
        writer.writeheader()
        writer.writerow(d)


def save_dict_pkl(d: dict[str, Any], fname: str) -> None:
    """
    Save dictionary to pickle file.

    Parameters
    ----------
    d : Dict[str, Any]
        Dictionary to save
    fname : str
        Output filename
    """
    with open(fname, "wb") as f:
        pickle.dump(d, f)


def save_dict_hdf5(d: dict[str, Any], fname: str) -> None:
    """
    Save dictionary to HDF5 file.

    Parameters
    ----------
    d : Dict[str, Any]
        Dictionary to save
    fname : str
        Output filename
    """
    with h5py.File(fname, "a") as f:
        n = len(f)
        dts = f.create_group(f"{n:02d}")
        for k, v in d.items():
            dts.create_dataset(k, data=v)


def save_dict_yaml(d: dict[str, Any], fname: str) -> None:
    """
    Save dictionary to YAML file.

    Parameters
    ----------
    d : Dict[str, Any]
        Dictionary to save
    fname : str
        Output filename
    """
    with open(fname, "w", encoding="UTF-8") as f:
        yaml.safe_dump(d, f)


def load_dict(fname: str, fmt: Optional[str] = None) -> dict[str, Any]:
    """
    Load a dictionary from a file in various formats.

    Parameters
    ----------
    fname : str
        Input filename
    fmt : Optional[str], optional
        File format (json, csv, pkl, hdf5, yaml). If None, inferred from extension.

    Returns
    -------
    Dict[str, Any]
        Loaded dictionary

    Raises
    ------
    KeyError
        If the format is not supported
    NotImplementedError
        If HDF5 format is requested (not implemented)
    """
    if fmt is None:
        fmt = os.path.splitext(fname)[1][1:]
    try:
        return globals()[f"load_dict_{fmt}"](fname)
    except KeyError as exc:
        raise KeyError(f"Unknown format {fmt}") from exc


def load_dict_json(fname: str) -> dict[str, Any]:
    """
    Load dictionary from JSON file.

    Parameters
    ----------
    fname : str
        Input filename

    Returns
    -------
    Dict[str, Any]
        Loaded dictionary
    """
    with open(fname, encoding="UTF-8") as f:
        d = json.load(f)
    return d


def load_dict_csv(fname: str) -> dict[str, str]:
    """
    Load dictionary from CSV file.

    Parameters
    ----------
    fname : str
        Input filename

    Returns
    -------
    Dict[str, str]
        Loaded dictionary
    """
    with open(fname, encoding="UTF-8") as f:
        data = csv.reader(f)
        d = {rows[0]: rows[1] for rows in data}
    return d


def load_dict_pkl(fname: str) -> dict[str, Any]:
    """
    Load dictionary from pickle file.

    Parameters
    ----------
    fname : str
        Input filename

    Returns
    -------
    Dict[str, Any]
        Loaded dictionary
    """
    with open(fname, "rb") as f:
        d = pickle.load(f)
    return d


def load_dict_yaml(fname: str) -> dict[str, Any]:
    """
    Load dictionary from YAML file.

    Parameters
    ----------
    fname : str
        Input filename

    Returns
    -------
    Dict[str, Any]
        Loaded dictionary
    """
    with open(fname, encoding="UTF-8") as f:
        return yaml.safe_load(f)


def load_dict_hdf5(fname: str) -> dict[str, Any]:
    """
    Load dictionary from HDF5 file (not implemented).

    Parameters
    ----------
    fname : str
        Input filename

    Raises
    ------
    NotImplementedError
        This function is not implemented
    """
    raise NotImplementedError


def get_efields(delta_v: float, layer_thicknesses: list[float], eps: list[float]) -> np.ndarray:
    """
    Calculate electric fields from voltage differences and dielectric properties.

    Parameters
    ----------
    delta_v : float
        Voltage difference
    layer_thicknesses : List[float]
        List of layer thicknesses
    eps : List[float]
        List of dielectric constants

    Returns
    -------
    np.ndarray
        Array of electric field values
    """
    r_field = 1.0 / np.array(eps)
    _delta_v = np.sum(np.array(layer_thicknesses) * r_field)
    v_coeff = delta_v / _delta_v
    return r_field * v_coeff


def safe_makedirs(dname: str) -> None:
    """
    Create directory if it doesn't exist.

    Parameters
    ----------
    dname : str
        Directory path to create
    """
    if not os.path.exists(dname):
        os.makedirs(dname)


def safe_symlink(src: str, dst: str, **kwargs) -> None:
    """
    Create symbolic link if it doesn't exist.

    Parameters
    ----------
    src : str
        Source file path
    dst : str
        Destination path
    **kwargs
        Additional arguments for os.symlink
    """
    with contextlib.suppress(OSError):
        os.symlink(src, dst, **kwargs)


def calc_density(
    n: Union[int, np.ndarray], v: Union[float, np.ndarray], mol_mass: float
) -> Union[float, np.ndarray]:
    """
    Calculate density (g/cm³) from the number of particles.

    Parameters
    ----------
    n : int or np.ndarray
        Number of particles
    v : float or np.ndarray
        Volume in Å³
    mol_mass : float
        Molar mass in g/mol

    Returns
    -------
    float or np.ndarray
        Density in g/cm³
    """
    rho = (n / constants.Avogadro * mol_mass) / (
        v * (constants.angstrom / constants.centi) ** 3
    )
    return rho


def calc_water_density(
    n: Union[int, np.ndarray], v: Union[float, np.ndarray]
) -> Union[float, np.ndarray]:
    """
    Calculate water density using water molar mass (18.015 g/mol).

    Parameters
    ----------
    n : int or np.ndarray
        Number of water molecules
    v : float or np.ndarray
        Volume in Å³

    Returns
    -------
    float or np.ndarray
        Water density in g/cm³
    """
    return calc_density(n, v, 18.015)


def calc_number(rho: Union[float, np.ndarray], v: Union[float, np.ndarray], mol_mass: float) -> int:
    """
    Calculate number of particles from density and volume.

    Parameters
    ----------
    rho : float or np.ndarray
        Density in g/cm³
    v : float or np.ndarray
        Volume in Å³
    mol_mass : float
        Molar mass in g/mol

    Returns
    -------
    int
        Number of particles (rounded to integer)
    """
    n = (
        rho
        * (v * (constants.angstrom / constants.centi) ** 3)
        * constants.Avogadro
        / mol_mass
    )
    return int(n)


def calc_water_number(rho: Union[float, np.ndarray], v: Union[float, np.ndarray]) -> int:
    """
    Calculate number of water molecules from density and volume.

    Parameters
    ----------
    rho : float or np.ndarray
        Water density in g/cm³
    v : float or np.ndarray
        Volume in Å³

    Returns
    -------
    int
        Number of water molecules
    """
    return calc_number(rho, v, 18.015)


def calc_coord_number(
    atoms: Atoms,
    c_ids: np.ndarray,
    neigh_ids: np.ndarray,
    cutoff: Optional[float] = None,
    voronoi: bool = False,
) -> np.ndarray:
    """
    Calculate coordination numbers for specified atoms.

    Parameters
    ----------
    atoms : Atoms
        ASE Atoms object
    c_ids : np.ndarray
        Indices of central atoms
    neigh_ids : np.ndarray
        Indices of neighbor atoms
    cutoff : Optional[float], optional
        Distance cutoff for coordination calculation
    voronoi : bool, optional
        If True, use Voronoi method instead of cutoff

    Returns
    -------
    np.ndarray
        Array of coordination numbers
    """
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


def calc_water_coord_number(atoms: Atoms) -> np.ndarray:
    """
    Calculate coordination numbers for water molecules (O-H coordination).

    Parameters
    ----------
    atoms : Atoms
        ASE Atoms object containing water molecules

    Returns
    -------
    np.ndarray
        Array of coordination numbers for oxygen atoms
    """
    atype = np.array(atoms.get_chemical_symbols())
    c_ids = np.where(atype == "O")[0]
    neigh_ids = np.where(atype == "H")[0]
    return calc_coord_number(atoms, c_ids, neigh_ids, 1.3)


def check_water(atoms: Atoms) -> bool:
    """
    Check if all water molecules have correct coordination (2 H per O).

    Parameters
    ----------
    atoms : Atoms
        ASE Atoms object containing water molecules

    Returns
    -------
    bool
        True if all water molecules are correct, False otherwise
    """
    cns = calc_water_coord_number(atoms)
    flags = cns == 2
    return False not in flags


def wrap_water(atoms: Atoms) -> Atoms:
    """
    Make water molecules whole (keep O and H atoms together).

    This function ensures that water molecules are not split across periodic
    boundaries. Works for pure water systems.

    Parameters
    ----------
    atoms : Atoms
        ASE Atoms object containing water molecules

    Returns
    -------
    Atoms
        New Atoms object with wrapped water molecules
    """
    atoms = atoms.copy()
    oxygen_mask = atoms.symbols == "O"
    hydrogen_mask = atoms.symbols == "H"
    other_mask = np.logical_not(np.logical_or(oxygen_mask, hydrogen_mask))

    coords = wrap_positions(atoms.get_positions(), atoms.get_cell())
    cellpar = atoms.cell.cellpar()
    oxygen_coords = coords[oxygen_mask]
    hydrogen_coords = coords[hydrogen_mask]
    new_atoms = Atoms(cell=atoms.cell, pbc=atoms.pbc)
    for oxygen_coord in oxygen_coords:
        oxygen_coord = np.reshape(oxygen_coord, (1, 3))
        ds = distance_array(oxygen_coord, hydrogen_coords, box=cellpar)
        mask = ds.reshape(-1) < 1.3
        cn = np.sum(mask)
        coords_rel = hydrogen_coords[mask].reshape(-1, 3) - oxygen_coord
        coords_rel = minimize_vectors(coords_rel, box=cellpar)
        _coords = np.concatenate(
            (
                oxygen_coord,
                oxygen_coord + coords_rel,
            ),
            axis=0,
        )
        new_atoms.extend(Atoms(f"OH{cn}", positions=_coords))
        # remove selected hydrogen atoms
        hydrogen_coords = hydrogen_coords[~mask]
    assert len(hydrogen_coords) == 0
    new_atoms.extend(atoms[other_mask])
    return new_atoms


def calc_molar_concentration(
    number_density: Union[int, float, np.ndarray], grid_volume: Union[float, np.ndarray]
) -> Union[float, np.ndarray]:
    """
    Calculate molar concentration from number density.

    Parameters
    ----------
    number_density : int, float, or np.ndarray
        Number of particles per grid
    grid_volume : float or np.ndarray
        Volume of each grid in Å³ (angstrom cubed)

    Returns
    -------
    float or np.ndarray
        Molar concentration in mol/L (moles of ions per liter of solution)
    """
    # Convert grid volume from Å³ to cm³ (1 Å³ = 1e-24 cm³)
    volume_cm3 = grid_volume * 1e-24

    # Convert volume from cm³ to liters (1 L = 1000 cm³)
    volume_liters = volume_cm3 / 1000.0

    # Calculate moles of ions in each grid
    moles = number_density / constants.Avogadro

    # Calculate molar concentration (mol/L)
    molar_concentration = moles / volume_liters

    return molar_concentration

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

def calc_lj_params(ks: list[str]) -> None:
    """
    Calculate and print Lennard-Jones parameters for element pairs.

    This function calculates mixed Lennard-Jones parameters using the
    Lorentz-Berthelot combining rules and prints them to stdout.

    Parameters
    ----------
    ks : List[str]
        List of element symbols
    """
    n_elements = len(ks)
    _sigma = []
    _epsilon = []
    for i in range(n_elements):
        try:
            sigma_i = lj_params[ks[i]]["sigma"]
            epsilon_i = lj_params[ks[i]]["epsilon"]
            _sigma.append(sigma_i)
        except KeyError as exc:
            raise KeyError(f"sigma for {ks[i]} not found") from exc
        for j in range(i, n_elements):
            try:
                sigma_j = lj_params[ks[j]]["sigma"]
                epsilon_j = lj_params[ks[j]]["epsilon"]
            except KeyError as exc:
                raise KeyError(f"sigma for {ks[j]} not found") from exc
            print(ks[i], ks[j], (sigma_i + sigma_j) / 2, np.sqrt(epsilon_i * epsilon_j))


def get_bins_from_bin_edge(bin_edges: Union[np.ndarray, List[float]]) -> np.ndarray:
    """
    Get bin centers from bin edges.

    Parameters
    ----------
    bin_edges : np.ndarray or List[float]
        Array of bin edges

    Returns
    -------
    np.ndarray
        Array of bin centers
    """
    bin_edges = np.reshape(bin_edges, (-1,))
    bins = bin_edges[:-1] + np.diff(bin_edges) / 2
    return bins
