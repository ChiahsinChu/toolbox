# SPDX-License-Identifier: LGPL-3.0-or-later
"""Molecular system building module.

This module provides classes and functions for building molecular systems,
including solution boxes, interfaces, and adding ions or molecules
to existing systems.
"""

import glob
import os
from typing import List, Optional, Union

import MDAnalysis as mda
import mdapackmol
import numpy as np
from ase import Atoms, build, io
from MDAnalysis.lib.distances import distance_array

from ..utils.utils import calc_water_number


class SolutionBox:
    """Base class for solution box configuration.
    
    This class defines the boundary conditions for solution
    systems with optional slit geometry.
    
    Parameters
    ----------
    boundary : Union[List, np.ndarray]
        Box dimensions as [a, b, c] or [[a1, a2], [b1, b2], [c1, c2]]
    slit : Union[float, List, np.ndarray], optional
        Slit distance or dimensions, by default 1.0
    """

    def __init__(
        self,
        boundary: Union[List, np.ndarray],
        slit: Union[float, List, np.ndarray] = 1.0,
    ) -> None:
        """Initialize SolutionBox.
        
        Parameters
        ----------
        boundary : Union[List, np.ndarray]
            Box dimensions
        slit : Union[float, List, np.ndarray], optional
            Slit distance or dimensions, by default 1.0
        """
        boundary = np.array(boundary)
        assert len(boundary) == 3
        _boundary = np.reshape(boundary, (3, -1))
        # [[a1, a2], [b1, b2], [c1, c2]]
        if len(_boundary[0]) == 1:
            self.boundary = np.concatenate([np.zeros((3, 1)), _boundary], axis=-1)
        elif len(_boundary[0]) == 2:
            self.boundary = _boundary
        else:
            raise AttributeError(
                "boundary must be [a, b, c] or [[a1, a2], [b1, b2], [c1, c2]]"
            )

        slit = np.reshape(slit, (-1))
        b = self.boundary.copy()
        if (len(slit) == 1) or (len(slit) == 3):
            b[:, 0] += slit
            b[:, 1] -= slit
        else:
            raise AttributeError("")
        self.boundary_string = np.array2string(np.transpose(b).flatten())[1:-1]


class WaterBox(SolutionBox):
    """Class for water box configuration.
    
    This class extends SolutionBox to specifically handle
    water molecules with density calculations.
    
    Parameters
    ----------
    boundary : Union[List, np.ndarray]
        Box dimensions as [a, b, c] or [[a1, a2], [b1, b2], [c1, c2]]
    slit : Union[float, List, np.ndarray], optional
        Slit distance or dimensions, by default 1.0
    rho : float
        Water density in g/cm³
    """

    def __init__(
        self,
        boundary: Union[List, np.ndarray],
        slit: Union[float, List, np.ndarray] = 1.0,
        rho: float = 1.0,
    ) -> None:
        """Initialize WaterBox.
        
        Parameters
        ----------
        boundary : Union[List, np.ndarray]
            Box dimensions
        slit : Union[float, List, np.ndarray], optional
            Slit distance or dimensions, by default 1.0
        rho : float, optional
            Water density in g/cm³, by default 1.0
        """
        super().__init__(boundary, slit)

        volume = np.prod(np.diff(self.boundary, axis=-1))
        self.n_wat = calc_water_number(rho, volume)

    def write(
        self,
        fname,
        n_wat: Optional[int] = None,
        seed: int = -1,
        verbose: bool = False,
        **kwargs,
    ):
        """Write water box configuration to file.
        
        Parameters
        ----------
        fname : str
            Output filename
        n_wat : Optional[int], optional
            Number of water molecules. If None, uses calculated value from density
        seed : int, optional
            Random seed for reproducible molecule placement. Use -1 for random behavior
        verbose : bool, optional
            If True, keeps temporary files, by default False
        **kwargs
            Additional keyword arguments passed to ase.io.write
        """
        atoms = build.molecule("H2O")
        io.write("water.pdb", atoms)
        water = mda.Universe("water.pdb")

        if n_wat is None:
            n_wat = self.n_wat

        system = mdapackmol.packmol(
            [
                mdapackmol.PackmolStructure(
                    water,
                    number=n_wat,
                    instructions=[f"inside box {self.boundary_string}", f"seed {seed}"],
                )
            ]
        )
        system.atoms.write(fname, **kwargs)
        if not verbose:
            os.remove("water.pdb")
            os.remove("packmol.stdout")


class ElectrolyteBox(WaterBox):
    """Class for electrolyte solution box configuration.
    
    This class extends WaterBox to handle uniform electrolyte
    solutions with specified solutes and concentrations.
    
    Parameters
    ----------
    boundary : Union[List, np.ndarray]
        Box dimensions as [a, b, c] or [[a1, a2], [b1, b2], [c1, c2]]
    solutes : Dict[str, float]
        Dictionary mapping solute names to concentrations in mol/L
    """

    def __init__(
        self,
        boundary: Union[List, np.ndarray],
        solutes: dict[str, float],
        **kwargs,
    ) -> None:
        """Initialize ElectrolyteBox.
        
        Parameters
        ----------
        boundary : Union[List, np.ndarray]
            Box dimensions
        solutes : Dict[str, float]
            Dictionary mapping solute names to concentrations in mol/L
        **kwargs
            Additional keyword arguments passed to parent class
        """
        super().__init__(boundary, **kwargs)
        self.solutes = solutes

    def write(
        self,
        fname,
        n_wat: Optional[int] = None,
        seed: int = -1,
        verbose=False,
        **kwargs,
    ):
        """Write electrolyte box configuration to file.
        
        Parameters
        ----------
        fname : str
            Output filename
        n_wat : Optional[int], optional
            Number of water molecules. If None, uses calculated value from density
        seed : int, optional
            Random seed for reproducible molecule placement. Use -1 for random behavior
        verbose : bool, optional
            If True, keeps temporary files, by default False
        **kwargs
            Additional keyword arguments passed to ase.io.write
        """
        atoms = build.molecule("H2O")
        io.write("tmp.pdb", atoms)
        water = mda.Universe("tmp.pdb")

        if n_wat is None:
            n_wat = self.n_wat

        structures = [
            mdapackmol.PackmolStructure(
                water,
                number=n_wat,
                instructions=[f"inside box {self.boundary_string}", f"seed {seed}"],
            )
        ]
        for k, v in self.solutes.items():
            fnames = glob.glob(os.path.join(os.path.dirname(__file__), f"ion_structure_lib/*/{k}.*"))
            if len(fnames) > 0:
                atoms = io.read(fnames[0])
            else:
                atoms = Atoms(k, positions=[[0, 0, 0]])
            io.write("tmp.pdb", atoms)
            solute = mda.Universe("tmp.pdb")
            structures.append(
                mdapackmol.PackmolStructure(
                    solute,
                    number=int(n_wat / 55.6 * v),
                    instructions=[f"inside box {self.boundary_string}", f"seed {seed}"],
                )
            )

        system = mdapackmol.packmol(structures)
        system.atoms.write(fname, **kwargs)
        if not verbose:
            os.remove("tmp.pdb")
            os.remove("packmol.stdout")


class Interface:
    """Class for creating interface systems.
    
    This class creates interface systems between a slab
    and water box for electrochemical simulations.
    """
    
    def __init__(
        self,
        slab: Atoms,
        l_water: float = 30,
    ) -> None:
        """Initialize Interface.
        
        Parameters
        ----------
        slab : ase.Atoms
            Slab atoms object
        l_water : float, optional
            Length of water box in Å, by default 30
        """
        # shift half cell and wrap
        coord = slab.get_positions()
        z = coord[:, 2]
        l_slab = z.max() - z.min()
        a = slab.get_cell()[0][0]
        b = slab.get_cell()[1][1]
        c = l_slab + l_water
        new_cell = [a, b, c]
        shift_z = -z.min() - l_slab / 2
        coord[:, 2] += shift_z
        slab = Atoms(
            slab.get_chemical_symbols(), positions=coord, cell=new_cell, pbc=True
        )
        slab.wrap()

        self.slab = slab
        self.boundary = [[0, a], [0, b], [l_slab / 2, l_slab / 2 + l_water]]

    def run(
        self,
        rho: float = 1.0,
        n_wat: Optional[int] = None,
        seed: int = -1,
        sol: Optional[SolutionBox] = None,
        verbose=False,
    ) -> Atoms:
        """Create interface system.
        
        Parameters
        ----------
        rho : float, optional
            Water density in g/cm³, by default 1.0
        n_wat : Optional[int], optional
            Number of water molecules, by default None
        seed : int, optional
            Random seed, by default -1
        sol : Optional[SolutionBox], optional
            Solution box object, by default None
        verbose : bool, optional
            If True, keeps temporary files, by default False
            
        Returns
        -------
        ase.Atoms
            Combined interface system
        """
        if sol is None:
            sol = WaterBox(
                rho=rho,
                boundary=self.boundary,
                slit=[1.0, 1.0, 2.5],
            )
        sol.write("waterbox.xyz", n_wat=n_wat, verbose=verbose, seed=seed)
        waterbox = io.read("waterbox.xyz")
        waterbox.set_cell(self.slab.get_cell())
        waterbox.set_pbc(True)
        waterbox.center(axis=2)

        self.atoms = waterbox + self.slab
        # self.atoms.set_cell(self.slab.get_cell())
        self.atoms.set_pbc(True)
        if not verbose:
            os.remove("waterbox.xyz")
        return self.atoms


def add_ion(atoms, ion, region, cutoff=2.0, max_trial=500):
    """Add an ion to atoms system in specified region.
    
    This function attempts to place an ion in the specified region
    without overlapping with existing atoms.
    
    Parameters
    ----------
    atoms : ase.Atoms
        Existing atoms system
    ion : ase.Atoms
        Ion to add (single atom)
    region : list
        Region boundaries as [x_min, x_max, y_min, y_max, z_min, z_max]
    cutoff : float, optional
        Minimum distance from existing atoms, by default 2.0
    max_trial : int, optional
        Maximum number of placement attempts, by default 500
        
    Returns
    -------
    ase.Atoms or None
        New atoms system with ion added, or None if placement failed
        
    Examples
    --------
    >>> from ase import io
    >>> atoms = io.read("coord.xyz")
    >>> ion = Atoms("K")
    >>> atoms = add_ion(atoms, ion, [8.0, 13.0])
    >>> ion = io.read("ClO4.pdb")
    >>> atoms = add_ion(atoms, ion, [8.0, 13.0])
    >>> print(atoms)
    """
    coords = ion.get_positions()
    cog = np.mean(coords, axis=0)
    coords -= cog.reshape(1, 3)
    rotation_matrix = random_rotation_matrix()
    coords = np.dot(coords, rotation_matrix)

    flag = False
    for _ in range(max_trial):
        random_positions = get_region_random_location(atoms, region)
        random_positions = coords + random_positions.reshape(1, 3)
        ds = distance_array(
            random_positions, atoms.get_positions(), box=atoms.cell.cellpar()
        )
        if ds.min() > cutoff:
            flag = True
            break
    if flag:
        ion.set_positions(random_positions)
        new_atoms = atoms.copy()
        new_atoms.extend(ion)
        return new_atoms
    else:
        raise Warning("Failed to add ion")
        return None


def add_water(atoms, region, cutoff=2.0, max_trial=500):
    """Add a water molecule to atoms system in specified region.
    
    Parameters
    ----------
    atoms : ase.Atoms
        Existing atoms system
    region : list
        Region boundaries as [x_min, x_max, y_min, y_max, z_min, z_max]
    cutoff : float, optional
        Minimum distance from existing atoms, by default 2.0
    max_trial : int, optional
        Maximum number of placement attempts, by default 500
        
    Returns
    -------
    ase.Atoms or None
        New atoms system with water added, or None if placement failed
    """
    water = build.molecule("H2O")
    return add_ion(atoms, water, region, cutoff=cutoff, max_trial=max_trial)


def get_region_random_location(atoms, region, extent=0.9):
    """Generate random location within specified region.
    
    Parameters
    ----------
    atoms : ase.Atoms
        Atoms system for box dimensions
    region : list
        Region boundaries as [x_min, x_max, y_min, y_max, z_min, z_max]
    extent : float, optional
        Fraction of box extent to use, by default 0.9
        
    Returns
    -------
    np.ndarray
        Random position [x, y, z] within region
    """
    x_region = [atoms.get_cell()[0][0] * (1 - extent), atoms.get_cell()[0][0] * extent]
    y_region = [atoms.get_cell()[1][1] * (1 - extent), atoms.get_cell()[1][1] * extent]

    location_x = np.random.uniform(x_region[0], x_region[1])
    location_y = np.random.uniform(y_region[0], y_region[1])
    location_z = np.random.uniform(region[0], region[1])

    return np.array([location_x, location_y, location_z])


def is_atom_overlap(atoms, index, atom_num=1, r=1.5):
    """Check if atoms overlap with existing atoms.
    
    This function checks if specified atoms in a molecule or ion
    overlap with existing atoms in the system.
    
    Parameters
    ----------
    atoms : ase.Atoms
        Existing atoms system
    index : int
        Index of first atom in the added molecule/ion
    atom_num : int, optional
        Number of atoms in the added molecule/ion, by default 1
    r : float, optional
        Minimum distance between atoms, recommended 1.0-1.8, by default 1.5
        
    Returns
    -------
    bool
        True if atoms overlap, False otherwise
    """
    distance = atoms.get_all_distances(mic=True)[index]
    min_index = np.argsort(distance)[atom_num]

    overlap = distance[min_index] < r
    return overlap


def random_rotation_matrix():
    """Generate a random 3D rotation matrix.
    
    Returns
    -------
    np.ndarray
        3x3 rotation matrix
    """
    random_rotvec = np.random.randn(3)
    random_rotvec /= np.linalg.norm(random_rotvec)
    angle = np.random.rand() * 2 * np.pi
    cos_theta = np.cos(angle)
    sin_theta = np.sin(angle)
    one_minus_cos_theta = 1 - cos_theta

    rot_matrix = np.zeros((3, 3))
    rot_matrix[0, 0] = cos_theta + random_rotvec[0] ** 2 * one_minus_cos_theta
    rot_matrix[1, 1] = cos_theta + random_rotvec[1] ** 2 * one_minus_cos_theta
    rot_matrix[2, 2] = cos_theta + random_rotvec[2] ** 2 * one_minus_cos_theta

    rot_matrix[0, 1] = (
        random_rotvec[0] * random_rotvec[1] * one_minus_cos_theta
        - random_rotvec[2] * sin_theta
    )
    rot_matrix[1, 0] = (
        random_rotvec[0] * random_rotvec[1] * one_minus_cos_theta
        + random_rotvec[2] * sin_theta
    )

    rot_matrix[0, 2] = (
        random_rotvec[0] * random_rotvec[2] * one_minus_cos_theta
        + random_rotvec[1] * sin_theta
    )
    rot_matrix[2, 0] = (
        random_rotvec[0] * random_rotvec[2] * one_minus_cos_theta
        - random_rotvec[1] * sin_theta
    )

    rot_matrix[1, 2] = (
        random_rotvec[1] * random_rotvec[2] * one_minus_cos_theta
        - random_rotvec[0] * sin_theta
    )
    rot_matrix[2, 1] = (
        random_rotvec[1] * random_rotvec[2] * one_minus_cos_theta
        + random_rotvec[0] * sin_theta
    )

    return rot_matrix
