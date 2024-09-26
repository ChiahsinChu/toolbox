# SPDX-License-Identifier: LGPL-3.0-or-later
import MDAnalysis as mda
import mdapackmol
from ase import Atoms, build

from ..utils import *
from ..utils.utils import calc_water_number


class WaterBox:
    def __init__(self, box=None, rho=1.0, slit=1.0) -> None:
        box = np.array(box)
        assert len(box) == 3
        _box = np.reshape(box, (3, -1))
        if len(_box[0]) == 1:
            box = np.concatenate([np.zeros((3, 1)), _box], axis=-1)

        volume = np.prod(np.diff(box, axis=-1))
        self.n_wat = calc_water_number(rho, volume)

        slit = np.reshape(slit, (-1))
        if (len(slit) == 1) or (len(slit) == 3):
            box[:, 0] += slit
            box[:, 1] -= slit
        else:
            raise AttributeError("")
        self.box_string = np.array2string(np.transpose(box).flatten())[1:-1]

    def run(self, fname, **kwargs):
        atoms = build.molecule("H2O")
        io.write("water.pdb", atoms)
        water = mda.Universe("water.pdb")

        system = mdapackmol.packmol(
            [
                mdapackmol.PackmolStructure(
                    water,
                    number=self.n_wat,
                    instructions=["inside box %s" % self.box_string],
                )
            ]
        )
        system.atoms.write(fname, **kwargs)
        os.remove("water.pdb")


class Interface:
    def __init__(self, slab) -> None:
        self.slab = slab

    def run(self, l_water: float = 30):
        coord = self.slab.get_positions()
        z = coord[:, 2]
        l_slab = z.max() - z.min()
        a = self.slab.get_cell()[0][0]
        b = self.slab.get_cell()[1][1]
        c = l_slab + l_water
        new_cell = [a, b, c]
        shift_z = -z.min() - l_slab / 2
        coord[:, 2] += shift_z
        slab = Atoms(
            self.slab.get_chemical_symbols(), positions=coord, cell=new_cell, pbc=True
        )
        slab.wrap()

        WaterBox(
            box=[[0, a], [0, b], [l_slab / 2, l_slab / 2 + l_water]],
            slit=[1.0, 1.0, 2.5],
        ).run("waterbox.xyz")
        waterbox = io.read("waterbox.xyz")

        self.atoms = waterbox + slab
        self.atoms.set_cell(new_cell)
        self.atoms.set_pbc(True)
        os.remove("waterbox.xyz")

    def relax(self):
        # TODO: add relaxation step with lammps

        pass


def add_ion(atoms, ion, region, cutoff=2.0):
    """
    Example

        atoms = io.read("coord.xyz")
        ion = Atoms("K")
        atoms = add_ion(atoms, ion, [8.0, 13.0])
        print(atoms)
    """
    is_overlap = 1
    while is_overlap != 0:
        random_positions = get_region_random_location(atoms, region)
        # print(random_positions)
        ion.set_positions(random_positions.reshape(1, 3))
        new_atoms = atoms.copy()
        new_atoms.extend(ion)
        is_overlap = is_atom_overlap(new_atoms, -1, r=cutoff)
    return new_atoms


def add_water(atoms, region, cutoff=2.0):
    """
    Example

        atoms = io.read("coord.xyz")
        ion = Atoms("K")
        atoms = add_ion(atoms, ion, [8.0, 13.0])
        print(atoms)
    """
    water = build.molecule("H2O")
    coords = water.get_positions()
    coords -= coords[0].reshape(1, 3)
    rotation_matrix = random_rotation_matrix()
    coords = np.dot(coords, rotation_matrix)

    is_overlap = 1
    while is_overlap != 0:
        random_positions = get_region_random_location(atoms, region)
        # print(random_positions)
        water.set_positions(coords + random_positions.reshape(1, 3))
        new_atoms = atoms.copy()
        new_atoms.extend(water)
        is_overlap = sum(
            [is_atom_overlap(new_atoms, -i, 3, r=cutoff) for i in range(1, 3 + 1)]
        )
    return new_atoms


def get_region_random_location(atoms, region, extent=0.9):
    x_region = [atoms.get_cell()[0][0] * (1 - extent), atoms.get_cell()[0][0] * extent]
    y_region = [atoms.get_cell()[1][1] * (1 - extent), atoms.get_cell()[1][1] * extent]

    location_x = np.random.uniform(x_region[0], x_region[1])
    location_y = np.random.uniform(y_region[0], y_region[1])
    location_z = np.random.uniform(region[0], region[1])

    return np.array([location_x, location_y, location_z])


def is_atom_overlap(atoms, index, atom_num=1, r=1.5):
    """
    judge the overlap situation for such atom and atoms in molecule or ions by index
    Args:
        atoms: atoms for adding
        index: index in added atoms. example: O and H1 and H2 index for such added H2O
        atom_num: Number of atoms of an added single molecule, such as 3 for added H2O, 1 for added F
        r: recommended: 1.0-1.8: set up this factor to specify the minimum distance between two atoms

    Returns:

    """
    distance = atoms.get_all_distances(mic=True)[index]
    min_index = np.argsort(distance)[atom_num]

    if distance[min_index] < r:
        overlap = True
    else:
        overlap = False
    return overlap


def random_rotation_matrix():
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
