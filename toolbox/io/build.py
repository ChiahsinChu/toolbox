import MDAnalysis as mda
import mdapackmol
from ase import build, Atoms

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

        system = mdapackmol.packmol([
            mdapackmol.PackmolStructure(
                water,
                number=self.n_wat,
                instructions=["inside box %s" % self.box_string])
        ])
        system.atoms.write(fname, **kwargs)
        os.remove("water.pdb")


class Interface:
    def __init__(self, slab) -> None:
        self.slab = slab

    def run(self, l_water=30):
        coord = self.slab.get_positions()
        z = coord[:, 2]
        l_slab = z.max() - z.min()
        a = self.slab.get_cell()[0][0]
        b = self.slab.get_cell()[1][1]
        c = l_slab + l_water
        new_cell = [a, b, c]
        shift_z = -z.min() - l_slab / 2
        coord[:, 2] += shift_z
        slab = Atoms(self.slab.get_chemical_symbols(),
                     positions=coord,
                     cell=new_cell,
                     pbc=True)
        slab.wrap()

        WaterBox(box=[[0, a], [0, b], [l_slab / 2, l_slab / 2 + l_water]],
                 slit=[1., 1., 2.5]).run("waterbox.xyz")
        waterbox = io.read("waterbox.xyz")

        self.atoms = waterbox + slab
        self.atoms.set_cell(new_cell)
        self.atoms.set_pbc(True)
        os.remove("waterbox.xyz")
        
    def relax(self):
        # TODO: add relaxation step with lammps
        
        pass    

