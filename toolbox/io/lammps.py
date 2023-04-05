from ase import io, Atoms
import numpy as np


def read_dump(dump_file="dump.lammpstrj"):
    traj = io.read(dump_file, index=":")
    coords = []
    forces = []
    boxs = []
    for atoms in traj:
        forces.append(atoms.get_forces())
        coords.append(atoms.get_positions())
        boxs.append(atoms.get_cell())
    coords = np.reshape(coords, (len(traj), -1))
    forces = np.reshape(forces, (len(traj), -1))
    boxs = np.reshape(boxs, (len(traj), -1))
    type_list = atoms.get_array('numbers')
    type_list = np.array(type_list) - 1
    return coords, forces, boxs, type_list


class LammpsData:
    def __init__(self, atoms) -> None:
        self.atoms = atoms
        self._setup()

        self.angles = None
        self.bonds = None
        self.dihedrals = None

    def write(self, out_file="system.data", **kwargs):
        specorder = kwargs.get("specorder", None)
        if specorder is not None:
            self.set_atype_from_specorder(specorder)
        atom_style = kwargs.get("atom_style", "full")

        with open(out_file, "w", encoding='utf-8') as f:
            header = self._make_header(out_file)
            f.write(header)
            body = self._make_atoms(atom_style)
            f.write(body)
            if self.bonds is not None:
                f.write("\nBonds\n\n")
                np.savetxt(f, self.bonds, fmt="%d")
            if self.angles is not None:
                f.write("\nAngles\n\n")
                np.savetxt(f, self.angles, fmt="%d")
            if self.dihedrals is not None:
                f.write("\nDihedrals\n\n")
                np.savetxt(f, self.dihedrals, fmt="%d")

    def _make_header(self, out_file):
        cell = self.atoms.cell.cellpar()
        nat = len(self.atoms)
        s = "%s (written by IntDielec)\n\n" % out_file
        s += "%d atoms\n" % nat
        s += "%d atom types\n" % len(np.unique(self.atoms.numbers))
        if self.bonds is not None:
            s += "%d bonds\n" % len(self.bonds)
            s += "%d bond types\n" % len(np.unique(self.bonds[:, 1]))
        if self.angles is not None:
            s += "%d angles\n" % len(self.angles)
            s += "%d angle types\n" % len(np.unique(self.angles[:, 1]))
        if self.dihedrals is not None:
            s += "%d dihedrals\n" % len(self.dihedrals)
            s += "%d dihedral types\n" % len(np.unique(self.dihedrals[:, 1]))
        s += "%.4f %.4f xlo xhi\n%.4f %.4f ylo yhi\n%.4f %.4f zlo zhi\n\n\n" % (
            0.0, cell[0], 0.0, cell[1], 0.0, cell[2])
        return s

    def _make_atoms(self, atom_style):
        """
        full atom_id res_id type q x y z
        atomic ...
        """
        return getattr(self, "_make_atoms_%s" % atom_style)()

    def _make_atoms_full(self):
        """
        full atom_id res_id type q x y z
        """
        s = "Atoms\n\n"
        for atom in self.atoms:
            ii = atom.index
            s += "%d %d %d %.6f %.6f %.6f %.6f\n" % (
                ii + 1, self.res_id[ii], self.atype[ii], self.charges[ii],
                self.positions[ii][0], self.positions[ii][1],
                self.positions[ii][2])
        return s

    def _make_atoms_atomic(self):
        pass

    def set_res_id(self, res_id):
        self.res_id = np.reshape(res_id, (-1))

    def set_atype(self, atype):
        self.atype = np.reshape(atype, (-1))

    def set_atype_from_specorder(self, specorder):
        atype = []
        for ii in self.atoms.get_chemical_symbols():
            atype.append(specorder.index(ii))
        self.atype = np.array(atype, dtype=np.int32) + 1

    def set_bonds(self, bonds):
        self.bonds = np.reshape(bonds, (-1, 4))

    def set_angles(self, angles):
        self.angles = np.reshape(angles, (-1, 5))

    def set_dihedral(self, dihedrals):
        self.dihedrals = np.reshape(dihedrals, (-1, 6))

    def _setup(self):
        self.positions = self.atoms.get_positions()
        if len(self.atoms.get_initial_charges()) > 0:
            self.charges = self.atoms.get_initial_charges().reshape(-1)
        else:
            self.charges = np.zeros(len(self.atoms))

        if hasattr(self, "res_id"):
            assert len(self.res_id) == len(self.atoms)
        else:
            self.res_id = np.zeros(len(self.atoms), dtype=np.int32)


class LammpsDump:
    def __init__(self, traj, type_map) -> None:
        self.traj = traj
        self._set_atype(type_map)
        # self.type_map = type_map

    def write(
        self,
        start=0,
        step=1,
        out_file="out.lammpstrj",
        append=False,
    ):
        if isinstance(self.traj, Atoms):
            self._write_dump(self.traj, ts, out_file, append)
        else:
            nframe = len(self.traj)
            _ts = np.arange(start, step * nframe, step)
            for ts, atoms in zip(_ts, self.traj):
                self._write_dump(atoms, ts, out_file, append=True)

    def _set_atype(self, type_map):
        if isinstance(type_map, list):
            self.atype_dict = {}
            for ii, atype in enumerate(type_map, start=1):
                self.atype_dict[atype] = {}
                self.atype_dict[atype]["type"] = ii
                self.atype_dict[atype]["element"] = atype
        elif isinstance(type_map, dict):
            self.atype_dict = type_map
        else:
            raise AttributeError("Unknown type of type_map")

    def _write_dump(self, atoms, ts, out_file, append):
        if append:
            with open(out_file, "a", encoding='utf-8') as f:
                header = self.make_header(atoms, ts)
                f.write(header)
                body = self.make_body(atoms, self.atype_dict)
                f.write(body)
        else:
            with open(out_file, "w", encoding='utf-8') as f:
                header = self.make_header(atoms, ts)
                f.write(header)
                body = self.make_body(atoms, self.atype_dict)
                f.write(body)

    @staticmethod
    def make_header(atoms, ts):
        bc_dict = {True: "pp", False: "ff"}

        cell = atoms.cell.cellpar()
        bc = atoms.get_pbc()
        nat = len(atoms)
        s = "ITEM: TIMESTEP\n%d\n" % ts
        s += "ITEM: NUMBER OF ATOMS\n"
        s += "%d\n" % nat
        s += "ITEM: BOX BOUNDS %s %s %s\n" % (bc_dict[bc[0]], bc_dict[bc[1]],
                                              bc_dict[bc[2]])
        s += "%.4f %.4f\n%.4f %.4f\n%.4f %.4f\n" % (0.0, cell[0], 0.0, cell[1],
                                                    0.0, cell[2])
        if len(atoms.get_initial_charges()) > 0:
            s += "ITEM: ATOMS id type element x y z q\n"
        else:
            s += "ITEM: ATOMS id type element x y z\n"
        return s

    @staticmethod
    def make_body(atoms, atype_dict):
        if len(atoms.get_initial_charges()) > 0:
            q_flag = True
            charges = atoms.get_initial_charges()
        else:
            q_flag = False
        ps = atoms.get_positions()

        s = ""
        for atom in atoms:
            ii = atom.index
            if q_flag:
                s += "%d %d %s %.6f %.6f %.6f %.6f\n" % (
                    ii + 1, atype_dict[atom.symbol]["type"],
                    atype_dict[atom.symbol]["element"], ps[ii][0], ps[ii][1],
                    ps[ii][2], charges[ii])
            else:
                s += "%d %d %s %.6f %.6f %.6f\n" % (
                    ii + 1, atype_dict[atom.symbol]["type"],
                    atype_dict[atom.symbol]["element"], ps[ii][0], ps[ii][1],
                    ps[ii][2])
        return s


def write_dump(traj,
               type_map,
               start=0,
               step=1,
               out_file="out.lammpstrj",
               append=False):
    if isinstance(type_map, list):
        atype_dict = {}
        for ii, atype in enumerate(type_map, start=1):
            atype_dict[atype] = {}
            atype_dict[atype]["type"] = ii
            atype_dict[atype]["element"] = atype
    elif isinstance(type_map, dict):
        atype_dict = type_map
    else:
        raise AttributeError("Unknown type of type_map")

    if isinstance(traj, Atoms):
        _write_dump(traj, atype_dict, ts, out_file, append)
    else:
        nframe = len(traj)
        _ts = np.arange(start, step * nframe, step)
        for ts, atoms in zip(_ts, traj):
            _write_dump(atoms, atype_dict, ts, out_file, append=True)


def _write_dump(atoms, atype_dict, ts, out_file, append):
    if append:
        with open(out_file, "a", encoding='utf-8') as f:
            header = make_dump_header(atoms, ts)
            f.write(header)
            body = make_dump_body(atoms, atype_dict)
            f.write(body)
    else:
        with open(out_file, "w", encoding='utf-8') as f:
            header = make_dump_header(atoms, ts)
            f.write(header)
            body = make_dump_body(atoms, atype_dict)
            f.write(body)


def make_dump_header(atoms, ts):
    cell = atoms.cell.cellpar()
    nat = len(atoms)
    s = "ITEM: TIMESTEP\n%d\n" % ts
    s += "ITEM: NUMBER OF ATOMS\n"
    s += "%d\n" % nat
    s += "ITEM: BOX BOUNDS pp pp pp\n"
    s += "%.4f %.4f\n%.4f %.4f\n%.4f %.4f\n" % (0.0, cell[0], 0.0, cell[1],
                                                0.0, cell[2])
    if len(atoms.get_initial_charges()) > 0:
        s += "ITEM: ATOMS id type element x y z q\n"
    else:
        s += "ITEM: ATOMS id type element x y z\n"
    return s


def make_dump_body(atoms, atype_dict):
    if len(atoms.get_initial_charges()) > 0:
        q_flag = True
        charges = atoms.get_initial_charges()
    else:
        q_flag = False
    ps = atoms.get_positions()

    s = ""
    for atom in atoms:
        ii = atom.index
        if q_flag:
            s += "%d %d %s %.6f %.6f %.6f %.6f\n" % (
                ii + 1, atype_dict[atom.symbol]["type"],
                atype_dict[atom.symbol]["element"], ps[ii][0], ps[ii][1],
                ps[ii][2], charges[ii])
        else:
            s += "%d %d %s %.6f %.6f %.6f\n" % (
                ii + 1, atype_dict[atom.symbol]["type"],
                atype_dict[atom.symbol]["element"], ps[ii][0], ps[ii][1],
                ps[ii][2])
    return s
