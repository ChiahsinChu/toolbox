# SPDX-License-Identifier: LGPL-3.0-or-later
import re

import numpy as np
from ase import Atoms, io
from ase.calculators.lammps import Prism, convert
from ase.units import fs


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
    type_list = atoms.get_array("numbers")
    type_list = np.array(type_list) - 1
    return coords, forces, boxs, type_list


class LammpsData:
    def __init__(self, atoms) -> None:
        self.atoms = atoms
        self._setup()

        self.angles = None
        self.bonds = None
        self.dihedrals = None
        self.velocities = None

    def write(self, out_file="system.data", **kwargs):
        specorder = kwargs.get("specorder", None)
        if specorder is not None:
            self.set_atype_from_specorder(specorder)
            n_atype = len(specorder)
        else:
            n_atype = len(np.unique(self.atoms.numbers))
        atom_style = kwargs.get("atom_style", "full")

        with open(out_file, "w", encoding="utf-8") as f:
            header = self._make_header(out_file, n_atype)
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
            if self.velocities is not None:
                f.write("\nVelocities\n\n")
                np.savetxt(f, self.velocities, fmt=["%d", "%.16f", "%.16f", "%.16f"])

    def _make_header(self, out_file, n_atype):
        # cell = self.atoms.cell.cellpar()
        nat = len(self.atoms)
        s = "%s (written by toolbox by Jia-Xin Zhu)\n\n" % out_file
        s += "%d atoms\n" % nat
        s += "%d atom types\n" % n_atype
        if self.bonds is not None:
            s += "%d bonds\n" % len(self.bonds)
            s += "%d bond types\n" % len(np.unique(self.bonds[:, 1]))
        if self.angles is not None:
            s += "%d angles\n" % len(self.angles)
            s += "%d angle types\n" % len(np.unique(self.angles[:, 1]))
        if self.dihedrals is not None:
            s += "%d dihedrals\n" % len(self.dihedrals)
            s += "%d dihedral types\n" % len(np.unique(self.dihedrals[:, 1]))
        # s += "%.4f %.4f xlo xhi\n%.4f %.4f ylo yhi\n%.4f %.4f zlo zhi\n\n\n" % (
        #     0.0, cell[0], 0.0, cell[1], 0.0, cell[2])
        prismobj = Prism(self.atoms.get_cell())
        xhi, yhi, zhi, xy, xz, yz = convert(
            prismobj.get_lammps_prism(), "distance", "ASE", "metal"
        )
        s += "0.0 %.6f xlo xhi\n" % xhi
        s += "0.0 %.6f ylo yhi\n" % yhi
        s += "0.0 %.6f zlo zhi\n" % zhi
        if prismobj.is_skewed():
            s += "%.6f %.6f %.6f xy xz yz\n" % (xy, xz, yz)
        s += "\n"
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
            s += "%d %d %d %.16f %.16f %.16f %.16f\n" % (
                ii + 1,
                self.res_id[ii],
                self.atype[ii],
                self.charges[ii],
                self.positions[ii][0],
                self.positions[ii][1],
                self.positions[ii][2],
            )
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

    def set_dihedrals(self, dihedrals):
        self.dihedrals = np.reshape(dihedrals, (-1, 6))

    def set_charges(self, charges):
        self.charges = np.reshape(charges, (-1))

    def set_velocities(self, velocities):
        self.velocities = np.reshape(velocities, (-1, 4))

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


class DPLRLammpsData(LammpsData):
    """
    Example:

    atoms = io.read("coord.xyz")
    sel_type = ["O"]
    sys_charge_dict = {
        "O": 6.0,
        "H": 1.0,
    }
    lmp_data = DPLRLammpsData(atoms, sel_type=sel_type, sys_charge_dict=sys_charge_dict)
    lmp_data.write("system.data", format="lammps-data", specorder=["O", "H"], atom_style="full")
    """

    def __init__(self, atoms, sel_type, sys_charge_dict) -> None:
        atoms, center_ids = self._make_extended_atoms(
            atoms.copy(), sel_type, sys_charge_dict
        )
        super().__init__(atoms)
        # set bonds between real atoms and wannier atoms
        nbonds = len(center_ids)
        bonds = np.ones((nbonds, 4), dtype=int)
        np.copyto(bonds[:, 0], np.arange(nbonds) + 1)
        wannier_ids = np.arange(len(atoms) - nbonds, len(atoms)) + 1
        np.copyto(bonds[:, 2], center_ids)
        np.copyto(bonds[:, 3], wannier_ids)
        self.set_bonds(bonds)

    def _make_extended_atoms(self, atoms, sel_type, sys_charge_dict):
        dummy_type_map = ["He", "Ne", "Ar", "Kr", "Xe", "Rn"]
        # extend Wannier atoms
        center_ids = []
        for ii, _atype in enumerate(sel_type):
            sel_ids = np.where(atoms.symbols == _atype)[0]
            sel_atoms = atoms[sel_ids]
            while dummy_type_map[ii] in atoms.get_chemical_symbols():
                dummy_type_map.pop(ii)
            sel_atoms.symbols[:] = dummy_type_map[ii]
            atoms.extend(sel_atoms)
            center_ids.append(sel_ids + 1)
        center_ids = np.concatenate(center_ids)

        charges = np.full(len(atoms), -8.0)
        for k, v in sys_charge_dict.items():
            charges[atoms.symbols == k] = v
        atoms.set_initial_charges(charges)

        self.dummy_type_map = dummy_type_map[: len(sel_type)]
        return atoms, center_ids

    def write(self, out_file="system.data", specorder=None, **kwargs):
        if specorder is not None:
            specorder.extend(self.dummy_type_map)
        super().write(out_file, specorder=specorder, **kwargs)


class DPLRRestartLammpsData(LammpsData):
    """
    Example:

    atoms = io.read("after_system.data",
                format="lammps-data",
                sort_by_id=True,
                Z_of_type={1: 8, 2: 1, 3: 2}
               )
    sel_type = ["O"]
    lmp_data = DPLRRestartLammpsData(atoms, sel_type=sel_type)
    lmp_data.write("restart_system.data", format="lammps-data", specorder=["O", "H", "He"], atom_style="full")
    """

    def __init__(self, atoms, sel_type) -> None:
        n_wannier = np.count_nonzero(np.isin(atoms.symbols, sel_type))
        n_real = len(atoms) - n_wannier
        center_ids = []
        count = n_real
        coords = atoms.get_positions()
        for _atype in sel_type:
            sel_ids = np.where(atoms.symbols == _atype)[0]
            np.copyto(coords[count : count + len(sel_ids)], coords[sel_ids])
            center_ids.append(sel_ids + 1)
            count += len(sel_ids)
        center_ids = np.concatenate(center_ids)
        atoms.set_positions(coords)
        super().__init__(atoms)
        # set velocities
        # from ase internal unit to lammps metal unit (A/ps)
        # https://wiki.fysik.dtu.dk/ase/ase/units.html#units
        vs = atoms.get_velocities() * fs * 1e3
        # vs[n_real:] = 0.0
        self.set_velocities(
            np.concatenate((np.arange(1, len(atoms) + 1).reshape(-1, 1), vs), axis=1)
        )
        # set bonds
        nbonds = len(center_ids)
        bonds = np.ones((nbonds, 4), dtype=int)
        np.copyto(bonds[:, 0], np.arange(nbonds) + 1)
        wannier_ids = np.arange(len(atoms) - nbonds, len(atoms)) + 1
        np.copyto(bonds[:, 2], center_ids)
        np.copyto(bonds[:, 3], wannier_ids)
        self.set_bonds(bonds)


class OHHWaterLammpsData(LammpsData):
    def __init__(self, atoms) -> None:
        super().__init__(atoms)

        n_water = np.count_nonzero(atoms.symbols == "O")
        oxygen_ids = np.where(atoms.symbols == "O")[0] + 1
        hydrogen_ids = np.where(atoms.symbols == "H")[0] + 1

        nbonds = n_water * 2
        bonds = np.ones((nbonds, 4), dtype=int)
        np.copyto(bonds[:, 0], np.arange(nbonds) + 1)
        np.copyto(bonds[::2, 2], oxygen_ids)
        np.copyto(bonds[1::2, 2], oxygen_ids)
        np.copyto(bonds[:, 3], hydrogen_ids)

        angles = np.ones((n_water, 5), dtype=int)
        np.copyto(angles[:, 0], np.arange(n_water) + 1)
        np.copyto(angles[:, 2], hydrogen_ids[::2])
        np.copyto(angles[:, 3], oxygen_ids)
        np.copyto(angles[:, 4], hydrogen_ids[1::2])

        res_id = np.arange(n_water) + 1
        res_id = np.tile(res_id.reshape(-1, 1), [1, 3]).reshape(-1)

        self.set_bonds(bonds)
        self.set_angles(angles)
        self.set_res_id(res_id)


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
            with open(out_file, "a", encoding="utf-8") as f:
                header = self.make_header(atoms, ts)
                f.write(header)
                body = self.make_body(atoms, self.atype_dict)
                f.write(body)
        else:
            with open(out_file, "w", encoding="utf-8") as f:
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
        s += "ITEM: BOX BOUNDS %s %s %s\n" % (
            bc_dict[bc[0]],
            bc_dict[bc[1]],
            bc_dict[bc[2]],
        )
        s += "%.4f %.4f\n%.4f %.4f\n%.4f %.4f\n" % (
            0.0,
            cell[0],
            0.0,
            cell[1],
            0.0,
            cell[2],
        )
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
                s += "%d %d %s %.16f %.16f %.16f %.16f\n" % (
                    ii + 1,
                    atype_dict[atom.symbol]["type"],
                    atype_dict[atom.symbol]["element"],
                    ps[ii][0],
                    ps[ii][1],
                    ps[ii][2],
                    charges[ii],
                )
            else:
                s += "%d %d %s %.16f %.16f %.16f\n" % (
                    ii + 1,
                    atype_dict[atom.symbol]["type"],
                    atype_dict[atom.symbol]["element"],
                    ps[ii][0],
                    ps[ii][1],
                    ps[ii][2],
                )
        return s


class LammpsLog:
    def __init__(self, fname="log.lammps") -> None:
        self.log_file = fname
        with open(fname, "r") as f:
            self.content = f.readlines()
        # self.string = "".join(self.content)
        self.setup()

    def setup(self):
        for line in self.content:
            if re.search("MPI tasks", line):
                self.cpu_util = float(line.split()[0][:-1])
                self.n_mpi = int(line.split()[4])
                self.n_thread = int(line.split()[-3])
            if re.match("Loop time of", line):
                self.n_atoms = int(line.split()[-2])
                self.n_step = int(line.split()[-5])
                self.n_proc = int(line.split()[5])
                self.wall_time = float(line.split()[3])
            if re.match("Performance", line):
                out = line.split()
                self.performance = {
                    "ns_per_d": float(out[1]),
                    "h_per_ns": float(out[3]),
                    "ts_per_s": float(out[5]),
                }

    @property
    def timing_breakdown(self):
        start = False
        timing_breakdown = []
        for line in self.content:
            if start:
                timing_breakdown.append(line)
            if re.match("MPI task timing breakdown", line):
                start = True
            if start and re.match("Nlocal:", line):
                break

        ldata = [line.split("|") for line in timing_breakdown[2:-2]]
        # read data into a dict
        data = {
            d[0].strip().lower(): {
                "time": float(d[2].strip()),
                "percentage": float(d[-1].strip()),
            }
            for d in ldata
        }
        return data


def write_dump(traj, type_map, start=0, step=1, out_file="out.lammpstrj", append=False):
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
        with open(out_file, "a", encoding="utf-8") as f:
            header = make_dump_header(atoms, ts)
            f.write(header)
            body = make_dump_body(atoms, atype_dict)
            f.write(body)
    else:
        with open(out_file, "w", encoding="utf-8") as f:
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
    s += "%.4f %.4f\n%.4f %.4f\n%.4f %.4f\n" % (
        0.0,
        cell[0],
        0.0,
        cell[1],
        0.0,
        cell[2],
    )
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
            s += "%d %d %s %.16f %.16f %.16f %.16f\n" % (
                ii + 1,
                atype_dict[atom.symbol]["type"],
                atype_dict[atom.symbol]["element"],
                ps[ii][0],
                ps[ii][1],
                ps[ii][2],
                charges[ii],
            )
        else:
            s += "%d %d %s %.16f %.16f %.16f\n" % (
                ii + 1,
                atype_dict[atom.symbol]["type"],
                atype_dict[atom.symbol]["element"],
                ps[ii][0],
                ps[ii][1],
                ps[ii][2],
            )
    return s
