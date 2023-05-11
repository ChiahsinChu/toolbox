"""
Ref:
- https://userguide.mdanalysis.org/stable/formats/index.html
- https://manual.gromacs.org/current/reference-manual/file-formats.html#xtc
"""
import MDAnalysis as mda

topo = "interface.psf"

traj = "dump.lammpstrj"
u_lmp = mda.Universe(topo, traj, topology_format="PSF", format="LAMMPSDUMP")
traj = "dump.xtc"
u_xtc = mda.Universe(topo, traj, topology_format="PSF", format="XTC")

# ASE: atoms.get_chemical_symbols()
u_xtc.atoms.types
# ASE: atoms.get_positions()
u_xtc.atoms.positions
"""
ASE:
for atoms in traj:
    atoms.get_positions()
"""
for atoms in u_xtc.trajectory:
    atoms.positions
