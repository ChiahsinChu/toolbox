# SPDX-License-Identifier: LGPL-3.0-or-later
from ase import build

from toolbox.calculator.cp2k import Cp2kCalculator
from toolbox.io.cp2k import Cp2kHartreeCube, Cp2kInput, Cp2kOutput

# prepare input files for CP2K calculation
atoms = build.fcc111("Al", (2, 2, 6), vacuum=10.0, orthogonal=True, periodic=True)
cp2k_inp = Cp2kInput(atoms, cutoff=200, hartree=True)
cp2k_inp.write(".")

# run CP2K for single point energy calculation
Cp2kCalculator(work_dir=".").run(stdout="output")

cube = Cp2kHartreeCube("cp2k-v_hartree-1_0.cube")
out = cube.get_ave_cube()
# get the Hartree potential in vacuum [eV]
ref_hartree = out[1][-1]
# get the Fermi level [eV]
e_fermi = Cp2kOutput("output").fermi

print("Workfunction: %.4f [eV]" % (e_fermi - ref_hartree))
