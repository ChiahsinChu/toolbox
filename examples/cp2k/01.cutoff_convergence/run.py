# SPDX-License-Identifier: LGPL-3.0-or-later
"""
Find a cutoff value for converged total energy
"""

import os

import matplotlib.pyplot as plt
import numpy as np
from ase import build

from toolbox import plot
from toolbox.calculator.cp2k import Cp2kCalculator
from toolbox.io.cp2k import Cp2kInput, Cp2kOutput

E_CONVERGENCE = 1e-2
plot.use_style("pub")

# prepare input files for CP2K calculation
atoms = build.fcc111("Al", (2, 2, 6), vacuum=10.0, orthogonal=True, periodic=True)

energy = []
_cutoff = []
for cutoff in np.arange(400, 1300, 100):
    dname = "%d" % cutoff
    cp2k_inp = Cp2kInput(atoms, cutoff=cutoff)
    cp2k_inp.write(dname)
    Cp2kCalculator(work_dir=dname).run(stdout="output")
    # total energy in eV
    e = Cp2kOutput(os.path.join(dname, "output")).energy[0]
    if len(energy) > 0:
        if np.abs(e - energy[-1]) <= E_CONVERGENCE:
            energy.append(e)
            _cutoff.append(cutoff)
            break
    energy.append(e)
    _cutoff.append(cutoff)

# make plot
fig, ax = plt.subplots()
xlabel = "Cutoff [Ry]"
ylabel = "Total energy [eV]"
ax.plot(_cutoff, energy, "-o", color="blue")
plot.ax_setlabel(ax, xlabel, ylabel)
fig.savefig("cutoff_test.png")
