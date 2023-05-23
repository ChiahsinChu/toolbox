# CP2K

> For pre-/post-processing of CP2K, you might also be interested in [pycp2k](https://github.com/SINGROUP/pycp2k) and [cp2kdata](https://github.com/robinzyb/cp2kdata).

## Generate input files: `Cp2kInput`

Create cp2k input files.

```python
from toolbox.utils import *
from toolbox.io.cp2k import Cp2kInput

atoms = io.read("coord.xyz")
# remember to set cell parameters and pbc!
cp2k_inp = Cp2kInput(atoms=atoms, input_type="energy", cutoff=1000)
fp_params = {
    "FORCE_EVAL": {
        "DFT": {
            "SCF": {
                "MAX_SCF": 3,
                "EPS_DIIS": 1e-15,
                "DIAGONALIZATION": {
                    "_": ".FALSE."
                },
                "MIXING": {
                    "ALPHA": 0.,
                    "METHOD": "DIRECT_P_MIXING"
                }
            }
        }
    }
}
cp2k_inp.write(output_dir="test", fp_params=fp_params)

```

You can modify the internal template by setting the variables when initializing the `Cp2kInput` or via the dict `fp_params`.

## Run CP2K calculation in python

You can perform CP2K calculations in python script via [`Cp2kCalculator`](calculator.md).

## Post-processing of output files

### `Cp2kOutput`

Read data from cp2k standard output (for single frame).
All output values have been converted to eV (from Hartree).

```python
from toolbox.io.cp2k import Cp2kOutput

cp2k_out = Cp2kOutput("output.out")
print(cp2k_out.fermi)
print(cp2k_out.energy)

```

### `Cp2kCube` and `Cp2kHartreeCube`

Read and do the average of cp2k cube file.
In `Cp2kCube`, no unit conversion is taken. In `Cp2kHartreeCube`, atmoic unit of energy (i.e. Hartree) is converted to electron volt (eV).

```python
from toolbox.utils import *
from toolbox.io.cp2k import Cp2kCube, Cp2kHartreeCube

e_cube = Cp2kCube("cp2k-TOTAL_DENSITY-1_0.cube")
output = e_cube.get_ave_cube(axis=2, gaussian_sigma=0.)
plt.plot(output[0], output[1])

v_cube = Cp2kHartreeCube("cp2k-v_hartree-1_0.cube")
output = v_cube.get_ave_cube(axis=2, gaussian_sigma=0.)
plt.plot(output[0], output[1])

plt.show()
```
