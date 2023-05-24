# Toolbox

## Introduction

Some python codes used for my project. Please feel free to reach me via email: **jiaxinzhu@stu.xmu.edu.cn** .

## To-do list

- [ ] MultiFrameCp2kOutput
- [ ] [unittest](doc/unittest.md)
- [ ] write class for workflow plotting

## Installation

```bash
git clone --recursive https://github.com/ChiahsinChu/toolbox.git
cd toolbox
pip install .
```

If you don't use `--recursive`, you will miss the submodule(s). You can, alternatively, download the submodule(s) in two steps:

```bash
git clone https://github.com/ChiahsinChu/toolbox.git
cd toolbox
git submodule update --init --recursive
```

You can check if the installation succeeds by running:

```bash
python -c "import toolbox; toolbox.test()"
```

## User's guide

### pre-/post-processing in commonly used packages

- [CP2K](doc/cp2k.md)
- [LAMMPS](doc/lammps.md)
- [DP](doc/dp.md)

### routines in projects

- [do calculation in python scripts](doc/calculator.md)
- [make plots](doc/plot.md)
