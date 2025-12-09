# Toolbox

## Introduction

Some python codes used for my project. Please feel free to reach me via email: **jiaxinzhu@stu.xmu.edu.cn** .

## To-do list

- [ ] MultiFrameCp2kOutput
- [ ] [Unit Testing](unittest.md)
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

- [CP2K](cp2k.md)
- [LAMMPS](lammps.md)
- [Deep Potential](dp.md)

### routines in projects

- [do calculation in python scripts](calculator.md)
- [make plots](plot.md)