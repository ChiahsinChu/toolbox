# Toolbox

## Introduction

Some python codes used for my project. Please feel free to reach me via email: **jiaxinzhu@stu.xmu.edu.cn** .

## Documentation

📚 **Live Documentation**: [https://ChiahsinChu.github.io/toolbox](https://ChiahsinChu.github.io/toolbox)

The documentation is automatically built and deployed to GitHub Pages whenever changes are pushed to the main branch.

## To-do list

- [ ] MultiFrameCp2kOutput
- [ ] [Unit Testing](doc/unittest.md)
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

## Documentation Development

This project uses MkDocs for documentation generation. The documentation is automatically built and deployed using GitHub Actions.

### Local Development

To build and serve the documentation locally:

```bash
# Install documentation dependencies
pip install -e ".[docs]"

# Build the documentation
mkdocs build

# Serve the documentation locally (for development)
mkdocs serve
```

The documentation will be available at `http://127.0.0.1:8000`.

### Documentation Structure

- `docs/` - Contains all documentation source files
- `mkdocs.yml` - MkDocs configuration file
- `.github/workflows/docs.yml` - GitHub Actions workflow for auto-deployment
