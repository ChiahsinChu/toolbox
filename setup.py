from setuptools import setup, find_packages

setup(
    name="toolbox",
    version="1.0",
    # include all packages in src
    packages=find_packages(),
    # exclude by
    #find_packages(exclude=["*.tests", "*.tests.*", "tests.*", "tests"])
    # all packages required for this package is in src
    #package_dir = {'':'zjxpack'},
    include_package_data=True,
    install_requires=[
        "monty", "regex", "ase", "scipy", "scikit-learn", "h5py"
    ])
