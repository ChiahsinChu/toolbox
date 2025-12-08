# SPDX-License-Identifier: LGPL-3.0-or-later
"""Template configuration for dpdispatcher.

This module contains default machine and resource configurations
for running calculations on HPC systems.
"""

machine = {
    "batch_type": "Slurm",
    "context_type": "LocalContext",
    "local_root": "./",
    "remote_root": "/data/jxzhu/nnp/dp_workdir",
}

resources = {
    "number_node": 1,
    "cpu_per_node": 24,
    "gpu_per_node": 0,
    "kwargs": {"gpu_usage": False},
    "custom_flags": ["#SBATCH -J cp2k", "#SBATCH -t 24:00:00"],
    "queue_name": "c51-large",
    "group_size": 80,
    "module_list": ["mpi/intel/2017.5.239", "intel/17.5.239", "gcc/7.4.0", "cp2k/7.1"],
}
