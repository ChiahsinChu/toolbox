params = {
    "type_map": ["O", "H", "Pt"],
    "mass_map": [15.999, 1.0079, 195.08],
    "init_data_prefix":
    "/data",
    "init_data_sys": ["init.000"],
    "init_batch_size": [1],
    "sys_configs": [["/data"]],
    "sys_batch_size": [1],
    "numb_models":
    4,
    "train_param":
    "input.json",
    "default_training_param": {
        "model": {
            "descriptor": {
                "type": "se_a",
                "_comment": "modify according your system",
                "sel": [68, 136, 64],
                "rcut_smth": 0.5,
                "rcut": 6.0,
                "neuron": [25, 50, 100],
                "resnet_dt": false,
                "axis_neuron": 16,
                "seed": 1
            },
            "fitting_net": {
                "neuron": [240, 240, 240],
                "resnet_dt": true,
                "seed": 1
            }
        },
        "learning_rate": {
            "type": "exp",
            "start_lr": 0.005,
            "decay_steps": 2000,
            "_comment": "last 20000 or 400000"
        },
        "loss": {
            "start_pref_e": 0.02,
            "limit_pref_e": 1,
            "start_pref_f": 1000,
            "limit_pref_f": 1,
            "start_pref_v": 0,
            "limit_pref_v": 0
        },
        "training": {
            "systems": [],
            "set_prefix": "set",
            "numb_steps": 200000,
            "batch_size": "auto",
            "seed": 1,
            "disp_file": "lcurve.out",
            "disp_freq": 100,
            "numb_test": 4,
            "save_freq": 1000
        }
    },
    "model_devi_dt":
    0.0005,
    "model_devi_skip":
    0,
    "model_devi_f_trust_lo":
    0.1,
    "model_devi_f_trust_hi":
    0.5,
    "model_devi_clean_traj":
    false,
    "model_devi_jobs": [{
        "temps": [330, 430, 530],
        "sys_idx": [1],
        "trj_freq": 10,
        "nsteps": 24000,
        "ensemble": "nvt",
        "_idx": 0
    }],
    "fp_style":
    "cp2k",
    "shuffle_poscar":
    false,
    "fp_task_max":
    100,
    "fp_task_min":
    5,
    "fp_pp_path":
    ".",
    "fp_pp_files": [],
    "fp_params": {}
}

machine = {
    "api_version":
    "1.0",
    "train": [{
        "command": "dp",
        "machine": {
            "batch_type": "LSF",
            "context_type": "LocalContext",
            "local_root": "./",
            "remote_root": "/data/jxzhu/nnp/dp_workdir"
        },
        "resources": {
            "number_node": 1,
            "cpu_per_node": 4,
            "gpu_per_node": 1,
            "queue_name": "gpu3",
            "group_size": 1,
            "kwargs": {
                "gpu_usage": true,
                "gpu_new_syntax": true,
                "gpu_exclusive": true
            },
            "custom_flags": ["#BSUB -J train", "#BSUB -W 24:00"],
            "module_list": ["deepmd/2.0"]
        }
    }],
    "model_devi": [{
        "command": "lmp_mpi",
        "machine": {
            "batch_type": "Slurm",
            "context_type": "LocalContext",
            "local_root": "./",
            "remote_root": "/data/jxzhu/nnp/dp_workdir"
        },
        "resources": {
            "number_node": 1,
            "cpu_per_node": 1,
            "gpu_per_node": 1,
            "queue_name": "gpu3",
            "group_size": 20,
            "custom_flags": ["#SBATCH -J md", "#SBATCH -t 48:00:00"],
            "module_list": ["deepmd/2.0"]
        }
    }],
    "fp": [{
        "command":
        "mpiexec.hydra cp2k.popt input.inp",
        "machine": {
            "batch_type": "Slurm",
            "context_type": "LocalContext",
            "local_root": "./",
            "remote_root": "/data/jxzhu/nnp/dp_workdir"
        },
        "resources": {
            "number_node": 1,
            "cpu_per_node": 32,
            "gpu_per_node": 0,
            "custom_flags": ["#SBATCH -J cp2k"],
            "queue_name": "c53-medium",
            "group_size": 120,
            "module_list":
            ["mkl/latest", "mpi/latest", "gcc/9.3.0", "cp2k/9.1"]
        },
        "user_backward_files": [
            "cp2k-v_hartree-1_0.cube", "cp2k-ELECTRON_DENSITY-1_0.cube",
            "cp2k-TOTAL_DENSITY-1_0.cube", "cp2k-RESTART.wfn"
        ]
    }]
}
