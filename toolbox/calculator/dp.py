from dpdispatcher import Machine, Resources, Submission, Task
from toolbox.utils import *


class DPDispatcher:
    def __init__(self, work_dir=".") -> None:
        self.work_dir = work_dir
    
    def _setup(self, machine_setup={}, resources_setup={}):
        _machine_setup = {
            "batch_type": "Slurm",
            "context_type": "LocalContext",
            "local_root": "./",
            "remote_root": "/data/jxzhu/nnp/dp_workdir"
        }
        _machine_setup.update(machine_setup)
        _resources_setup = {
            "number_node": 1,
            "cpu_per_node": 24,
            "gpu_per_node": 0,
            "kwargs": {
                "gpu_usage": False
            },
            "queue_name": "c51-large",
            "group_size": 80,
            "module_purge": True
        }
        _resources_setup.update(resources_setup)
        self._check_node(_resources_setup)

        self.machine = Machine(**_machine_setup)
        self.resources = Resources(**_resources_setup)

    def run(self, dnames, machine_setup={}, resources_setup={}, task_setup={}):
        self._setup(machine_setup, resources_setup)

        task_list = []
        for dname in dnames:
            task = Task(task_work_path=dname, **task_setup)
            task_list.append(task)

        submission = Submission(work_base=self.work_dir,
                                machine=self.machine,
                                resources=self.resources,
                                task_list=task_list)
        submission.run_submission()

    @staticmethod
    def _check_node(resources_setup):
        queue = resources_setup["queue_name"] 

        if "c51" in queue:
            resources_setup["cpu_per_node"] = 24
        elif "c52" in queue:
            resources_setup["cpu_per_node"] = 28
        elif "c53" in queue:
            resources_setup["cpu_per_node"] = 32
        else:
            raise AttributeError("Unknown queue: %s" % queue)


#     def collect(self, type_map, out_dir="deepmd"):
#         safe_makedirs(os.path.join(self.work_dir, out_dir))

#         ms = dpdata.MultiSystems.from_dir(dir_name=os.path.join(
#                 self.work_dir, "iter.000000/02.fp"), file_name="output", fmt="cp2k/output")
#         ms.to_deepmd_npy(os.path.join(self.work_dir, "tmp_dpdata"))

#         fname_type_map = glob.glob("tmp_dpdata/**/type_map.raw", recursive=True)[0]
#         _type_map = np.loadtxt(fname_type_map, dtype=str)
#         dname_dpdata = os.path.dirname(fname_type_map)
#         _atype = np.loadtxt(os.path.join(dname_dpdata, "type.raw"), dtype=np.int32)

#         # npy data
#         shutil.move(os.path.join(dname_dpdata, "set.000"), os.path.join(out_dir, "set.000"))
#         # type_map.raw
#         np.savetxt(os.path.join(out_dir, "type_map.raw"), type_map, fmt="%s")
#         # type.raw
#         atype = [type_map.index(_type_map[ii]) for ii in _atype]
#         np.savetxt(os.path.join(out_dir, "type.raw"), atype, fmt="%d")
        
#         safe_makedirs(os.path.join(self.work_dir, "tmp_dpdata"))


class CP2KDPDispatcher(DPDispatcher):
    def __init__(self, work_dir=".", v9:bool=False) -> None:
        super().__init__(work_dir)
        self.v9 = v9

    def _setup(self, machine_setup={}, resources_setup={}):
        if self.v9:
            _resources_setup = {
                "custom_flags": [
                    "#SBATCH -J cp2k"
                ],
                "module_list": [
                    "mkl/latest",
                    "mpi/latest",
                    "gcc/9.3.0",
                    "cp2k/9.1"
                ]
            }
        else:
            _resources_setup = {
                "custom_flags": [
                    "#SBATCH -J cp2k"
                ],
                "module_list": [
                    "mpi/intel/2017.5.239",
                    "intel/17.5.239",
                    "gcc/7.4.0",
                    "cp2k/7.1"
                ]
            }
        _resources_setup.update(resources_setup)
        super()._setup(machine_setup, _resources_setup)

    def run(self, dnames, machine_setup={}, resources_setup={}, task_setup={}):
        _task_setup = {
            "command": "mpiexec.hydra cp2k.popt input.inp", 
            "forward_files": ["input.inp", "coord.xyz"], 
            "backward_files": [
                "cp2k-v_hartree-1_0.cube",
                "cp2k-ELECTRON_DENSITY-1_0.cube",
                "cp2k-TOTAL_DENSITY-1_0.cube", "cp2k-RESTART.wfn"
            ],
            "outlog":"output"
        }
        _task_setup.update(task_setup)
        super().run(dnames, machine_setup, resources_setup, _task_setup)