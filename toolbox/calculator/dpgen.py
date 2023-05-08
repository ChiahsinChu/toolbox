import dpdata
import shutil

from toolbox.io.cp2k import Cp2kInput
from toolbox.utils import *
from toolbox.utils.utils import save_dict, safe_makedirs

from .template import machine, params


class FPTask:

    def __init__(self, work_dir=".") -> None:
        self.work_dir = work_dir

    def run(self, configs, **kwargs):
        if isinstance(configs, str):
            self.configs = io.read(configs, index=":")
        elif isinstance(configs, list):
            self.configs = []
            for config in configs:
                self.configs.append(io.read(config))
        else:
            raise AttributeError("Unknown type of configs: %s" % type(configs))

        dip_cor = kwargs.pop("dip_cor", False)
        for ii, atoms in enumerate(self.configs):
            cp2k_inp = Cp2kInput(atoms,
                                 hartree=True,
                                 eden=True,
                                 totden=True,
                                 dip_cor=dip_cor)
            cp2k_inp.write(os.path.join(
                self.work_dir, "iter.000000/02.fp/task.000.%06d" % ii),
                           save_dict=True)

        self._setup(**kwargs)

    def collect(self, type_map, out_dir="deepmd"):
        safe_makedirs(out_dir)

        ms = dpdata.MultiSystems.from_dir(dir_name=os.path.join(
                self.work_dir, "iter.000000/02.fp"), file_name="output", fmt="cp2k/output")
        ms.to_deepmd_npy("tmp_dpdata")

        fname_type_map = glob.glob("tmp_dpdata/**/type_map.raw", recursive=True)[0]
        _type_map = np.loadtxt(fname_type_map, dtype=str)
        dname_dpdata = os.path.dirname(fname_type_map)
        _atype = np.loadtxt(os.path.join(dname_dpdata, "type.raw"), dtype=np.int32)

        # npy data
        shutil.move(os.path.join(dname_dpdata, "set.000"), os.path.join(out_dir, "set.000"))
        # type_map.raw
        np.savetxt(type_map, os.path.join(out_dir, "type_map.raw"), fmt="%s")
        # type.raw
        atype = [type_map.index(_type_map[ii]) for ii in _atype]
        np.savetxt(atype, os.path.join(out_dir, "type.raw"), fmt="%d")
        
        os.removedirs("tmp_dpdata")


    def _setup(self, **kwargs):
        # record.dpgen
        a = np.zeros((7, 2))
        a[:, 1] = np.arange(7)
        np.savetxt(os.path.join(self.work_dir, "record.dpgen"), fmt="%d")
        # params.json
        save_dict(params, os.path.join(self.work_dir, "params.json"))
        # machine.json
        self._set_queue(machine, **kwargs)
        save_dict(machine, os.path.join(self.work_dir, "machine.json"))

    @staticmethod
    def _set_queue(machine, queue="c52-medium"):
        machine["fp"][0]["resources"]["queue_name"] = queue

        if "c51" in queue:
            machine["fp"][0]["resources"]["cpu_per_node"] = 24
        elif "c52" in queue:
            machine["fp"][0]["resources"]["cpu_per_node"] = 28
        elif "c53" in queue:
            machine["fp"][0]["resources"]["cpu_per_node"] = 32
        else:
            raise AttributeError("Unknown queue: %s" % queue)
