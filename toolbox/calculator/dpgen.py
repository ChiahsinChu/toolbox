from toolbox.utils import *
from toolbox.utils.utils import save_dict
from toolbox.io.cp2k import Cp2kInput

from .template import machine, params


class FPTask:
    def __init__(self, configs) -> None:

        if isinstance(configs, str):
            self.configs = io.read(configs, index=":")
        elif isinstance(configs, list):
            self.configs = []
            for config in configs:
                self.configs.append(io.read(config))
        else:
            raise AttributeError("Unknown type of configs: %s" % type(configs))

    def run(self, dst, **kwargs):
        dip_cor = kwargs.pop("dip_cor", False)
        for ii, atoms in enumerate(self.configs):
            cp2k_inp = Cp2kInput(atoms,
                                 hartree=True,
                                 eden=True,
                                 totden=True,
                                 dip_cor=dip_cor)
            cp2k_inp.write(os.path.join(
                dst, "iter.000000/02.fp/task.000.%06d" % ii),
                           save_dict=True)

        self._setup(dst, **kwargs)

    def _setup(self, dname, **kwargs):
        # record.dpgen
        a = np.zeros((7, 2))
        a[:, 1] = np.arange(7)
        np.savetxt(os.path.join(dname, "record.dpgen"), fmt="%d")
        # params.json
        save_dict(params, os.path.join(dname, "params.json"))
        # machine.json
        self._set_queue(machine, **kwargs)
        save_dict(machine, os.path.join(dname, "machine.json"))

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
