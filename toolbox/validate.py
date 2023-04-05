import os
import glob
import numpy as np

from deepmd.infer import DeepPot

class DPTest:
    def __init__(self, 
                 system, 
                 dp_model="graph.pb") -> None:
        self.system = system
        self.model = DeepPot(dp_model)
    
    def single_run(self, set_dir):
        coord = np.load()
        cell = np.load()
        hardness = np.load()

        return info_dict

    def run(self):
        atype = np.loadtxt("type.raw", dtype=np.int32)
        set_dirs = glob.glob(os.path.join(self.system, "set.*"))
        for set_dir in set_dirs:
            info_dict = self.single_run(set_dir=set_dir)
            self.print_info(info_dict)

    @ staticmethod
    def print_info(info_dict):
        pass

def validate(system, dp_model="graph.pb"):
    
    
        

# 
# e, f, v = dp.eval(coord, cell, atom_types)