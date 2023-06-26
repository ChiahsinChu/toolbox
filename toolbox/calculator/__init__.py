from .cp2k import Cp2kCalculator
from .lammps import LammpsCalculator
from .dp import DPDispatcher, CP2KDPDispatcher

__all__ = ["Cp2kCalculator", "LammpsCalculator", "DPDispatcher", "CP2KDPDispatcher"]

def get_template_json(name):
    from . import template
    from ..utils.utils import save_dict
    d = getattr(template, name)
    save_dict(d, name + ".json")
   