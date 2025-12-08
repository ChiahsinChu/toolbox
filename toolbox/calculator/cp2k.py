# SPDX-License-Identifier: LGPL-3.0-or-later
"""CP2K quantum chemistry calculator interface.

This module provides a calculator interface for running CP2K quantum
chemistry calculations with proper error handling and output checking.
"""

import logging
import os
import sys

from ..io.cp2k import Cp2kOutput
from .bash import BashCalculator


class Cp2kCalculator(BashCalculator):
    """A calculator for running CP2K quantum chemistry calculations.

    This class extends BashCalculator to provide specific functionality for
    CP2K calculations, including error handling and output checking.
    """

    def __init__(self, work_dir) -> None:
        """Initialize the Cp2kCalculator.

        Parameters
        ----------
        work_dir : str
            The working directory where CP2K calculations will be executed
        """
        super().__init__(work_dir)

    def run(
        self,
        command="cp2k.popt",
        stdin="input.inp",
        stdout="output.out",
        stderr="cp2k.stderr",
        mpi_command: str = "mpiexec.hydra",
        ignore_finished_tag: bool = False,
        ignore_err: bool = True,
    ):
        """Run a CP2K calculation.

        Parameters
        ----------
        command : str, optional
            CP2K executable command, by default "cp2k.popt"
        stdin : str, optional
            Input file name, by default "input.inp"
        stdout : str, optional
            Standard output file name, by default "output.out"
        stderr : str, optional
            Standard error file name, by default "cp2k.stderr"
        mpi_command : str, optional
            MPI command to use, by default "mpiexec.hydra"
        ignore_finished_tag : bool, optional
            Whether to ignore existing finished tag, by default False
        ignore_err : bool, optional
            Whether to ignore calculation errors, by default True
        """
        self.ignore_err = ignore_err
        super().run(command, stdin, stdout, stderr, mpi_command, ignore_finished_tag)

    def _make_finished_tag(self, stdout):
        """Create a finished tag after successful CP2K calculation.

        Parameters
        ----------
        stdout : str
            Name of the CP2K output file to check for completion

        Raises
        ------
        SystemExit
            If CP2K calculation did not finish and ignore_err is False
        """
        try:
            Cp2kOutput(stdout)
            open(os.path.join("finished_tag"), "w").close()
        except Exception:
            warning_msg = "CP2K calculation does not finish!"
            if self.ignore_err:
                logging.warning(warning_msg)
            else:
                sys.exit(warning_msg)


# from ase import io
# from ase.calculators.cp2k import CP2K
# from ..utils.utils import iterdict, load_dict

# class DeprecatedCp2kCalculator:
#     def __init__(self, work_dir) -> None:
#         self.work_dir = work_dir

#     def run(self, type="bash", ignore_err=False, **kwargs):
#         self.ignore_err_tag = ignore_err
#         root_dir = os.getcwd()
#         os.chdir(self.work_dir)
#         logging.info("{:=^50}".format(" Start: CP2K calculation "))
#         getattr(self, "run_%s" % type)(**kwargs)
#         logging.info("{:=^50}".format(" End: CP2K calculation "))
#         os.chdir(root_dir)

#     def run_bash(self,
#                  command="mpiexec.hydra cp2k.popt",
#                  stdin="input.inp",
#                  stdout="output.out",
#                  stderr="cp2k.stderr"):
#         logging.info("Path: %s" % os.getcwd())
#         if stdin is not None:
#             command += " %s " % stdin
#         if stdout is not None:
#             command += " 1> %s " % stdout
#         if stderr is not None:
#             command += " 2> %s " % stderr
#         os.system(command=command)

#         # args = shlex.split(command)
#         # args.append(stdin)
#         # f = open(stdout, "w")
#         # Popen(args, stdout=f)
#         # f.close()

#         try:
#             cp2k_out = Cp2kOutput(stdout)
#             with open("finished_tag", 'w') as f:
#                 pass
#         except:
#             if self.ignore_err_tag:
#                 logging.warning("Calculation does not finish!")
#             else:
#                 sys.exit("Calculation does not finish!")

#     def run_ase(self):
#         atoms = io.read("coord.xyz")
#         inp_dict = load_dict("input.json")
#         label = inp_dict["GLOBAL"].pop("PROJECT", "cp2k")
#         inp_dict["FORCE_EVAL"]["SUBSYS"].pop("TOPOLOGY", None)
#         inp = self._dict_to_cp2k_input(inp_dict)

#         calc = CP2K(command=self.command,
#                     inp=inp,
#                     label=label,
#                     force_eval_method=None,
#                     print_level=None,
#                     stress_tensor=None,
#                     basis_set=None,
#                     pseudo_potential=None,
#                     basis_set_file=None,
#                     potential_file=None,
#                     cutoff=None,
#                     max_scf=None,
#                     xc=None,
#                     uks=None,
#                     charge=None,
#                     poisson_solver=None)
#         atoms.calc = calc

#         atoms.get_potential_energy()

#     @staticmethod
#     def _dict_to_cp2k_input(input_dict):
#         input_str = iterdict(input_dict, out_list=["\n"], loop_idx=0)
#         s = "\n".join(input_str)
#         s = s.strip("\n")
#         return s
