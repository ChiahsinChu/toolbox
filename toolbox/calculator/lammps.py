# SPDX-License-Identifier: LGPL-3.0-or-later
"""LAMMPS calculator module.

This module provides a calculator for running LAMMPS molecular dynamics
simulations, extending the BashCalculator class.
"""

import logging
import os
import re

from .bash import BashCalculator


class LammpsCalculator(BashCalculator):
    """A calculator for running LAMMPS molecular dynamics simulations.
    
    This class extends BashCalculator to provide specific functionality for
    LAMMPS calculations, including completion checking.
    """
    
    def __init__(self, work_dir) -> None:
        """Initialize LammpsCalculator.
        
        Parameters
        ----------
        work_dir : str
            The working directory where LAMMPS calculations will be executed
        """
        super().__init__(work_dir)

    def run(
        self,
        command: str = "lmp",
        stdin: str = "input.lmp",
        stdout: str = "lammps.stdout",
        stderr: str = "lammps.stderr",
        mpi_command: str = "mpiexec.hydra",
        ignore_finished_tag=False,
        modifier: str = None,
    ):
        """Run a LAMMPS calculation.
        
        Parameters
        ----------
        command : str, optional
            LAMMPS executable command, by default "lmp"
        stdin : str, optional
            Input file name, by default "input.lmp"
        stdout : str, optional
            Standard output file name, by default "lammps.stdout"
        stderr : str, optional
            Standard error file name, by default "lammps.stderr"
        mpi_command : str, optional
            MPI command to use, by default "mpiexec.hydra"
        ignore_finished_tag : bool, optional
            Whether to ignore existing finished tag, by default False
        modifier : str, optional
            Additional command modifiers, by default None
            
        Notes
        -----
        For more information on LAMMPS run options, see:
        https://docs.lammps.org/Run_options.html
        """
        stdin = "-i " + stdin
        if modifier is not None:
            stdin += f" {modifier} "
        super().run(command, stdin, stdout, stderr, mpi_command, ignore_finished_tag)

    @staticmethod
    def _make_finished_tag(stdout):
        """Create a finished tag after successful LAMMPS calculation.
        
        Parameters
        ----------
        stdout : str
            Name of LAMMPS output file to check for completion
        """
        with open(stdout, "rb") as f:
            offset = -50
            while True:
                f.seek(offset, 2)
                lines = f.readlines()
                if len(lines) >= 2:
                    last_line = lines[-1]
                    break
                offset *= 2

        pattern = re.compile(r"Total wall time")
        if pattern.search(last_line.decode()) is not None:
            with open(os.path.join("finished_tag"), "w") as f:
                pass
        else:
            warning_msg = "LAMMPS calculation does not finish!"
            logging.warning(warning_msg)
