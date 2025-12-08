# SPDX-License-Identifier: LGPL-3.0-or-later
"""Bash calculator for running external commands.

This module provides a calculator interface for running bash commands
and external programs with proper error handling and logging.
"""

import logging
import os


class BashCalculator:
    """A calculator for running bash commands in a specified working directory.
    
    This class provides functionality to execute shell commands with proper
    directory management and logging.
    """
    
    def __init__(self, work_dir) -> None:
        """Initialize the BashCalculator.
        
        Parameters
        ----------
        work_dir : str
            The working directory where commands will be executed
        """
        self.work_dir = work_dir

    def run(
        self,
        command: str = None,
        stdin: str = None,
        stdout: str = "job.stdout",
        stderr: str = "job.stderr",
        mpi_command: str = "mpiexec.hydra",
        ignore_finished_tag: bool = False,
    ):
        """Execute a command in the working directory.
        
        Parameters
        ----------
        command : str, optional
            Command to execute, by default None
        stdin : str, optional
            Standard input file, by default None
        stdout : str, optional
            Standard output file, by default "job.stdout"
        stderr : str, optional
            Standard error file, by default "job.stderr"
        mpi_command : str, optional
            MPI command to use, by default "mpiexec.hydra"
        ignore_finished_tag : bool, optional
            Whether to ignore existing finished tag, by default False
        """
        if (ignore_finished_tag) or (
            not os.path.exists(os.path.join(self.work_dir, "finished_tag"))
        ):
            if mpi_command is not None:
                command = f"{mpi_command} {command}"

            root_dir = os.getcwd()
            os.chdir(self.work_dir)
            logging.info("{:=^50}".format(" Start calculation "))

            logging.info(f"Path: {os.getcwd()}")
            if stdin is not None:
                command += f" {stdin} "
            if stdout is not None:
                command += f" 1> {stdout} "
            if stderr is not None:
                command += f" 2> {stderr} "
            os.system(command=command)

            self._make_finished_tag(stdout)

            logging.info("{:=^50}".format(" End calculation "))
            os.chdir(root_dir)

    @staticmethod
    def _make_finished_tag(stdout):
        """Create a finished tag file after successful execution.
        
        Parameters
        ----------
        stdout : str
            Name of the standard output file to check
        """
        pass
