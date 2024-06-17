# SPDX-License-Identifier: LGPL-3.0-or-later
import logging
import os


class BashCalculator:
    def __init__(self, work_dir) -> None:
        self.work_dir = work_dir

    def run(
        self,
        command: str = None,
        stdin: str = None,
        stdout: str = "job.stdout",
        stderr: str = "job.stderr",
        mpi_command: str = "mpiexec.hydra",
        ignore_finished_tag=False,
    ):
        if (ignore_finished_tag == True) or (
            os.path.exists(os.path.join(self.work_dir, "finished_tag")) == False
        ):
            if mpi_command is not None:
                command = "%s %s" % (mpi_command, command)

            root_dir = os.getcwd()
            os.chdir(self.work_dir)
            logging.info("{:=^50}".format(" Start calculation "))

            logging.info("Path: %s" % os.getcwd())
            if stdin is not None:
                command += " %s " % stdin
            if stdout is not None:
                command += " 1> %s " % stdout
            if stderr is not None:
                command += " 2> %s " % stderr
            os.system(command=command)

            self._make_finished_tag(stdout)

            logging.info("{:=^50}".format(" End calculation "))
            os.chdir(root_dir)

    @staticmethod
    def _make_finished_tag(stdout):
        pass
