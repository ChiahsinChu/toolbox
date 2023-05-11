import logging
import os


class BashCalculator:
    def __init__(self, work_dir) -> None:
        self.work_dir = work_dir

    def run(self,
            command: str = None,
            stdin: str = None,
            stdout: str = "job.stdout",
            stderr: str = "job.stderr",
            mpi_command: str = "mpiexec.hydra"):
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

        logging.info("{:=^50}".format(" End calculation "))
        os.chdir(root_dir)


class LammpsCalculator(BashCalculator):
    def __init__(self, work_dir) -> None:
        super().__init__(work_dir)

    def run(self,
            command: str = "lmp",
            stdin: str = "input.lmp",
            stdout: str = "lammps.stdout",
            stderr: str = "lammps.stderr",
            mpi_command: str = "mpiexec.hydra",
            modifier: str = None):
        """
        https://docs.lammps.org/Run_options.html
        """
        stdin = "-i " + stdin
        if modifier is not None:
            stdin += " %s " % modifier
        super().run(command, stdin, stdout, stderr, mpi_command)