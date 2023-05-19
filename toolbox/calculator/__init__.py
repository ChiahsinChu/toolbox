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
            mpi_command: str = "mpiexec.hydra", 
            ignore_finished_tag=False):
        if (ignore_finished_tag == True) or (os.path.exists(os.path.join(self.work_dir, "finished_tag")) == False):
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


class LammpsCalculator(BashCalculator):
    def __init__(self, work_dir) -> None:
        super().__init__(work_dir)

    def run(self,
            command: str = "lmp",
            stdin: str = "input.lmp",
            stdout: str = "lammps.stdout",
            stderr: str = "lammps.stderr",
            mpi_command: str = "mpiexec.hydra",
            ignore_finished_tag=False,
            modifier: str = None):
        """
        https://docs.lammps.org/Run_options.html
        """
        stdin = "-i " + stdin
        if modifier is not None:
            stdin += " %s " % modifier
        super().run(command, stdin, stdout, stderr, mpi_command, ignore_finished_tag)

    @staticmethod
    def _make_finished_tag(stdout):
        # TODO: check if lammps task finished successfully 
        pass
