from typing import IO, Any, Union, Dict

from ase.md import MDLogger as _MDLogger
from ase import Atoms


class MDLogger(_MDLogger):
    def __init__(
        self,
        dyn: Any,  # not fully annotated so far to avoid a circular import
        atoms: Atoms,
        logfile: Union[IO, str],
        fname: str,
        write_kwargs: Dict = {},
        **kwargs,
    ):
        super(MDLogger, self).__init__(dyn, atoms, logfile, **kwargs)
        self.fname = fname
        self.write_kwargs = write_kwargs

    def __call__(self):
        super().__call__()
        self.atoms.write(self.fname, append=True, **self.write_kwargs)
