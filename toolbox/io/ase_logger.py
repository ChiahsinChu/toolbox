from typing import IO, Any, Union, Optional

from ase.md import MDLogger as _MDLogger
from ase import Atoms


class MDLogger(_MDLogger):
    def __init__(
        self,
        dyn: Any,  # not fully annotated so far to avoid a circular import
        atoms: Atoms,
        logfile: Union[IO, str],
        fname: str,
        fmt: Optional[str] = None,
        **kwargs,
    ):
        super(MDLogger, self).__init__(dyn, atoms, logfile, **kwargs)
        self.fname = fname
        self.fmt = fmt

    def __call__(self):
        super().__call__()
        self.atoms.write(self.fname, format=self.fmt, append=True)
