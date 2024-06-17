# SPDX-License-Identifier: LGPL-3.0-or-later
from ..utils import *
from .figure import Figure


class DPTrainFigure(Figure):
    def __init__(self, ax) -> None:
        super().__init__(ax)

    def setup(self, x, y, **kwargs):
        kwargs.update({"alpha": 0.2})
        self.ax.set_yscale("log")
        super().setup(x, y, **kwargs)
        # def setup(self, x, y, **kwargs):
        #     self.ax.plot(x, y, **kwargs)
        #     self.ax.set_xlim(np.min(x), np.max(x))
