import numpy as np
from . import *


class Figure:
    def __init__(self, ax) -> None:
        self.ax = ax

    def setup(self, x, y, **kwargs):
        self.ax.plot(x, y, **kwargs)
        self.ax.set_xlim(np.min(x), np.max(x))
        # self.ax.set_ylim(bottom=0.)


class DensityFigure(Figure):
    def __init__(self, ax) -> None:
        super().__init__(ax)

    def setup(self, x, y, **kwargs):
        super().setup(x, y, **kwargs)
        xlabel = r"$\rho [g/cm^3]$"
        ylabel = r"z [A]"
        ax_setlabel(self.ax, xlabel, ylabel)