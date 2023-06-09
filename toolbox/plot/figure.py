from ..utils import *
from .core import *


class Figure:
    def __init__(self, **kwargs) -> None:
        self.fig, self.ax = plt.subplots(**kwargs)

    def setup(self, x, y, **kwargs):
        self.ax.plot(x, y, **kwargs)

    def set_labels(self, kw, xlabel=None, ylabel=None, **kwargs):
        try:
            labels = label_dict.get(kw)
            xlabel = labels[0]
            ylabel = labels[1]
        except:
            assert (xlabel is not None) and (ylabel is not None)   
        ax_setlabel(self.ax, xlabel, ylabel, **kwargs)


class FullCellFigure(Figure):
    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)

    def setup(self, x, y, xlim=None, ylim=None, z_surfs=None, **kwargs):
        if xlim is None:
            xlim = (np.min(x), np.max(x))
        if ylim is None:
            ylim = (np.min(y), np.max(y))
        self.ax.set_xlim(xlim)
        self.ax.set_ylim(ylim)

        super().setup(x, y, **kwargs)
        ax = self.ax
        ax.axvline(x=z_surfs[0])
        ax.axvline(x=z_surfs[1])
        ax.axvline(x=(z_surfs[0] + z_surfs[1]) / 2, ls="--")
        ax.fill_between(x, ylim[0], ylim[1], where=x <= z_surfs[0], facecolor='gray', alpha=.5)
        ax.fill_between(x, ylim[0], ylim[1], where=x >= z_surfs[1], facecolor='gray', alpha=.5)

class HalfCellFigure(Figure):
    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)

    def setup(self, x, y, z_surf=None, **kwargs):
        if xlim is None:
            xlim = (np.min(x), np.max(x))
        if ylim is None:
            ylim = (np.min(y), np.max(y))
        self.ax.set_xlim(xlim)
        self.ax.set_ylim(ylim)

        super().setup(x, y, xlim, ylim, **kwargs)
        ax = self.ax
        ax.axvline(x=z_surf)
        ax.fill_between(x, ylim[0], ylim[1], where=x < z_surf, facecolor='gray', alpha=.5)
        

label_dict = {
    "density": [r"z [A]", r"$\rho [g/cm^3]$"],
    "orient_inveps": [r"z [A]", r"$\epsilon_{ori,\perp}^{-1}$"],
    "orient_eps": [r"z [A]", r"$\epsilon_{ori,\perp}$"],
    "elec_inveps": [r"z [A]", r"$\epsilon_{elec,\perp}^{-1}$"],
    "elec_eps": [r"z [A]", r"$\epsilon_{elec,\perp}$"],
    "hartree": [r"z [A]", r"$V_H$ [eV]"],
    "rho": [r"z [A]", r"$\rho$"]
}