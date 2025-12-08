# SPDX-License-Identifier: LGPL-3.0-or-later
import matplotlib.pyplot as plt
import numpy as np

from .core import ax_setlabel


class Figure:
    """Base class for matplotlib figure creation.
    
    This class provides a framework for creating
    matplotlib figures with consistent styling and layout.
    """
    
    def __init__(self, **kwargs) -> None:
        """Initialize Figure.
        
        Parameters
        ----------
        **kwargs
            Additional keyword arguments passed to plt.subplots
        """
        self.fig, self.ax = plt.subplots(**kwargs)

    def setup(self, x, y, xlim=None, ylim=None, **kwargs):
        """Set up basic plot with data.
        
        Parameters
        ----------
        x : array_like
            X data values
        y : array_like
            Y data values
        xlim : tuple, optional
            X-axis limits as (min, max), by default None
        ylim : tuple, optional
            Y-axis limits as (min, max), by default None
        **kwargs
            Additional keyword arguments passed to ax.plot
        """
        self.ax.plot(x, y, **kwargs)
        if xlim is None:
            xlim = (np.min(x), np.max(x))
        if ylim is None:
            ylim = (
                np.min(y) - (np.max(y) - np.min(y)) * 0.1,
                np.max(y) + (np.max(y) - np.min(y)) * 0.1,
            )
        self.xlim = xlim
        self.ylim = ylim
        self.ax.set_xlim(xlim)
        self.ax.set_ylim(ylim)

    def set_labels(self, kw, xlabel=None, ylabel=None, **kwargs):
        """Set axis labels from keyword dictionary.
        
        Parameters
        ----------
        kw : str
            Keyword to look up in label dictionary
        xlabel : str, optional
            X-axis label, by default None
        ylabel : str, optional
            Y-axis label, by default None
        **kwargs
            Additional keyword arguments passed to ax_setlabel
        """
        try:
            labels = label_dict.get(kw)
            xlabel = labels[0]
            ylabel = labels[1]
        except KeyError:
            assert (xlabel is not None) and (ylabel is not None)
            ax_setlabel(self.ax, xlabel, ylabel, **kwargs)


class FullCellFigure(Figure):
    """Figure class for full cell visualization.
    
    This class extends Figure to provide specific
    functionality for visualizing systems with full periodic
    boundary conditions.
    """
    
    def __init__(self, **kwargs) -> None:
        """Initialize FullCellFigure.
        
        Parameters
        ----------
        **kwargs
            Additional keyword arguments passed to parent class
        """
        super().__init__(**kwargs)

    def setup(self, x, y, xlim=None, ylim=None, z_surfs=None, **kwargs):
        super().setup(x, y, xlim, ylim, **kwargs)

        ax = self.ax
        ax.axvline(x=z_surfs[0])
        ax.axvline(x=z_surfs[1])
        ax.axvline(x=(z_surfs[0] + z_surfs[1]) / 2, ls="--")
        ax.fill_between(
            x,
            self.ylim[0],
            self.ylim[1],
            where=x <= z_surfs[0],
            facecolor="gray",
            alpha=0.5,
        )
        ax.fill_between(
            x,
            self.ylim[0],
            self.ylim[1],
            where=x >= z_surfs[1],
            facecolor="gray",
            alpha=0.5,
        )


class HalfCellFigure(Figure):
    """Figure class for half cell visualization.
    
    This class extends Figure to provide specific
    functionality for visualizing systems with half periodic
    boundary conditions (slab geometry).
    """
    
    def __init__(self, **kwargs) -> None:
        """Initialize HalfCellFigure.
        
        Parameters
        ----------
        **kwargs
            Additional keyword arguments passed to parent class
        """
        super().__init__(**kwargs)

    def setup(self, x, y, xlim=None, ylim=None, z_surf=None, **kwargs):
        super().setup(x, y, xlim, ylim, **kwargs)

        ax = self.ax
        ax.axvline(x=z_surf)
        ax.fill_between(
            x, self.ylim[0], self.ylim[1], where=x < z_surf, facecolor="gray", alpha=0.5
        )


label_dict = {
    "density": [r"z [Å]", r"$\rho [g/cm^3]$"],
    "orient_inveps": [r"z [Å]", r"$\epsilon_{ori,\perp}^{-1}$"],
    "orient_eps": [r"z [Å]", r"$\epsilon_{ori,\perp}$"],
    "elec_inveps": [r"z [Å]", r"$\epsilon_{elec,\perp}^{-1}$"],
    "elec_eps": [r"z [Å]", r"$\epsilon_{elec,\perp}$"],
    "hartree": [r"z [Å]", r"$V_H$ [eV]"],
    "rho": [r"z [Å]", r"$\rho$"],
    "rho_pol": [r"z [Å]", r"$\rho_{pol}$"],
    "polarization": [r"z [Å]", r"P [eA$^{-2}$"],
}
