# SPDX-License-Identifier: LGPL-3.0-or-later
"""Deep Potential plotting module.

This module provides specialized plotting classes for visualizing
Deep Potential (DP) training data and results.
"""

from .figure import Figure


class DPTrainFigure(Figure):
    """Figure class for DP training visualization.
    
    This class extends Figure to provide specific
    functionality for visualizing Deep Potential training data.
    """
    
    def __init__(self, ax) -> None:
        """Initialize DPTrainFigure.
        
        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Matplotlib axes object
        """
        super().__init__(ax)

    def setup(self, x, y, **kwargs):
        """Set up DP training plot.
        
        Parameters
        ----------
        x : array_like
            Training data features
        y : array_like
            Training data targets
        **kwargs
            Additional keyword arguments
        """
        kwargs.update({"alpha": 0.2})
        self.ax.set_yscale("log")
        super().setup(x, y, **kwargs)
