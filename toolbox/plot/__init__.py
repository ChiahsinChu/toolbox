# SPDX-License-Identifier: LGPL-3.0-or-later
"""Plotting and visualization utilities for computational chemistry.

This module provides tools for:
- Creating publication-quality figures
- Customizing plot styles
- Visualizing training data and results
- Statistical analysis plots

Examples
--------
>>> from toolbox.plot import use_style
>>> use_style("pub")  # Use publication style
"""
from .core import (
    ax_bin_stats,
    ax_colormap_lines,
    ax_rmse,
    ax_setlabel,
    plot_bin_stats,
    plot_colormap_lines,
    plot_lcurve,
    plot_rmse,
)
from .style import use_style

__all__ = [
    "use_style",
    "ax_setlabel",
    "ax_rmse",
    "ax_bin_stats",
    "ax_colormap_lines",
    "plot_lcurve",
    "plot_rmse",
    "plot_bin_stats",
    "plot_colormap_lines",
]
