# SPDX-License-Identifier: LGPL-3.0-or-later
"""Matplotlib style module.

This module provides functionality for using custom matplotlib styles
and color maps for scientific plotting.

References
----------
- https://matplotlib.org/stable/tutorials/introductory/customizing.html
- https://matplotlib.org/stable/users/prev_whats_new/dflt_style_changes.html
- https://matplotlib.org/stable/gallery/style_sheets/style_sheets_reference.html
"""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import colors

MODULE_DIR = Path(__file__).resolve().parent


def use_style(style_name):
    """Use custom matplotlib style.
    
    Parameters
    ----------
    style_name : str
        Name of the style to use
        
    Notes
    -----
    If the specified style is not found, falls back to fivethirtyeight
    style and prints a warning.
    """
    plt.style.use("fivethirtyeight")
    fname = str(MODULE_DIR / (f"mplstyle/{style_name}.mplstyle"))
    # print(fname)
    try:
        plt.style.use(fname)
    except OSError:
        print(
            f"Warning: no style {style_name} is found. Use matplotlib default style."
        )
    """
    Set colors:
    https://stackoverflow.com/questions/68664116/is-there-a-way-to-change-the-color-names-color-in-mplstyle-file
    """
    color_map = colors.get_named_colors_mapping()
    color_map["red"] = "#D7422A"
    color_map["blue"] = "#003C88"
