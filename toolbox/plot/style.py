"""
References:
- https://matplotlib.org/stable/tutorials/introductory/customizing.html
- https://matplotlib.org/stable/users/prev_whats_new/dflt_style_changes.html
- https://matplotlib.org/stable/gallery/style_sheets/style_sheets_reference.html
"""
import matplotlib.pyplot as plt
from matplotlib import colors
from pathlib import Path

MODULE_DIR = Path(__file__).resolve().parent


def use_style(style_name):
    plt.style.use("fivethirtyeight")
    fname = str(MODULE_DIR / ("mplstyle/%s.mplstyle" % style_name))
    # print(fname)
    try:
        plt.style.use(fname)
    except:
        print("Warning: no style %s is found. Use matplotlib default style." %
              style_name)
    """
    Set colors:
    https://stackoverflow.com/questions/68664116/is-there-a-way-to-change-the-color-names-color-in-mplstyle-file
    """
    color_map = colors.get_named_colors_mapping()
    color_map["red"] = "#D7422A"
    color_map["blue"] = "#003C88"
