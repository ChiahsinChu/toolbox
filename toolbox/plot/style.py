import matplotlib.pyplot as plt
from matplotlib import colors
from pathlib import Path

MODULE_DIR = Path(__file__).resolve().parent


def use_style(style_name):
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
