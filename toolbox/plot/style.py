import matplotlib.pyplot as plt
from pathlib import Path

MODULE_DIR = Path(__file__).resolve().parent

def use_style(style_name):
    fname = str(MODULE_DIR / ("mplstyle/%s.mplstyle" % style_name))
    # print(fname)
    try:
        plt.style.use(fname)
    except:
        print("Warning: no style %s is found. Use matplotlib default style." % style_name)