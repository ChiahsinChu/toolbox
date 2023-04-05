import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from sklearn import metrics

from .style import use_style

use_style("pub")


def ax_setlabel(ax, xlabel, ylabel, **kwargs):
    # set axis label
    ax.set_xlabel(xlabel, **kwargs)
    ax.set_ylabel(ylabel, **kwargs)


def ax_rmse(ax, x, y):
    # scatter
    ax.scatter(x, y, color='steelblue', alpha=0.2)
    # ref line
    ref = np.arange(x.min(), x.max(), (x.max() - x.min()) / 100)
    ax.plot(ref, ref, color='firebrick', lw=1.5)


def plot_lcurve(fname, col, xlabel=None, ylabel=None, **kwargs):
    fig, ax = plt.subplots()
    data = np.loadtxt(fname)
    x = data[:, 0]
    ax.plot(x, data[:, col], **kwargs)
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(bottom=0.)
    if xlabel is not None and ylabel is not None:
        ax_setlabel(ax, xlabel, ylabel)

    return fig, ax, data


def plot_rmse(x, y, xlabel, ylabel, **kwargs):
    """
    plot scatter/ref line
    return rmse
    """
    x = np.array(x)
    y = np.array(y)

    rmse = np.sqrt(metrics.mean_squared_error(x, y))

    fig, ax = plt.subplots(figsize=[4, 4])
    ax_rmse(ax, x, y)
    ax_setlabel(ax, xlabel, ylabel, **kwargs)

    return fig, ax, rmse


def plot_bin_stats(x, y, xlabel, ylabel, bins=None):
    """
    plot scatter/ref line
    return rmse
    """
    x = np.array(x).flatten()
    y = np.array(y).flatten()
    if bins is None:
        bins = len(x) // 10

    fig, ax = plt.subplots(figsize=[6, 4], dpi=200)
    ax_bin_stats(ax, x, y, bins=bins)
    ax_setlabel(ax, xlabel, ylabel)
    ax.legend(loc='center left', bbox_to_anchor=(1.1, 0.5))
    return fig, ax


def ax_bin_stats(ax, x, y, bins):
    # ref line (y = 0)
    ax.axhline(y=0, color='gray')
    # mean
    bin_means, bin_edges, binnumber = stats.binned_statistic(x, y, bins=bins)
    bin_means[np.isnan(bin_means)] = 0.0
    bin_centers = bin_edges[1:] - (bin_edges[1] - bin_edges[0]) / 2
    ax.plot(bin_centers, bin_means, color='black', lw=1.5, label='mean')
    # max/min
    bin_maxs, bin_edges, binnumber = stats.binned_statistic(x,
                                                            y,
                                                            statistic='max',
                                                            bins=bins)
    bin_mins, bin_edges, binnumber = stats.binned_statistic(x,
                                                            y,
                                                            statistic='min',
                                                            bins=bins)
    bin_maxs[np.isnan(bin_maxs)] = 0.0
    bin_mins[np.isnan(bin_maxs)] = 0.0
    ax.fill_between(bin_centers,
                    bin_mins,
                    bin_maxs,
                    color='silver',
                    label='[min, max]')
    # std
    bin_stds, bin_edges, binnumber = stats.binned_statistic(x,
                                                            y,
                                                            statistic='std',
                                                            bins=bins)
    bin_stds[np.isnan(bin_stds)] = 0.0
    ax.fill_between(bin_centers,
                    bin_means - bin_stds,
                    bin_means + bin_stds,
                    color='gray',
                    label='[mean-std, mean+std]')


def plot_colormap_lines(xs, ys, legends, xlabel, ylabel, colormap='GnBu'):
    """
    TBC
    """
    fig, ax = plt.subplots(figsize=[6, 4], dpi=200)
    ax_colormap_lines(ax, xs, ys, legends, colormap)
    ax_setlabel(ax, xlabel, ylabel)
    ax.legend(loc='center left', bbox_to_anchor=(1.1, 0.5))
    return fig, ax


def ax_colormap_lines(ax,
                      xs,
                      ys,
                      labels,
                      scale=(0., 1.),
                      colormap="GnBu",
                      **kwargs):
    cm_scales = (np.array(labels) - scale[0]) / (scale[1] - scale[0])
    for x, y, label, cm_scale in zip(xs, ys, labels, cm_scales):
        ax.plot(x,
                y,
                color=plt.get_cmap(colormap)(cm_scale),
                label=label,
                **kwargs)
    ax.set_xlim(np.min(x), np.max(x))
