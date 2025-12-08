# SPDX-License-Identifier: LGPL-3.0-or-later
"""Core plotting utilities module.

This module provides core plotting functions and utilities for creating
various types of plots including learning curves, RMSE plots,
binned statistics, and colored line plots.
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from sklearn import metrics


def ax_setlabel(ax, xlabel, ylabel, **kwargs):
    """Set axis labels for matplotlib axes.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Matplotlib axes object
    xlabel : str
        X-axis label
    ylabel : str
        Y-axis label
    **kwargs
        Additional keyword arguments passed to set_xlabel/set_ylabel
    """
    ax.set_xlabel(xlabel, **kwargs)
    ax.set_ylabel(ylabel, **kwargs)


def ax_rmse(ax, x, y):
    """Create RMSE scatter plot with reference line.
    
    This function creates a scatter plot of data points
    with a reference line (y=x) for visual comparison.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Matplotlib axes object
    x : array_like
        X values
    y : array_like
        Y values
    """
    # scatter
    ax.scatter(x, y, color="steelblue", alpha=0.2)
    # ref line
    ref = np.arange(x.min(), x.max(), (x.max() - x.min()) / 100)
    ax.plot(ref, ref, color="firebrick", lw=1.5)


def plot_lcurve(fname, col, xlabel=None, ylabel=None, **kwargs):
    """Plot learning curve from data file.
    
    This function reads a data file and plots a learning curve
    from the specified column.
    
    Parameters
    ----------
    fname : str
        Path to data file
    col : int
        Column index to plot (0-based)
    xlabel : str, optional
        X-axis label, by default None
    ylabel : str, optional
        Y-axis label, by default None
    **kwargs
        Additional keyword arguments passed to ax.plot()
        
    Returns
    -------
    tuple
        Tuple of (fig, ax, data) containing matplotlib figure,
        axes object, and loaded data
    """
    fig, ax = plt.subplots()
    data = np.loadtxt(fname)
    x = data[:, 0]
    ax.plot(x, data[:, col], **kwargs)
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(bottom=0.0)
    if xlabel is not None and ylabel is not None:
        ax_setlabel(ax, xlabel, ylabel)

    return fig, ax, data


def plot_rmse(x, y, xlabel, ylabel, **kwargs):
    """Plot RMSE scatter plot with reference line.
    
    This function creates a scatter plot of data points
    with a reference line and calculates RMSE.
    
    Parameters
    ----------
    x : array_like
        X values
    y : array_like
        Y values
    xlabel : str
        X-axis label
    ylabel : str
        Y-axis label
    **kwargs
        Additional keyword arguments passed to plotting functions
        
    Returns
    -------
    tuple
        Tuple of (fig, ax, rmse) containing matplotlib figure,
        axes object, and calculated RMSE
    """
    x = np.array(x)
    y = np.array(y)

    rmse = np.sqrt(metrics.mean_squared_error(x, y))

    fig, ax = plt.subplots(figsize=[4, 4])
    ax_rmse(ax, x, y)
    ax_setlabel(ax, xlabel, ylabel, **kwargs)

    return fig, ax, rmse


def plot_bin_stats(x, y, xlabel, ylabel, bins=None):
    """Plot binned statistics with confidence intervals.
    
    This function creates a plot showing mean, standard deviation,
    and min/max ranges for binned data.
    
    Parameters
    ----------
    x : array_like
        X values
    y : array_like
        Y values
    xlabel : str
        X-axis label
    ylabel : str
        Y-axis label
    bins : int, optional
        Number of bins, by default len(x)//10
        
    Returns
    -------
    tuple
        Tuple of (fig, ax) containing matplotlib figure and axes
    """
    x = np.array(x).flatten()
    y = np.array(y).flatten()
    if bins is None:
        bins = len(x) // 10

    fig, ax = plt.subplots(figsize=[6, 4], dpi=200)
    ax_bin_stats(ax, x, y, bins=bins)
    ax_setlabel(ax, xlabel, ylabel)
    ax.legend(loc="center left", bbox_to_anchor=(1.1, 0.5))
    return fig, ax


def ax_bin_stats(ax, x, y, bins):
    """Create binned statistics plot on given axes.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Matplotlib axes object
    x : array_like
        X values
    y : array_like
        Y values
    bins : int
        Number of bins
        
    Returns
    -------
    dict
        Dictionary containing binned statistics data
    """
    # ref line (y = 0)
    ax.axhline(y=0, color="gray")
    # mean
    bin_means, bin_edges, binnumber = stats.binned_statistic(x, y, bins=bins)
    bin_means[np.isnan(bin_means)] = 0.0
    bin_centers = bin_edges[1:] - (bin_edges[1] - bin_edges[0]) / 2
    ax.plot(bin_centers, bin_means, color="black", lw=1.5, label="mean")

    # max/min
    bin_maxs, bin_edges, binnumber = stats.binned_statistic(
        x, y, statistic="max", bins=bins
    )
    bin_mins, bin_edges, binnumber = stats.binned_statistic(
        x, y, statistic="min", bins=bins
    )
    bin_maxs[np.isnan(bin_maxs)] = 0.0
    bin_mins[np.isnan(bin_maxs)] = 0.0
    ax.fill_between(bin_centers, bin_mins, bin_maxs, color="silver", label="[min, max]")
    # std
    bin_stds, bin_edges, binnumber = stats.binned_statistic(
        x, y, statistic="std", bins=bins
    )
    bin_stds[np.isnan(bin_stds)] = 0.0
    ax.fill_between(
        bin_centers,
        bin_means - bin_stds,
        bin_means + bin_stds,
        color="gray",
        label="[mean-std, mean+std]",
    )
    data = {
        "grid": bin_centers,
        "mean": bin_means,
        "max": bin_maxs,
        "min": bin_mins,
        "std": bin_stds,
    }
    return data


def plot_colormap_lines(xs, ys, legends, xlabel, ylabel, colormap="GnBu"):
    """Plot multiple lines with colormap.
    
    This function creates a plot with multiple lines colored
    according to a colormap.
    
    Parameters
    ----------
    xs : list
        List of x data arrays
    ys : list
        List of y data arrays
    legends : list
        List of legend labels
    xlabel : str
        X-axis label
    ylabel : str
        Y-axis label
    colormap : str, optional
        Matplotlib colormap name, by default "GnBu"
        
    Returns
    -------
    tuple
        Tuple of (fig, ax) containing matplotlib figure and axes
    """
    fig, ax = plt.subplots(figsize=[6, 4], dpi=200)
    ax_colormap_lines(ax, xs, ys, legends, colormap)
    ax_setlabel(ax, xlabel, ylabel)
    ax.legend(loc="center left", bbox_to_anchor=(1.1, 0.5))
    return fig, ax


def ax_colormap_lines(
    ax, xs, ys, labels, scale=(0.0, 1.0), colormap="GnBu", fmt="%f", **kwargs
):
    """Plot multiple colored lines on axes.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Matplotlib axes object
    xs : list
        List of x data arrays
    ys : list
        List of y data arrays
    labels : list
        List of legend labels
    scale : tuple, optional
        Scale range for normalization, by default (0.0, 1.0)
    colormap : str, optional
        Matplotlib colormap name, by default "GnBu"
    fmt : str, optional
        Format string for legend labels, by default "%f"
    **kwargs
        Additional keyword arguments passed to ax.plot
        
    Returns
    -------
    None
        Function modifies the axes object in place
    """
    labels = np.array(labels)
    # normalization
    labels = (labels - labels.min()) / (labels.max() - labels.min())
    cm_scales = (labels - scale[0]) / (scale[1] - scale[0])
    for x, y, label, cm_scale in zip(xs, ys, labels, cm_scales, strict=False):
        ax.plot(
            x, y, color=plt.get_cmap(colormap)(cm_scale), label=fmt % label, **kwargs
        )
    ax.set_xlim(np.min(x), np.max(x))
