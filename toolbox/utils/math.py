# SPDX-License-Identifier: LGPL-3.0-or-later
import math

import numpy as np
from scipy import stats
from scipy.special import erf
from sklearn import metrics

from .utils import get_bins_from_bin_edge


def gaussian_func(x, mu=0, sigma=1):
    """Gaussian probability density function.
    
    Parameters
    ----------
    x : array_like
        Input values
    mu : float, optional
        Mean of the distribution, by default 0
    sigma : float, optional
        Standard deviation, by default 1
        
    Returns
    -------
    array_like
        Gaussian probability density values
    """
    return stats.norm.pdf(x, loc=mu, scale=sigma)


def gaussian_int(x, mu, sigma):
    """Gaussian cumulative distribution function.
    
    Parameters
    ----------
    x : array_like
        Input values
    mu : float
        Mean of the distribution
    sigma : float
        Standard deviation
        
    Returns
    -------
    array_like
        Gaussian cumulative distribution values
    """
    coeff = 1 / 2
    return coeff * erf((x - mu) / (np.sqrt(2) * sigma))


def block_ave(_x, _y, l_block):
    """Calculate block averages of data.
    
    Parameters
    ----------
    _x : array_like
        Input x values
    _y : array_like
        Input y values
    l_block : int
        Block size for averaging
        
    Returns
    -------
    tuple
        Tuple of (x_ave, y_ave) containing block-averaged values
    """
    assert len(_x) == len(_y)
    n_block = math.floor(len(_x) / l_block)
    x = _x[: (n_block * l_block)]
    y = _y[: (n_block * l_block)]
    x = np.reshape(x, (-1, l_block)).mean(axis=-1)
    y = np.reshape(y, (-1, l_block)).mean(axis=-1)
    return x, y


def get_dev(x, y):
    """Calculate derivative dy/dx.
    
    Parameters
    ----------
    x : array_like
        Input x values
    y : array_like
        Input y values
        
    Returns
    -------
    tuple
        Tuple of (x_mid, dy/dx) where x_mid are midpoints
    """
    x = np.array(x)
    y = np.array(y)
    delta_x = np.diff(x)
    delta_y = np.diff(y)
    out_x = (x[1:] + x[:-1]) / 2
    return out_x, delta_y / delta_x


def get_bin_ave(data, bin_size=10):
    """Calculate binned averages of data.
    
    Parameters
    ----------
    data : array_like
        Input data to bin
    bin_size : int, optional
        Size of each bin, by default 10
        
    Returns
    -------
    array_like
        Binned averages
    """
    total_size = len(data)
    out_size = bin_size * total_size // bin_size
    data = np.reshape(data[:out_size], (-1, bin_size))
    data = np.mean(data, axis=-1)
    return data


def get_int(x, y):
    """Calculate numerical integral using trapezoidal rule.
    
    Parameters
    ----------
    x : array_like
        Input x values
    y : array_like
        Input y values
        
    Returns
    -------
    tuple
        Tuple of (x_mid, integral) where x_mid are midpoints
    """
    x = np.array(x)
    y = np.array(y)
    ave_y = (y[1:] + y[:-1]) / 2
    delta_x = np.diff(x)
    x = (x[1:] + x[:-1]) / 2
    return x, np.cumsum(ave_y * delta_x)


def get_int_array(x, ys):
    """Calculate numerical integrals for multiple y arrays.
    
    Parameters
    ----------
    x : array_like
        Input x values
    ys : array_like
        Multiple y arrays with shape (n_arrays, n_points)
        
    Returns
    -------
    tuple
        Tuple of (x_mid, integrals) where x_mid are midpoints
    """
    x = np.array(x)
    ys = np.array(ys)
    ave_y = (ys[:, 1:] + ys[:, :-1]) / 2
    delta_x = np.diff(x)
    x = (x[1:] + x[:-1]) / 2
    return x, np.cumsum(ave_y * delta_x.reshape(1, -1), axis=-1)


def cumave(data):
    """Calculate cumulative average of data.
    
    Parameters
    ----------
    data : array_like
        Input data
        
    Returns
    -------
    array_like
        Cumulative average values
    """
    cum_sum = data.cumsum()
    cum_ave = cum_sum / (np.arange(len(data)) + 1)
    return cum_ave


def interp(x, dataset):
    """Interpolate y values from dataset.
    
    Parameters
    ----------
    x : array_like
        X values to interpolate at
    dataset : array_like
        Dataset with dataset[0] = x_values, dataset[1] = y_values
        
    Returns
    -------
    array_like
        Interpolated y values
    """
    y = np.interp(x, xp=dataset[0], fp=dataset[1])
    return y


def handle_zero_division(x, y, threshold=None):
    """Handle division by zero with optional threshold.
    
    This function safely divides x by y, handling zero
    division with optional threshold masking.
    
    Parameters
    ----------
    x : array_like
        Numerator values
    y : array_like
        Denominator values
    threshold : float, optional
        Threshold below which y values are set to zero
        
    Returns
    -------
    array_like
        Result of x/y with safe division
    """
    if threshold is not None:
        mask = np.abs(y) <= threshold
        y[np.nonzero(mask)] = 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        result = np.true_divide(x, y)
        # replace NaN and Inf values with 0
        result[~np.isfinite(result)] = 0
    return result


def error_test(y_true, y_pred):
    """Calculate regression error metrics.
    
    This function calculates various regression error metrics
    between true and predicted values.
    
    Parameters
    ----------
    y_true : array_like
        True values
    y_pred : array_like
        Predicted values
        
    Returns
    -------
    dict
        Dictionary containing error metrics:
        - max_err: Maximum error
        - mae: Mean absolute error
        - rmse: Root mean squared error
        - r2: R-squared score
        - mape: Mean absolute percentage error
        
    References
    ----------
    https://scikit-learn.org/stable/modules/model_evaluation.html#regression-metrics
    """
    results_dict = {}
    y_true = np.reshape(y_true, (-1,))
    y_pred = np.reshape(y_pred, (-1,))

    results_dict["max_err"] = metrics.max_error(y_true, y_pred)
    results_dict["mae"] = metrics.mean_absolute_error(y_true, y_pred)
    results_dict["rmse"] = np.sqrt(metrics.mean_squared_error(y_true, y_pred))
    results_dict["r2"] = metrics.r2_score(y_true, y_pred)
    results_dict["mape"] = metrics.mean_absolute_percentage_error(y_true, y_pred)
    return results_dict


def vec_project(vec, unit_vec):
    """Project vector onto unit vector.
    
    Parameters
    ----------
    vec : array_like
        Vector to project
    unit_vec : array_like
        Unit vector to project onto
        
    Returns
    -------
    array_like
        Projected vector
    """
    return np.dot(vec, unit_vec) * unit_vec


def gaussian_filter(data, bin_edge, sigma: float, weight=None):
    """Apply Gaussian filter to data.
    
    Parameters
    ----------
    data : array_like
        Input data to filter
    bin_edge : array_like
        Bin edges for filtering
    sigma : float
        Gaussian standard deviation
    weight : array_like, optional
        Weights for each data point
        
    Returns
    -------
    tuple
        Tuple of (bins, filtered_data)
    """
    data = np.reshape(data, (-1, 1))
    bins = get_bins_from_bin_edge(bin_edge)
    bins = np.reshape(bins, (1, -1))
    weight = np.ones_like(data) if weight is None else np.reshape(weight, (-1, 1))
    output = (
        np.exp(-(((bins - data) / sigma) ** 2)) / (np.sqrt(2 * np.pi) * sigma) * weight
    )
    return bins, output.sum(axis=0)
