import math

import numpy as np
from scipy.special import erf
from sklearn import metrics


def gaussian_func(x, mu, sigma):
    """
    Gaussian function
    """
    coeff = 1 / (sigma * np.sqrt(2 * np.pi))
    return coeff * np.exp(-np.power(x - mu, 2.) / (2 * np.power(sigma, 2.)))


def gaussian_int(x, mu, sigma):
    """
    Gaussian integral
    """
    coeff = 1 / 2
    return coeff * erf((x - mu) / (np.sqrt(2) * sigma))


def block_ave(_x, _y, l_block):
    assert len(_x) == len(_y)
    n_block = math.floor(len(_x) / l_block)
    x = _x[:(n_block * l_block)]
    y = _y[:(n_block * l_block)]
    x = np.reshape(x, (-1, l_block)).mean(axis=-1)
    y = np.reshape(y, (-1, l_block)).mean(axis=-1)
    return x, y


def get_dev(x, y):
    x = np.array(x)
    y = np.array(y)
    delta_x = np.diff(x)
    delta_y = np.diff(y)
    out_x = (x[1:] + x[:-1]) / 2
    return out_x, delta_y / delta_x


def get_bin_ave(data, bin_size=10):
    total_size = len(data)
    out_size = bin_size * total_size // bin_size
    data = np.reshape(data[:out_size], (-1, bin_size))
    data = np.mean(data, axis=-1)
    return data


def get_int(x, y):
    x = np.array(x)
    y = np.array(y)
    ave_y = (y[1:] + y[:-1]) / 2
    delta_x = np.diff(x)
    x = (x[1:] + x[:-1]) / 2
    return x, np.cumsum(ave_y * delta_x)


def get_int_array(x, ys):
    x = np.array(x)
    ys = np.array(ys)
    ave_y = (ys[:, 1:] + ys[:, :-1]) / 2
    delta_x = np.diff(x)
    x = (x[1:] + x[:-1]) / 2
    return x, np.cumsum(ave_y * delta_x.reshape(1, -1), axis=-1)


def cumave(data):
    cum_sum = data.cumsum()
    cum_ave = cum_sum / (np.arange(len(data)) + 1)
    return cum_ave


def interp(x, dataset):
    y = np.interp(x, xp=dataset[0], fp=dataset[1])
    return y


def handle_zero_division(x, y, threshold=None):
    """
    define a function to handle the runtime warning
    """
    if threshold is not None:
        mask = (np.abs(y) <= threshold)
        y[np.nonzero(mask)] = 0.
    with np.errstate(divide='ignore', invalid='ignore'):
        result = np.true_divide(x, y)
        # replace NaN and Inf values with 0
        result[~np.isfinite(result)] = 0
    return result


def error_test(y_true, y_pred):
    """
    https://scikit-learn.org/stable/modules/model_evaluation.html#regression-metrics
    """
    results_dict = {}
    y_true = np.reshape(y_true, (-1, ))
    y_pred = np.reshape(y_pred, (-1, ))

    results_dict["max_err"] = metrics.max_error(y_true, y_pred)
    results_dict["mae"] = metrics.mean_absolute_error(y_true, y_pred)
    results_dict["rmse"] = np.sqrt(metrics.mean_squared_error(y_true, y_pred))
    results_dict["r2"] = metrics.r2_score(y_true, y_pred)
    results_dict["mape"] = metrics.mean_absolute_percentage_error(
        y_true, y_pred)
    return results_dict


def vec_project(vec, unit_vec):
    return np.dot(vec, unit_vec) * unit_vec
