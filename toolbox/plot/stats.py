# SPDX-License-Identifier: LGPL-3.0-or-later
"""Statistical plotting module.

This module provides classes for creating statistical plots
and performing statistical tests, including finite difference
method testing.
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

# from .core import ax_setlabel


class FDMTest:
    """Finite difference method test for derivative.
    
    This class implements finite difference methods
    to test numerical derivatives against analytical solutions.
    """

    def __init__(self, **kwargs) -> None:
        """Initialize FDMTest.
        
        Parameters
        ----------
        **kwargs
            Additional keyword arguments passed to plt.subplots
        """
        self.fig, self.axs = plt.subplots(
            nrows=2, figsize=[4, 6], sharex="all", **kwargs
        )

    def setup(self, x, y, dydx):
        """Set up finite difference test.
        
        Parameters
        ----------
        x : array_like
            X values
        y : array_like
            Y values
        dydx : array_like
            Analytical derivative values
        """
        self.stats_dict = {}
        ax = self.axs[0]
        ax.scatter(x, dydx, color="blue", label="original data")
        fit_output = stats.linregress(x, dydx)
        ax.plot(
            x,
            x * fit_output.slope + fit_output.intercept,
            "--",
            color="red",
            label="fitted line",
        )
        self.stats_dict["slope"] = fit_output.slope
        self.stats_dict["intercept"] = fit_output.intercept
        self.stats_dict["rvalue"] = fit_output.rvalue

        ax = self.axs[1]
        ref_data_diff = np.diff(y, axis=0).reshape(-1)
        test_data_diff = 0.5 * np.diff(x, axis=0) * (dydx[1:] + dydx[:-1])
        ax.scatter(x[1:], ref_data_diff, color="blue", label="original data")
        ax.plot(x[1:], test_data_diff, "--", color="red", label="fitted line")
        self.stats_dict["mean error"] = np.mean(ref_data_diff - test_data_diff)
