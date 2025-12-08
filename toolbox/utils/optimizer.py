# SPDX-License-Identifier: LGPL-3.0-or-later
"""Optimizer module.

This module provides classes and functions for optimization
using various linear regression methods.
"""

import numpy as np
import statsmodels.api as sm
from scipy import stats


class Optimizer:
    """Optimizer for linear regression models.
    
    This class provides various optimization methods for
    fitting linear models to data.
    """
    
    def __init__(self, method="ols_cut") -> None:
        """Initialize Optimizer.
        
        Parameters
        ----------
        method : str, optional
            Optimization method to use, by default "ols_cut"
        """
        self.method = method

    def run(self, x, y, **kwargs):
        """Run optimization with specified method.
        
        Parameters
        ----------
        x : array_like
            Input features
        y : array_like
            Target values
        **kwargs
            Additional keyword arguments
            
        Returns
        -------
        array_like
            Optimization result
        """
        assert len(x) == len(y)
        output = getattr(self, f"_run_{self.method}")(x, y, **kwargs)
        return output

    def _run_ols(self, x, y):
        """Run ordinary least squares optimization.
        
        Parameters
        ----------
        x : array_like
            Input features
        y : array_like
            Target values
            
        Returns
        -------
        array_like
            Optimization result [-intercept/slope]
        """
        result = stats.linregress(x=x, y=y)
        output = -result.intercept / result.slope
        return output

    def _run_ols_cut(self, x, y, nstep=4):
        """Run OLS with cutoff on recent data.
        
        Parameters
        ----------
        x : array_like
            Input features
        y : array_like
            Target values
        nstep : int, optional
            Number of recent steps to use, by default 4
            
        Returns
        -------
        array_like
            Optimization result [-intercept/slope]
        """
        l_cut = min(len(x), nstep)
        result = stats.linregress(x=x[-l_cut:], y=y[-l_cut:])
        output = -result.intercept / result.slope
        return output

    def _run_wls(self, X, y):
        """Run weighted least squares optimization.
        
        Parameters
        ----------
        X : array_like
            Input features with constant term
        y : array_like
            Target values
            
        Returns
        -------
        array_like
            Optimization result [-intercept/slope]
        """
        # fit linear regression model
        X = sm.add_constant(self.x)
        wt = np.exp(-(np.array(self.y) ** 2) / 0.1)
        fit_wls = sm.WLS(self.y, X, weights=wt).fit()
        return -fit_wls.params[0] / fit_wls.params[1]

    def _run_wls_cut(self, X, y, nstep=4):
        """Run WLS with cutoff on recent data.
        
        Parameters
        ----------
        X : array_like
            Input features with constant term
        y : array_like
            Target values
        nstep : int, optional
            Number of recent steps to use, by default 4
            
        Returns
        -------
        array_like
            Optimization result [-intercept/slope]
        """
        l_cut = min(len(X), nstep)
        X = X[-l_cut:]
        y = y[-l_cut:]
        X = sm.add_constant(X)
        wt = np.exp(-(np.array(y) ** 2) / 0.1)
        fit_wls = sm.WLS(y, X, weights=wt).fit()
        return -fit_wls.params[0] / fit_wls.params[1]
