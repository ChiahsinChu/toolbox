# SPDX-License-Identifier: LGPL-3.0-or-later
import numpy as np
import statsmodels.api as sm
from scipy import stats


class Optimizer:
    def __init__(self, method="ols_cut") -> None:
        self.method = method

    def run(self, x, y, **kwargs):
        assert len(x) == len(y)
        output = getattr(self, "_run_%s" % self.method)(x, y, **kwargs)
        return output

    def _run_ols(self, x, y):
        result = stats.linregress(x=x, y=y)
        output = -result.intercept / result.slope
        return output

    def _run_ols_cut(self, x, y, nstep=4):
        l_cut = min(len(x), nstep)
        result = stats.linregress(x=x[-l_cut:], y=y[-l_cut:])
        output = -result.intercept / result.slope
        return output

    def _run_wls(self, X, y):
        # fit linear regression model
        X = sm.add_constant(x)
        wt = np.exp(-(np.array(y) ** 2) / 0.1)
        fit_wls = sm.WLS(y, X, weights=wt).fit()
        return -fit_wls.params[0] / fit_wls.params[1]

    def _run_wls_cut(self, X, y, nstep=4):
        l_cut = min(len(X), nstep)
        X = X[-l_cut:]
        y = y[-l_cut:]
        X = sm.add_constant(x)
        wt = np.exp(-(np.array(y) ** 2) / 0.1)
        fit_wls = sm.WLS(y, X, weights=wt).fit()
        return -fit_wls.params[0] / fit_wls.params[1]
