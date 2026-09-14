"""Response-surface modelling (RSM) for the AlaoLine process.

Fits a quadratic response surface to experimental data (coded factors) and
locates the stationary point, as used in WP3/T3.4 to identify optimal process
settings (T, t, citrate, N2 flow).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from itertools import combinations

# AlaoLine uses a resolution-V fractional factorial (2^(5-1)); for the RSM
# we treat the numeric factors as coded in [-1, +1].


def design_matrix(X: np.ndarray, intercept: bool = True) -> np.ndarray:
    """Expand X (n x p, coded) to a full quadratic basis.

    Columns: 1, x_i, x_i^2, x_i x_j (i<j).
    """
    n, p = X.shape
    cols = []
    if intercept:
        cols.append(np.ones(n))
    cols.append(X)
    cols.append(X**2)
    for i, j in combinations(range(p), 2):
        cols.append(X[:, i] * X[:, j])
    return np.column_stack(cols)


class ResponseSurface:
    """OLS quadratic response surface on coded factors."""

    def __init__(self) -> None:
        self.coef_: np.ndarray | None = None
        self.factor_names_: list[str] | None = None

    def fit(self, X: np.ndarray, y: np.ndarray, factor_names: list[str] | None = None) -> "ResponseSurface":
        Phi = design_matrix(X)
        self.coef_ = np.linalg.lstsq(Phi, y, rcond=None)[0]
        self.factor_names_ = factor_names or [f"x{i}" for i in range(X.shape[1])]
        return self

    def predict(self, X: np.ndarray) -> np.ndarray:
        if self.coef_ is None:
            raise RuntimeError("fit() must be called first")
        return design_matrix(X) @ self.coef_

    def r2(self, X: np.ndarray, y: np.ndarray) -> float:
        yhat = self.predict(X)
        ss_res = float(np.sum((y - yhat) ** 2))
        ss_tot = float(np.sum((y - np.mean(y)) ** 2))
        return 1.0 - ss_res / ss_tot

    def stationary_point(self) -> np.ndarray:
        """Solve for the stationary point of the fitted quadratic (coded space)."""
        if self.coef_ is None:
            raise RuntimeError("fit() must be called first")
        p = len(self.factor_names_)
        n_coef = 1 + p + p + p * (p - 1) // 2
        c = self.coef_
        # linear terms (indices 1..p)
        b = c[1 : 1 + p].astype(float)
        # quadratic terms (indices 1+p .. 1+2p-1)
        q = c[1 + p : 1 + 2 * p].astype(float)
        # cross terms
        cr = c[1 + 2 * p :].astype(float)
        B = np.diag(2.0 * q)
        idx = 0
        for i, j in combinations(range(p), 2):
            B[i, j] = cr[idx]
            B[j, i] = cr[idx]
            idx += 1
        return -np.linalg.solve(B + B.T, b + b.T) / 2.0

    def summary(self) -> str:
        names = self.factor_names_ or []
        return "ResponseSurface(factors={}, coefs={})".format(names, self.coef_)
