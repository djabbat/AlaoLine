import numpy as np
import pytest

from alaoline.process import AlaoProcess, RAW_FE2_FRACTION
from alaoline.design import fractional_factorial_2kminus1, add_center_points
from alaoline.rsm import ResponseSurface


def test_process_n2_preserves_fe2():
    p = AlaoProcess(atmosphere="n2")
    p.run(hours=36)
    assert p.fe2_fraction >= 0.60, f"N2 should preserve Fe2+ >=60%, got {p.fe2_fraction}"
    assert 0 < p.ip6_mg_100g <= p.ip6_0, f"IP6 should decrease or stay <= initial, got {p.ip6_mg_100g}"


def test_process_air_oxidizes_fe2():
    p = AlaoProcess(atmosphere="air")
    p.run(hours=36)
    assert p.fe2_fraction <= 0.10, f"air drying should drop Fe2+, got {p.fe2_fraction}"


def test_design_16_runs_resolution_v():
    df = fractional_factorial_2kminus1(k=5, seed=1)
    assert len(df) == 16
    assert list(df.columns) == ["Run", "A", "B", "C", "D", "E"]


def test_design_center_points():
    df = fractional_factorial_2kminus1(k=5, seed=1)
    df2 = add_center_points(df, n=3)
    assert len(df2) == 19


def test_rsm_fits_and_predicts():
    rng = np.random.default_rng(0)
    X = rng.uniform(-1, 1, size=(30, 3))
    y = 2.0 * X[:, 0] + 1.5 * X[:, 1] ** 2 - 0.8 * X[:, 0] * X[:, 2] + 1.0 + 0.01 * rng.normal(size=30)
    rsm = ResponseSurface().fit(X, y, factor_names=["A", "B", "C"])
    assert rsm.r2(X, y) > 0.95
    yhat = rsm.predict(X)
    assert yhat.shape == y.shape


def test_stationary_point_shape():
    rng = np.random.default_rng(1)
    X = rng.uniform(-1, 1, size=(40, 3))
    y = X.sum(axis=1) + 1.0
    rsm = ResponseSurface().fit(X, y, factor_names=["A", "B", "C"])
    sp = rsm.stationary_point()
    assert sp.shape == (3,)
