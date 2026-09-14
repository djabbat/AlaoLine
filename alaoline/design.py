"""Fractional-factorial experimental design (2^(5-1), resolution V).

Used in WP3 (process optimisation) to map the factor space:

    A: drying atmosphere   (Air / N2)
    B: citrate soaking      (Water / Citrate 50 mM pH 4.5)
    C: sprouting duration   (6 h / 12 h)
    D: N2 flow rate         (0.5 / 2.0 L/min)
    E: drying temperature   (30 / 42 C)

The 2^(5-1) design is the full 2^4 factorial in A, B, C, D (16 runs) with the
generating relation E = A*B*C*D. This is a resolution V design: it resolves
all main effects and two-way interactions. Centre points are added for
curvature/lack-of-fit testing.
"""

from __future__ import annotations

import itertools

import numpy as np
import pandas as pd

FACTORS = ["A", "B", "C", "D", "E"]
LEVELS = {
    "A": ("Air", "N2"),
    "B": ("Water", "Citrate"),
    "C": (6.0, 12.0),
    "D": (0.5, 2.0),
    "E": (30.0, 42.0),
}


def fractional_factorial_2kminus1(
    k: int = 5, seed: int | None = None, coded: bool = False
) -> pd.DataFrame:
    """Return the 2^(5-1) fractional-factorial design (16 runs).

    Full 2^4 factorial in A, B, C, D with E = A*B*C*D.
    """
    if k != 5:
        raise NotImplementedError("only the 2^(5-1) design is implemented")
    # full 2^4 factorial in A..D (16 rows), coded +/-1
    signs = np.array(list(itertools.product([1, -1], repeat=4)))  # (16,4)
    e = signs[:, 0] * signs[:, 1] * signs[:, 2] * signs[:, 3]      # E = ABCD
    design = np.column_stack([signs, e])                            # (16,5)

    rng = np.random.default_rng(seed)
    runs = rng.permutation(16)

    rows = []
    for r in runs:
        row = {}
        for i, f in enumerate(FACTORS):
            level = 0 if design[r, i] == 1 else 1
            row[f] = level if coded else LEVELS[f][level]
        rows.append(row)
    df = pd.DataFrame(rows)
    df.insert(0, "Run", np.arange(1, 17))
    return df


def add_center_points(df: pd.DataFrame, n: int = 3, coded: bool = False) -> pd.DataFrame:
    """Append n centre points (all factors at their midpoint)."""
    mid = {}
    for f in FACTORS:
        lo, hi = LEVELS[f]
        if isinstance(lo, (int, float)) and isinstance(hi, (int, float)):
            mid[f] = 0.0 if coded else (lo + hi) / 2.0
        else:
            mid[f] = 0.0 if coded else "N2" if f == "A" else "Citrate"
    center = pd.DataFrame([mid for _ in range(n)])
    center.insert(0, "Run", np.arange(len(df) + 1, len(df) + n + 1))
    return pd.concat([df, center], ignore_index=True)
