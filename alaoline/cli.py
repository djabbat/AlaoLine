"""Command-line interface for AlaoLine."""

from __future__ import annotations

import argparse
import json
import sys

import numpy as np
import pandas as pd

from .process import AlaoProcess
from .design import fractional_factorial_2kminus1, add_center_points
from .rsm import ResponseSurface


def cmd_simulate(args: argparse.Namespace) -> None:
    p = AlaoProcess(atmosphere=args.atmosphere, temperature=args.temp)
    p.run(hours=args.hours)
    print(f"atmosphere : {p.atmosphere}")
    print(f"final Fe2+ fraction : {p.fe2_fraction:.3f}")
    print(f"final IP6 (mg/100g): {p.ip6_mg_100g:.1f}")
    if args.plot is not None:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        plt.plot(p.t, p.fe2, label="Fe2+ fraction")
        plt.plot(p.t, p.ip6, label="IP6 mg/100g")
        plt.xlabel("time (h)"); plt.ylabel("value")
        plt.legend(); plt.tight_layout()
        plt.savefig(args.plot)


def cmd_design(args: argparse.Namespace) -> None:
    df = fractional_factorial_2kminus1(k=5, seed=args.seed, coded=args.coded)
    df = add_center_points(df, n=args.center, coded=args.coded)
    if args.out:
        df.to_csv(args.out, index=False)
        print(f"wrote {len(df)} runs to {args.out}")
    else:
        print(df.to_string(index=False))


def cmd_fit(args: argparse.Namespace) -> None:
    data = pd.read_csv(args.data)
    ycol = args.target
    xcols = [c for c in data.columns if c not in ("Run", ycol)]
    X = data[xcols].astype(float).to_numpy()
    y = data[ycol].astype(float).to_numpy()
    rsm = ResponseSurface().fit(X, y, factor_names=xcols)
    print("R2 =", round(rsm.r2(X, y), 4))
    print("stationary point (coded):", np.round(rsm.stationary_point(), 3).tolist())
    if args.out:
        with open(args.out, "w") as fh:
            json.dump({"r2": rsm.r2(X, y), "coef": rsm.coef_.tolist()}, fh, indent=2)
        print(f"wrote model to {args.out}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="alaoline", description="AlaoLine: N2 short-sprouting model")
    sub = parser.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("simulate", help="simulate the process")
    s.add_argument("--hours", type=float, default=36.0)
    s.add_argument("--atmosphere", choices=["n2", "air"], default="n2")
    s.add_argument("--temp", type=float, default=22.0)
    s.add_argument("--plot", default=None, help="save plot to path")
    s.set_defaults(func=cmd_simulate)

    d = sub.add_parser("design", help="generate 2^(5-1) design")
    d.add_argument("--runs", type=int, default=16)
    d.add_argument("--center", type=int, default=3)
    d.add_argument("--seed", type=int, default=None)
    d.add_argument("--coded", action="store_true")
    d.add_argument("--out", default=None)
    d.set_defaults(func=cmd_design)

    f = sub.add_parser("fit", help="fit RSM from CSV")
    f.add_argument("--data", required=True)
    f.add_argument("--target", required=True)
    f.add_argument("--out", default=None)
    f.set_defaults(func=cmd_fit)

    return parser


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
