"""Command-line interface: ``tls-finder``."""

from __future__ import annotations

import argparse
import sys

import pandas as pd

from . import __version__
from .core import (
    DEFAULT_B_LABELS,
    DEFAULT_FDC_LABELS,
    DEFAULT_HEV_LABELS,
    DEFAULT_T_LABELS,
    classify_cells,
    find_lymphoid_aggregates,
)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="tls-finder",
        description="Detect tertiary lymphoid structures from phenotyped cell coordinates.",
    )
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    p.add_argument("input", help="CSV file with one row per cell.")
    p.add_argument(
        "--mode",
        choices=["aggregates", "cells"],
        default="aggregates",
        help="'aggregates': discrete lymphoid aggregates (default). "
        "'cells': per-cell classification including HEV/fDC-associated TLS.",
    )
    p.add_argument("-r", "--radius", type=float, default=None,
                   help="Neighbourhood radius in coordinate units "
                        "(default: 160 for aggregates, 80 for cells).")
    p.add_argument("--min-cells", type=int, default=50,
                   help="Minimum B + T cells within the radius (default: 50).")
    p.add_argument("--b-labels", nargs="+", default=list(DEFAULT_B_LABELS),
                   help="Phenotype values counted as B cells (default: CD20).")
    p.add_argument("--t-labels", nargs="+", default=list(DEFAULT_T_LABELS),
                   help="Phenotype values counted as T cells (default: CD3 CD8).")
    p.add_argument("--hev-labels", nargs="+", default=list(DEFAULT_HEV_LABELS),
                   help="Phenotype values counted as HEV, cells mode only (default: HEV).")
    p.add_argument("--fdc-labels", nargs="+", default=list(DEFAULT_FDC_LABELS),
                   help="Phenotype values counted as follicular dendritic cells, "
                        "cells mode only (default: fDC).")
    p.add_argument("--x-col", default="X")
    p.add_argument("--y-col", default="Y")
    p.add_argument("--phenotype-col", default="Phenotype")
    p.add_argument("-o", "--output", help="Write results to this CSV file (default: print to stdout).")
    p.add_argument("--plot", help="Save a figure of the results to this image file.")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    df = pd.read_csv(args.input)
    cols = dict(x_col=args.x_col, y_col=args.y_col, phenotype_col=args.phenotype_col)

    try:
        if args.mode == "aggregates":
            radius = args.radius if args.radius is not None else 160.0
            result = find_lymphoid_aggregates(
                df, radius=radius, min_cells=args.min_cells,
                b_labels=args.b_labels, t_labels=args.t_labels, **cols,
            )
            summary = f"Found {len(result)} lymphoid aggregate(s) at radius {radius:g}."
        else:
            radius = args.radius if args.radius is not None else 80.0
            classes = classify_cells(
                df, radius=radius, min_cells=args.min_cells, b_labels=args.b_labels,
                t_labels=args.t_labels, hev_labels=args.hev_labels,
                fdc_labels=args.fdc_labels, **cols,
            )
            result = df.assign(TLS_Classification=classes)
            counts = classes.value_counts().to_dict()
            summary = f"Cell classification at radius {radius:g}: {counts}"
    except (ValueError, TypeError) as exc:
        print(f"tls-finder: error: {exc}", file=sys.stderr)
        return 2

    if args.output:
        result.to_csv(args.output, index=False)
    else:
        result.to_csv(sys.stdout, index=False)
    print(summary, file=sys.stderr)

    if args.plot:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        from .plotting import plot_aggregates, plot_cell_classes

        if args.mode == "aggregates":
            ax = plot_aggregates(df, result, x_col=args.x_col, y_col=args.y_col)
        else:
            ax = plot_cell_classes(df, classes, x_col=args.x_col, y_col=args.y_col)
        ax.figure.savefig(args.plot, dpi=150, bbox_inches="tight")
        plt.close(ax.figure)

    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
