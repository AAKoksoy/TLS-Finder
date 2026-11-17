"""Plotting helpers for TLS Finder results (requires matplotlib)."""

from __future__ import annotations

from typing import Mapping

import pandas as pd


def _axes(ax):
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(figsize=(8, 8))
    return ax


def _style(ax, title, invert_y):
    if invert_y:
        ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.set_xlabel("X (microns)")
    ax.set_ylabel("Y (microns)")
    ax.set_title(title)


def plot_phenotypes(
    df: pd.DataFrame,
    colors: Mapping[str, str] | None = None,
    ax=None,
    s: float = 1.0,
    invert_y: bool = True,
    x_col: str = "X",
    y_col: str = "Y",
    phenotype_col: str = "Phenotype",
):
    """Scatter all cells coloured by phenotype. Returns the matplotlib Axes."""
    ax = _axes(ax)
    phenos = list(colors) if colors else sorted(df[phenotype_col].astype(str).unique())
    for p in phenos:
        sub = df[df[phenotype_col] == p]
        kw = {"color": colors[p]} if colors else {}
        ax.scatter(sub[x_col], sub[y_col], s=s, label=p, **kw)
    _style(ax, "Cell distribution by phenotype", invert_y)
    ax.legend(markerscale=max(1, 10 / max(s, 0.1)), loc="center left", bbox_to_anchor=(1, 0.5))
    return ax


def plot_aggregates(
    df: pd.DataFrame,
    aggregates: pd.DataFrame,
    ax=None,
    s: float = 1.0,
    invert_y: bool = True,
    x_col: str = "X",
    y_col: str = "Y",
):
    """Plot all cells in grey with lymphoid aggregate centres highlighted.

    ``aggregates`` is the output of
    :func:`tls_finder.find_lymphoid_aggregates`.
    """
    ax = _axes(ax)
    ax.scatter(df[x_col], df[y_col], color="lightgray", s=s, label="Cells")
    ax.scatter(
        aggregates["center_x"],
        aggregates["center_y"],
        color="red",
        edgecolor="black",
        s=100,
        alpha=0.8,
        label="Lymphoid aggregate centres",
    )
    _style(ax, f"Lymphoid aggregates (n = {len(aggregates)})", invert_y)
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    return ax


_CLASS_COLORS = {
    "Not TLS": "lightgray",
    "Lymphoid Aggregate": "red",
    "Primary Follicle-like TLS": "cyan",
}


def plot_cell_classes(
    df: pd.DataFrame,
    classes: pd.Series,
    ax=None,
    s: float = 1.0,
    invert_y: bool = True,
    x_col: str = "X",
    y_col: str = "Y",
):
    """Plot cells coloured by :func:`tls_finder.classify_cells` output."""
    ax = _axes(ax)
    for label, color in _CLASS_COLORS.items():
        mask = (classes == label).to_numpy()
        if mask.any():
            ax.scatter(df.loc[mask, x_col], df.loc[mask, y_col], color=color, s=s, label=label)
    _style(ax, "Cells by TLS classification", invert_y)
    ax.legend(markerscale=max(1, 10 / max(s, 0.1)), loc="center left", bbox_to_anchor=(1, 0.5))
    return ax
