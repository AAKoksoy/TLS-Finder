"""Core TLS Finder algorithms.

Two complementary analyses are provided:

* :func:`find_lymphoid_aggregates` - a greedy, density-seeded search that
  reports one non-overlapping lymphoid aggregate per qualifying seed cell
  (the algorithm described in Koksoy et al., bioRxiv 2024).
* :func:`classify_cells` - a per-cell classification that labels each cell by
  the B/T-cell density and HEV presence in its neighbourhood.

All coordinates are assumed to be in the same unit as ``radius``
(typically microns).
"""

from __future__ import annotations

from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from scipy.spatial import KDTree

#: Default phenotype labels counted as B cells.
DEFAULT_B_LABELS: tuple[str, ...] = ("CD20",)
#: Default phenotype labels counted as T cells.
DEFAULT_T_LABELS: tuple[str, ...] = ("CD3", "CD8")
#: Default phenotype labels counted as high endothelial venules.
DEFAULT_HEV_LABELS: tuple[str, ...] = ("HEV",)
#: Default phenotype labels counted as follicular dendritic cells.
DEFAULT_FDC_LABELS: tuple[str, ...] = ("fDC",)

LYMPHOID_AGGREGATE = "Lymphoid Aggregate"
PRIMARY_FOLLICLE_LIKE_TLS = "Primary Follicle-like TLS"
NOT_TLS = "Not TLS"


def _as_tuple(labels: str | Iterable[str]) -> tuple[str, ...]:
    if isinstance(labels, str):
        return (labels,)
    labels = tuple(labels)
    if not labels:
        raise ValueError("Phenotype label lists must not be empty.")
    return labels


def validate_input(
    df: pd.DataFrame,
    x_col: str = "X",
    y_col: str = "Y",
    phenotype_col: str = "Phenotype",
) -> None:
    """Check that ``df`` has the columns and values TLS Finder needs.

    Raises
    ------
    TypeError
        If ``df`` is not a :class:`pandas.DataFrame`.
    ValueError
        If required columns are missing or coordinates are non-numeric or
        missing.
    """
    if not isinstance(df, pd.DataFrame):
        raise TypeError("Input must be a pandas DataFrame.")
    missing = [c for c in (x_col, y_col, phenotype_col) if c not in df.columns]
    if missing:
        raise ValueError(
            f"Input is missing required column(s): {missing}. "
            f"Available columns: {list(df.columns)}"
        )
    for col in (x_col, y_col):
        if not pd.api.types.is_numeric_dtype(df[col]):
            raise ValueError(f"Coordinate column '{col}' must be numeric.")
        if df[col].isna().any():
            raise ValueError(f"Coordinate column '{col}' contains missing values.")


def _check_params(radius: float, min_cells: int) -> None:
    if not np.isfinite(radius) or radius <= 0:
        raise ValueError("radius must be a positive number.")
    if int(min_cells) != min_cells or min_cells < 1:
        raise ValueError("min_cells must be a positive integer.")


def find_lymphoid_aggregates(
    df: pd.DataFrame,
    radius: float = 160.0,
    min_cells: int = 50,
    b_labels: str | Sequence[str] = DEFAULT_B_LABELS,
    t_labels: str | Sequence[str] = DEFAULT_T_LABELS,
    x_col: str = "X",
    y_col: str = "Y",
    phenotype_col: str = "Phenotype",
    return_labels: bool = False,
) -> pd.DataFrame | tuple[pd.DataFrame, pd.Series]:
    """Detect non-overlapping lymphoid aggregates.

    Cells are visited in input order. Each B or T cell not yet assigned to an
    aggregate is a candidate centre. If its circular neighbourhood of
    ``radius`` contains at least ``min_cells`` unassigned B or T cells (the
    candidate itself included), it is reported as a lymphoid aggregate centre
    and every cell in that neighbourhood, of any phenotype, is marked as
    assigned so it cannot seed or be counted towards a later aggregate.

    Parameters
    ----------
    df
        One row per cell, with coordinate and phenotype columns.
    radius
        Neighbourhood radius, in the same unit as the coordinates.
    min_cells
        Minimum number of unassigned B + T cells within ``radius``.
    b_labels, t_labels
        Phenotype values counted as B and T cells.
    x_col, y_col, phenotype_col
        Column names.
    return_labels
        If ``True``, also return a per-cell Series giving the aggregate each
        cell was assigned to (``-1`` for unassigned cells).

    Returns
    -------
    pandas.DataFrame
        One row per aggregate with columns ``aggregate_id``, ``center_index``
        (the index label of the seed cell in ``df``), ``center_x``,
        ``center_y``, ``classification``, ``n_cells`` (cells assigned to the
        aggregate), ``n_b_cells`` and ``n_t_cells``.
    pandas.Series, optional
        Per-cell aggregate assignment, indexed like ``df``.
    """
    validate_input(df, x_col, y_col, phenotype_col)
    _check_params(radius, min_cells)
    b_labels = _as_tuple(b_labels)
    t_labels = _as_tuple(t_labels)

    coords = df[[x_col, y_col]].to_numpy(dtype=float)
    phenotypes = df[phenotype_col].to_numpy()
    is_b = np.isin(phenotypes, b_labels)
    is_t = np.isin(phenotypes, t_labels)
    is_bt = is_b | is_t

    labels = np.full(len(df), -1, dtype=int)
    rows = []

    if len(df):
        tree = KDTree(coords)
        for pos in np.flatnonzero(is_bt):
            if labels[pos] != -1:
                continue
            neighbors = np.asarray(tree.query_ball_point(coords[pos], r=radius), dtype=int)
            free = neighbors[labels[neighbors] == -1]
            n_bt = int(is_bt[free].sum())
            if n_bt >= min_cells:
                agg_id = len(rows)
                labels[free] = agg_id
                rows.append(
                    {
                        "aggregate_id": agg_id,
                        "center_index": df.index[pos],
                        "center_x": coords[pos, 0],
                        "center_y": coords[pos, 1],
                        "classification": LYMPHOID_AGGREGATE,
                        "n_cells": int(free.size),
                        "n_b_cells": int(is_b[free].sum()),
                        "n_t_cells": int(is_t[free].sum()),
                    }
                )

    columns = [
        "aggregate_id",
        "center_index",
        "center_x",
        "center_y",
        "classification",
        "n_cells",
        "n_b_cells",
        "n_t_cells",
    ]
    result = pd.DataFrame(rows, columns=columns)
    if return_labels:
        return result, pd.Series(labels, index=df.index, name="aggregate_id")
    return result


def classify_cells(
    df: pd.DataFrame,
    radius: float = 80.0,
    min_cells: int = 50,
    b_labels: str | Sequence[str] = DEFAULT_B_LABELS,
    t_labels: str | Sequence[str] = DEFAULT_T_LABELS,
    hev_labels: str | Sequence[str] = DEFAULT_HEV_LABELS,
    fdc_labels: str | Sequence[str] = DEFAULT_FDC_LABELS,
    x_col: str = "X",
    y_col: str = "Y",
    phenotype_col: str = "Phenotype",
) -> pd.Series:
    """Classify every cell by the composition of its neighbourhood.

    For each cell, the B and T cells within ``radius`` (the cell itself
    included) are counted, and the presence of a high endothelial venule
    (HEV) and/or a follicular dendritic cell (fDC) within ``radius`` is
    checked:

    * at least ``min_cells`` B + T cells **and** an HEV and/or an fDC:
      ``"Primary Follicle-like TLS"``
    * at least ``min_cells`` B + T cells, neither HEV nor fDC:
      ``"Lymphoid Aggregate"``
    * otherwise: ``"Not TLS"``

    Neighbourhoods may overlap, so this is a per-cell map rather than a list
    of discrete structures. Use :func:`find_lymphoid_aggregates` for discrete
    aggregates.

    Returns
    -------
    pandas.Series
        Classification per cell, indexed like ``df``, named
        ``TLS_Classification``.
    """
    validate_input(df, x_col, y_col, phenotype_col)
    _check_params(radius, min_cells)
    b_labels = _as_tuple(b_labels)
    t_labels = _as_tuple(t_labels)
    hev_labels = _as_tuple(hev_labels)
    fdc_labels = _as_tuple(fdc_labels)

    coords = df[[x_col, y_col]].to_numpy(dtype=float)
    phenotypes = df[phenotype_col].to_numpy()
    is_bt = np.isin(phenotypes, b_labels + t_labels)
    is_maturity = np.isin(phenotypes, hev_labels + fdc_labels)

    def _counts(mask: np.ndarray) -> np.ndarray:
        if not mask.any():
            return np.zeros(len(df), dtype=int)
        return np.asarray(
            KDTree(coords[mask]).query_ball_point(coords, r=radius, return_length=True),
            dtype=int,
        )

    if len(df):
        n_bt = _counts(is_bt)
        has_maturity = _counts(is_maturity) > 0
    else:
        n_bt = np.zeros(0, dtype=int)
        has_maturity = np.zeros(0, dtype=bool)

    dense = n_bt >= min_cells
    out = np.where(
        dense & has_maturity,
        PRIMARY_FOLLICLE_LIKE_TLS,
        np.where(dense, LYMPHOID_AGGREGATE, NOT_TLS),
    )
    return pd.Series(out, index=df.index, name="TLS_Classification", dtype=object)


# ---------------------------------------------------------------------------
# Backwards-compatible wrappers matching the original notebook functions.
# ---------------------------------------------------------------------------


def measure_TLS(
    df: pd.DataFrame,
    radius: float,
    b_labels: str | Sequence[str] = ("CD20",),
    t_labels: str | Sequence[str] = ("CD8", "CD3"),
    min_cells: int = 50,
) -> pd.DataFrame:
    """Original notebook interface for :func:`find_lymphoid_aggregates`.

    Returns a DataFrame with the original column names
    ``Aggregate_Center_Index``, ``Aggregate_Center_Coordinates`` and
    ``Classification``.
    """
    res = find_lymphoid_aggregates(
        df, radius=radius, min_cells=min_cells, b_labels=b_labels, t_labels=t_labels
    )
    return pd.DataFrame(
        {
            "Aggregate_Center_Index": res["center_index"].tolist(),
            "Aggregate_Center_Coordinates": list(zip(res["center_x"], res["center_y"])),
            "Classification": res["classification"].tolist(),
        }
    )


def measure_TLS_HEV(
    df: pd.DataFrame,
    radius: float,
    b_labels: str | Sequence[str] = ("B-cell",),
    t_labels: str | Sequence[str] = ("T-cell",),
    hev_labels: str | Sequence[str] = ("HEV",),
    fdc_labels: str | Sequence[str] = ("fDC",),
    min_cells: int = 50,
) -> pd.DataFrame:
    """Original notebook interface for :func:`classify_cells`.

    Returns a copy of ``df`` with a ``TLS_Classification`` column added.
    Unlike the notebook version, the input DataFrame is not modified.
    """
    out = df.copy()
    out["TLS_Classification"] = classify_cells(
        df,
        radius=radius,
        min_cells=min_cells,
        b_labels=b_labels,
        t_labels=t_labels,
        hev_labels=hev_labels,
        fdc_labels=fdc_labels,
    )
    return out
