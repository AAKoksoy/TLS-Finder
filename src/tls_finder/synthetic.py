"""Synthetic tissue generator for examples, tutorials and tests."""

from __future__ import annotations

import numpy as np
import pandas as pd


def _ring(rng, cx, cy, r_min, r_max, n):
    radii = rng.uniform(r_min, r_max, n)
    angles = rng.uniform(0, 2 * np.pi, n)
    return cx + radii * np.cos(angles), cy + radii * np.sin(angles)


def make_synthetic_tissue(
    n_background: int = 5000,
    size: float = 10000.0,
    tls_centers: list[tuple[float, float]] | None = None,
    n_b_cells: int = 3000,
    n_t_cells: int = 2000,
    b_ring: tuple[float, float] = (100.0, 400.0),
    t_ring: tuple[float, float] = (400.0, 600.0),
    with_hev: bool = True,
    background_labels: tuple[str, ...] = ("Other", "CD20", "CD3"),
    b_label: str = "CD20",
    t_label: str = "CD3",
    hev_label: str = "HEV",
    seed: int | None = 0,
) -> pd.DataFrame:
    """Generate a square tissue with idealised TLS embedded in random cells.

    Each TLS is a central HEV (optional) surrounded by a ring of B cells,
    itself surrounded by a ring of T cells, as in the original TLS Finder
    example notebook. Background cells are scattered uniformly.

    Returns
    -------
    pandas.DataFrame
        Columns ``X``, ``Y``, ``Phenotype``.
    """
    rng = np.random.default_rng(seed)
    if tls_centers is None:
        tls_centers = [(size / 2, size / 2)]

    xs, ys, labels = [], [], []
    for cx, cy in tls_centers:
        if with_hev:
            xs.append(np.array([cx]))
            ys.append(np.array([cy]))
            labels.append(np.array([hev_label]))
        bx, by = _ring(rng, cx, cy, *b_ring, n_b_cells)
        tx, ty = _ring(rng, cx, cy, *t_ring, n_t_cells)
        xs += [bx, tx]
        ys += [by, ty]
        labels += [np.full(n_b_cells, b_label), np.full(n_t_cells, t_label)]

    xs.append(rng.uniform(0, size, n_background))
    ys.append(rng.uniform(0, size, n_background))
    labels.append(rng.choice(np.asarray(background_labels), n_background))

    return pd.DataFrame(
        {
            "X": np.concatenate(xs),
            "Y": np.concatenate(ys),
            "Phenotype": np.concatenate(labels).astype(object),
        }
    )
