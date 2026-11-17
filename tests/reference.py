"""Frozen copies of the original TLS Finder notebook functions (December 2024).

Used only as a regression reference: the packaged implementation must give the
same results. Unused graph construction and unused variables are removed; the
decision logic is unchanged. Do not edit.
"""

import pandas as pd
from scipy.spatial import KDTree


def measure_TLS_original(df, radius, bt_labels=("CD20", "CD8")):
    tree = KDTree(df[["X", "Y"]].values)
    used_cells = set()
    tls_classifications, aggregate_centers, aggregate_coordinates = [], [], []
    for idx, row in df.iterrows():
        if row["Phenotype"] in bt_labels and idx not in used_cells:
            neighbors = tree.query_ball_point([row["X"], row["Y"]], r=radius)
            bt_cell_count = sum(
                1 for i in neighbors
                if i not in used_cells and df.at[i, "Phenotype"] in bt_labels
            )
            if bt_cell_count >= 50:
                tls_classifications.append("Lymphoid Aggregate")
                aggregate_centers.append(idx)
                aggregate_coordinates.append((row["X"], row["Y"]))
                used_cells.update(neighbors)
    return pd.DataFrame(
        {
            "Aggregate_Center_Index": aggregate_centers,
            "Aggregate_Center_Coordinates": aggregate_coordinates,
            "Classification": tls_classifications,
        }
    )


def measure_TLS_HEV_original(df, radius):
    tree = KDTree(df[["X", "Y"]].values)
    out = []
    for _, row in df.iterrows():
        neighbors = tree.query_ball_point([row["X"], row["Y"]], r=radius)
        neighbor_phenotypes = df.iloc[neighbors]["Phenotype"]
        bt_cell_count = sum(neighbor_phenotypes.isin(["B-cell", "T-cell"]))
        contains_HEV = any(neighbor_phenotypes.isin(["HEV"]))
        if bt_cell_count >= 50 and contains_HEV:
            out.append("Primary Follicle-like TLS")
        elif bt_cell_count >= 50:
            out.append("Lymphoid Aggregate")
        else:
            out.append("Not TLS")
    return out
