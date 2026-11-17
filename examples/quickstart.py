"""TLS Finder quick start on synthetic data.

Run with:  python examples/quickstart.py
Requires:  pip install "tls-finder[plot]"
"""

import matplotlib.pyplot as plt

from tls_finder import classify_cells, find_lymphoid_aggregates, make_synthetic_tissue
from tls_finder.plotting import plot_aggregates, plot_cell_classes, plot_phenotypes

# Three idealised TLS (central HEV, B-cell ring, T-cell ring) in a sparse background.
cells = make_synthetic_tissue(
    n_background=20000,
    size=10000,
    tls_centers=[(2500, 2500), (7000, 3000), (5000, 7500)],
    seed=42,
)
print(cells["Phenotype"].value_counts(), "\n")

aggregates = find_lymphoid_aggregates(
    cells, radius=160, b_labels=["CD20"], t_labels=["CD3"]
)
print(f"{len(aggregates)} lymphoid aggregates found")
print(aggregates.head(), "\n")

classes = classify_cells(
    cells, radius=80, b_labels=["CD20"], t_labels=["CD3"], hev_labels=["HEV"]
)
print(classes.value_counts())

fig, axes = plt.subplots(1, 3, figsize=(20, 6))
plot_phenotypes(cells, colors={"Other": "lightgray", "CD20": "blue", "CD3": "green", "HEV": "red"}, ax=axes[0])
plot_aggregates(cells, aggregates, ax=axes[1])
plot_cell_classes(cells, classes, ax=axes[2])
fig.tight_layout()
fig.savefig("tls_finder_quickstart.png", dpi=150)
print("\nFigure saved to tls_finder_quickstart.png")
