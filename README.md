<img width="1254" height="1254" alt="image" src="https://github.com/user-attachments/assets/1da54564-1b3f-4dfa-978e-78d598527fac" />

# TLS Finder

[![tests](https://github.com/AAKoksoy/TLS-Finder/actions/workflows/tests.yml/badge.svg)](https://github.com/AAKoksoy/TLS-Finder/actions/workflows/tests.yml)

**TLS Finder** is an open-source Python package for automated identification of tertiary lymphoid structures (TLS) from spatial immune-cell coordinate data.

📄 **Original publication:** Koksoy et al., bioRxiv (2024)

🔗 **DOI:** https://doi.org/10.1101/2024.12.26.630405

📅 **First public release:** December 26, 2024

## 📖 Citation

If you use TLS Finder in your research, please cite the original
publication:

> Koksoy, A.A., Salvatierra, M.E., MOSAIC Team, Solis Soto, L.M.,
> Haymaker, C. (2024). *TLS_Finder: An algorithm for Identifying
> Tertiary Lymphoid Structures Using Immune Cell Spatial
> Coordinates.* bioRxiv.
> https://doi.org/10.1101/2024.12.26.630405

A `CITATION.cff` file is included in this repository — GitHub will
show a **"Cite this repository"** button on the repo's main page
(top right, under "About") that generates both the software citation
and the article citation above in APA and BibTeX formats
automatically.

**BibTeX (article):**

```bibtex
@article{koksoy2024tlsfinder,
  author  = {Koksoy, Ayse A. and Salvatierra, Maria Esther and {MOSAIC Team} and Solis Soto, Luisa Maren and Haymaker, Cara},
  title   = {TLS\_Finder: An algorithm for Identifying Tertiary Lymphoid Structures Using Immune Cell Spatial Coordinates},
  journal = {bioRxiv},
  year    = {2024},
  doi     = {10.1101/2024.12.26.630405},
  url     = {https://doi.org/10.1101/2024.12.26.630405}
}
```

**BibTeX (software):**

```bibtex
@software{koksoy_tlsfinder_software,
  author  = {Koksoy, Ayse A.},
  title   = {TLS Finder},
  version = {1.0.0},
  year    = {2026},
  url     = {https://github.com/AAKoksoy/TLS-Finder}
}
```

## Project Information

| Item | Details |
|------|---------|
| **Published** | December 2024 (bioRxiv) |
| **DOI** | https://doi.org/10.1101/2024.12.26.630405 |
| **First Public Release** | December 26, 2024 |
| **Current version** | 1.0.0 (see [CHANGELOG](CHANGELOG.md)) |
| **Purpose** | Automated detection of tertiary lymphoid structures from spatial immune-cell coordinates |
| **Status** | Stable |

TLS Finder is part of the MOSAIC ecosystem for graph-native spatial biology.

---

## 🔬 Overview

TLS Finder identifies **tertiary lymphoid structures (TLS)** and lymphoid
aggregates from phenotyped single-cell coordinates produced by spatial
platforms such as multiplex immunofluorescence (mIF), CODEX/PhenoCycler,
imaging mass cytometry, or cell-segmented spatial transcriptomics.

It needs only a table of cell coordinates and phenotype labels, so it runs
downstream of any segmentation and phenotyping pipeline.

## 📦 Installation

```bash
pip install "tls-finder[plot] @ git+https://github.com/AAKoksoy/TLS-Finder"
```

Requires Python ≥ 3.10 with NumPy, pandas and SciPy; matplotlib is needed only
for plotting (`[plot]`).

## 🚀 Quick start

```python
import pandas as pd
from tls_finder import find_lymphoid_aggregates, classify_cells
from tls_finder.plotting import plot_aggregates

cells = pd.read_csv("cells.csv")          # columns: X, Y, Phenotype

# 1. Discrete, non-overlapping lymphoid aggregates
aggregates = find_lymphoid_aggregates(
    cells, radius=160, min_cells=50,
    b_labels=["CD20"], t_labels=["CD3", "CD8"],
)
print(aggregates[["center_x", "center_y", "n_b_cells", "n_t_cells"]])
plot_aggregates(cells, aggregates)

# 2. Per-cell classification, including HEV- and/or fDC-associated TLS
cells["TLS_Classification"] = classify_cells(
    cells, radius=80, b_labels=["CD20"], t_labels=["CD3", "CD8"], hev_labels=["HEV"], fdc_labels=["fDC"],
)
```

No data at hand? `tls_finder.make_synthetic_tissue()` generates a tissue with
idealised TLS (central HEV, B-cell ring, T-cell ring) embedded in random cells.
See [`examples/quickstart.py`](examples/quickstart.py).

### Command line

```bash
tls-finder cells.csv --radius 160 -o aggregates.csv --plot aggregates.png
tls-finder cells.csv --mode cells --radius 80 --hev-labels HEV --fdc-labels fDC -o classified.csv
tls-finder --help
```

## ⚙️ Method

**`find_lymphoid_aggregates`** (the algorithm in Koksoy et al., 2024):

1. Cells are treated as nodes of a spatial proximity graph in which two
   cells are linked if they lie within `radius` of each other. A cell's
   neighbourhood is its node plus all directly linked nodes (one hop).
   Neighbourhoods are retrieved with a KD-tree radius query, which returns
   exactly this set without building the full graph in memory.
2. Cells are visited in input order. Each B or T cell that is not yet part of
   an aggregate is a candidate centre.
3. If at least `min_cells` (default 50) B + T cells that are not yet part of
   an aggregate lie within `radius` of the candidate (itself included), the
   candidate is recorded as a **lymphoid aggregate** centre.
4. All cells in that neighbourhood, of any phenotype, are assigned to the
   aggregate and excluded from later candidates and counts, so aggregates
   never share cells.

**`classify_cells`** labels every cell by its own neighbourhood:
`Primary Follicle-like TLS` if it has ≥ `min_cells` B + T cells and at least
one high endothelial venule (HEV) **and/or** follicular dendritic cell (fDC)
within `radius`; `Lymphoid Aggregate` if it has the B + T density but neither
HEV nor fDC; otherwise `Not TLS`. Neighbourhoods overlap, so this is a per-cell
map rather than a list of structures.

**Choosing parameters.** `radius` and `min_cells` together set a density
threshold (B + T cells per area), so appropriate values depend on tissue
type, imaging resolution and phenotyping depth. Defaults follow the original
publication (radius 160 µm for aggregates; 80 µm for per-cell
classification). Inspect results visually with the plotting helpers when
applying TLS Finder to a new dataset.

## 📂 Input

A pandas DataFrame (or CSV for the command line) with one row per cell:

| Column      | Description                                  |
| ----------- | -------------------------------------------- |
| `X`         | x coordinate (same unit as `radius`, e.g. µm) |
| `Y`         | y coordinate                                 |
| `Phenotype` | cell phenotype label (e.g. `CD20`, `CD8`, `HEV`, `fDC`) |

Column names and phenotype labels are configurable (`x_col`, `y_col`,
`phenotype_col`, `b_labels`, `t_labels`, `hev_labels`, `fdc_labels`).

## 📊 Output

`find_lymphoid_aggregates` returns one row per aggregate:

| Column | Description |
| ------ | ----------- |
| `aggregate_id` | 0-based aggregate number |
| `center_index` | index label of the centre cell in the input |
| `center_x`, `center_y` | centre coordinates |
| `classification` | `Lymphoid Aggregate` |
| `n_cells`, `n_b_cells`, `n_t_cells` | cells assigned to the aggregate |

With `return_labels=True` it also returns a per-cell `aggregate_id` Series
(`-1` = not in an aggregate).

`classify_cells` returns a per-cell `TLS_Classification` Series.

The original notebook functions `measure_TLS` and `measure_TLS_HEV` remain
available with their original output format.

## 📓 Notebooks

The original analysis notebooks are kept in this repository for reference.
New work should use the installable package above.

## 🧪 Tests

```bash
pip install -e ".[test]"
pytest
```

The test suite includes regression tests that check the package reproduces
the results of the original 2024 notebook implementation exactly.

## 🤝 Contributing

Contributions, bug reports and questions are welcome — see
[CONTRIBUTING.md](CONTRIBUTING.md).

## 📜 License

© 2024–2026 Ayse A. Koksoy

This project is licensed under the Apache License 2.0.
See the [LICENSE](LICENSE) file for details.

## 👩‍🔬 Author

**Ayse A. Koksoy, MD, PhD**
Computational scientist in spatial omics and multiplex imaging

"I build AI-assisted pipelines that integrate spatial, transcriptomic, and imaging data into interpretable biological insights."
