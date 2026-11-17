# Changelog

All notable changes to TLS Finder are documented here. This project uses
[semantic versioning](https://semver.org/).

## [Unreleased]

## [1.0.0]

First packaged release. The detection algorithm is unchanged from the 2024
notebooks and is checked against them by regression tests.

### Added
- Installable Python package `tls_finder` (`pip install`), with a
  `tls-finder` command-line tool.
- `find_lymphoid_aggregates()`: aggregate detection with configurable radius,
  cell-count threshold, phenotype labels and column names; reports B/T cell
  counts per aggregate and optional per-cell aggregate membership.
- `classify_cells()`: per-cell classification (Lymphoid Aggregate /
  Primary Follicle-like TLS / Not TLS). A neighbourhood is Primary
  Follicle-like when an HEV **and/or** a follicular dendritic cell (fDC) is
  present (configurable `hev_labels`, `fdc_labels`); the notebooks used HEV only.
- `make_synthetic_tissue()` for examples and testing.
- Plotting helpers: `plot_phenotypes`, `plot_aggregates`, `plot_cell_classes`.
- Input validation with clear error messages.
- Test suite and continuous integration on Linux, macOS and Windows.
- `CONTRIBUTING.md`, `CHANGELOG.md`; citation file renamed to `CITATION.cff`
  so GitHub detects it.

### Changed
- README method description corrected: TLS Finder evaluates each candidate
  cell's one-hop neighbourhood in the radius proximity graph, not connected
  components.
- Results no longer depend on the DataFrame having a default 0..n-1 index.
- `measure_TLS_HEV` no longer modifies the input DataFrame in place.
- Per-cell classification is vectorised and aggregate detection skips the
  unused NetworkX graph construction, making both much faster on large
  images (≈120,000 cells in under a second).

## [0.1.0] - 2024-12-26

- Initial public release as Jupyter notebooks accompanying the bioRxiv
  preprint (doi:10.1101/2024.12.26.630405).
