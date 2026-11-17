import matplotlib

matplotlib.use("Agg")

import pandas as pd

from tls_finder import classify_cells, find_lymphoid_aggregates, make_synthetic_tissue
from tls_finder.cli import main
from tls_finder.plotting import plot_aggregates, plot_cell_classes, plot_phenotypes


def _tissue():
    return make_synthetic_tissue(
        n_background=500, size=3000, tls_centers=[(1500, 1500)],
        n_b_cells=300, n_t_cells=200, seed=4,
    )


def test_cli_aggregates(tmp_path, capsys):
    inp, out, fig = tmp_path / "cells.csv", tmp_path / "agg.csv", tmp_path / "agg.png"
    _tissue().to_csv(inp, index=False)
    assert main([str(inp), "-r", "160", "-o", str(out), "--plot", str(fig)]) == 0
    res = pd.read_csv(out)
    assert len(res) >= 1 and "center_x" in res.columns
    assert fig.stat().st_size > 0
    assert "lymphoid aggregate" in capsys.readouterr().err


def test_cli_cells_mode(tmp_path):
    inp, out = tmp_path / "cells.csv", tmp_path / "cls.csv"
    _tissue().to_csv(inp, index=False)
    assert main([str(inp), "--mode", "cells", "-r", "150", "-o", str(out)]) == 0
    res = pd.read_csv(out)
    assert "TLS_Classification" in res.columns
    assert (res["TLS_Classification"] == "Primary Follicle-like TLS").any()


def test_cli_reports_bad_input(tmp_path, capsys):
    inp = tmp_path / "bad.csv"
    pd.DataFrame({"X": [1], "Y": [2]}).to_csv(inp, index=False)
    assert main([str(inp)]) == 2
    assert "missing required column" in capsys.readouterr().err


def test_plots_return_axes():
    df = _tissue()
    agg = find_lymphoid_aggregates(df)
    classes = classify_cells(df, radius=150)
    for ax in (
        plot_phenotypes(df),
        plot_phenotypes(df, colors={"CD20": "blue", "CD3": "green"}),
        plot_aggregates(df, agg),
        plot_cell_classes(df, classes),
    ):
        assert ax.get_xlabel() == "X (microns)"
        assert ax.yaxis_inverted()
