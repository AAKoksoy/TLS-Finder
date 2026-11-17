import numpy as np
import pandas as pd
import pytest

from tls_finder import (
    LYMPHOID_AGGREGATE,
    NOT_TLS,
    PRIMARY_FOLLICLE_LIKE_TLS,
    classify_cells,
    find_lymphoid_aggregates,
    make_synthetic_tissue,
    measure_TLS,
    measure_TLS_HEV,
)

from reference import measure_TLS_HEV_original, measure_TLS_original


@pytest.fixture
def tissue():
    # Two TLS in a sparse background; small enough for the slow reference code.
    return make_synthetic_tissue(
        n_background=1500,
        size=4000,
        tls_centers=[(1000, 1000), (3000, 2800)],
        n_b_cells=300,
        n_t_cells=200,
        background_labels=("Other", "CD20", "CD8"),
        b_label="CD20",
        t_label="CD8",
        seed=1,
    )


# --- find_lymphoid_aggregates ------------------------------------------------


def test_detects_both_synthetic_tls(tissue):
    res = find_lymphoid_aggregates(tissue, radius=160)
    assert len(res) >= 2
    centers = res[["center_x", "center_y"]].to_numpy()
    for cx, cy in [(1000, 1000), (3000, 2800)]:
        d = np.hypot(centers[:, 0] - cx, centers[:, 1] - cy)
        assert d.min() < 700, "no aggregate found near a synthetic TLS"


def test_no_aggregates_in_sparse_tissue():
    df = make_synthetic_tissue(n_background=2000, size=10000, tls_centers=[], seed=3)
    res = find_lymphoid_aggregates(df, radius=160)
    assert res.empty
    assert list(res.columns)[:5] == [
        "aggregate_id", "center_index", "center_x", "center_y", "classification",
    ]


@pytest.mark.parametrize("radius", [120, 160, 250])
def test_matches_original_notebook_algorithm(tissue, radius):
    expected = measure_TLS_original(tissue, radius)
    got = measure_TLS(tissue, radius, b_labels=["CD20"], t_labels=["CD8"])
    pd.testing.assert_frame_equal(got, expected, check_dtype=False)


def test_aggregates_do_not_share_cells(tissue):
    res, labels = find_lymphoid_aggregates(tissue, radius=160, return_labels=True)
    assert labels.index.equals(tissue.index)
    counts = labels[labels >= 0].value_counts().sort_index()
    assert counts.tolist() == res["n_cells"].tolist()
    assert (res["n_b_cells"] + res["n_t_cells"] >= 50).all()


def test_non_default_index_is_handled(tissue):
    shuffled = tissue.sample(frac=1, random_state=0)
    shuffled.index = [f"cell_{i}" for i in range(len(shuffled))]
    res = find_lymphoid_aggregates(shuffled, radius=160)
    assert res["center_index"].isin(shuffled.index).all()
    reset = find_lymphoid_aggregates(shuffled.reset_index(drop=True), radius=160)
    assert len(res) == len(reset)


def test_input_not_modified(tissue):
    before = tissue.copy()
    find_lymphoid_aggregates(tissue, radius=160)
    classify_cells(tissue, radius=80)
    pd.testing.assert_frame_equal(tissue, before)


def test_custom_column_names(tissue):
    renamed = tissue.rename(columns={"X": "x_um", "Y": "y_um", "Phenotype": "cell_type"})
    a = find_lymphoid_aggregates(tissue, radius=160)
    b = find_lymphoid_aggregates(
        renamed, radius=160, x_col="x_um", y_col="y_um", phenotype_col="cell_type"
    )
    pd.testing.assert_frame_equal(a, b)


def test_min_cells_threshold():
    # 10 B cells on top of each other: aggregate at min_cells=10, not at 11.
    df = pd.DataFrame({"X": [0.0] * 10, "Y": [0.0] * 10, "Phenotype": ["CD20"] * 10})
    assert len(find_lymphoid_aggregates(df, radius=1, min_cells=10)) == 1
    assert len(find_lymphoid_aggregates(df, radius=1, min_cells=11)) == 0


def test_empty_input():
    df = pd.DataFrame({"X": [], "Y": [], "Phenotype": []}, dtype=float)
    assert find_lymphoid_aggregates(df).empty
    assert classify_cells(df).empty


@pytest.mark.parametrize(
    "kwargs",
    [{"radius": 0}, {"radius": -5}, {"radius": float("nan")}, {"min_cells": 0}, {"min_cells": 2.5}],
)
def test_invalid_parameters(tissue, kwargs):
    with pytest.raises(ValueError):
        find_lymphoid_aggregates(tissue, **kwargs)


def test_missing_column():
    with pytest.raises(ValueError, match="missing required column"):
        find_lymphoid_aggregates(pd.DataFrame({"X": [1.0], "Y": [1.0]}))


def test_non_numeric_or_nan_coordinates():
    with pytest.raises(ValueError, match="numeric"):
        find_lymphoid_aggregates(pd.DataFrame({"X": ["a"], "Y": [1.0], "Phenotype": ["CD20"]}))
    with pytest.raises(ValueError, match="missing values"):
        find_lymphoid_aggregates(pd.DataFrame({"X": [np.nan], "Y": [1.0], "Phenotype": ["CD20"]}))


def test_rejects_non_dataframe():
    with pytest.raises(TypeError):
        find_lymphoid_aggregates([[0, 0, "CD20"]])


# --- classify_cells ----------------------------------------------------------


@pytest.fixture
def hev_tissue():
    return make_synthetic_tissue(
        n_background=1500,
        size=4000,
        tls_centers=[(2000, 2000)],
        n_b_cells=600,
        n_t_cells=400,
        background_labels=("Others", "B-cell", "T-cell"),
        b_label="B-cell",
        t_label="T-cell",
        seed=2,
    )


@pytest.mark.parametrize("radius", [60, 80, 150])
def test_classify_matches_original_notebook(hev_tissue, radius):
    expected = measure_TLS_HEV_original(hev_tissue, radius)
    got = measure_TLS_HEV(hev_tissue, radius)
    assert got["TLS_Classification"].tolist() == expected


def test_classify_finds_hev_tls(hev_tissue):
    classes = classify_cells(
        hev_tissue, radius=150, b_labels="B-cell", t_labels="T-cell"
    )
    hev_pos = np.flatnonzero(hev_tissue["Phenotype"].to_numpy() == "HEV")[0]
    assert classes.iloc[hev_pos] == PRIMARY_FOLLICLE_LIKE_TLS
    assert set(classes.unique()) <= {PRIMARY_FOLLICLE_LIKE_TLS, LYMPHOID_AGGREGATE, NOT_TLS}
    assert (classes == NOT_TLS).any()


def test_classify_without_any_hev(hev_tissue):
    no_hev = hev_tissue[hev_tissue["Phenotype"] != "HEV"]
    classes = classify_cells(no_hev, radius=150, b_labels="B-cell", t_labels="T-cell")
    assert not (classes == PRIMARY_FOLLICLE_LIKE_TLS).any()
    assert (classes == LYMPHOID_AGGREGATE).any()


def test_fdc_alone_marks_primary_follicle_like_tls(hev_tissue):
    with_fdc = hev_tissue.copy()
    with_fdc.loc[with_fdc["Phenotype"] == "HEV", "Phenotype"] = "fDC"
    classes = classify_cells(with_fdc, radius=150, b_labels="B-cell", t_labels="T-cell")
    centre = np.flatnonzero(with_fdc["Phenotype"].to_numpy() == "fDC")[0]
    assert classes.iloc[centre] == PRIMARY_FOLLICLE_LIKE_TLS


def test_hev_and_fdc_together_same_as_either(hev_tissue):
    both = pd.concat(
        [hev_tissue, pd.DataFrame({"X": [2000.0], "Y": [2000.0], "Phenotype": ["fDC"]})],
        ignore_index=True,
    )
    a = classify_cells(both, radius=150, b_labels="B-cell", t_labels="T-cell")
    hev_pos = np.flatnonzero(both["Phenotype"].to_numpy() == "HEV")[0]
    assert a.iloc[hev_pos] == PRIMARY_FOLLICLE_LIKE_TLS


def test_custom_fdc_labels(hev_tissue):
    df = hev_tissue[hev_tissue["Phenotype"] != "HEV"].copy()
    df = pd.concat([df, pd.DataFrame({"X": [2000.0], "Y": [2000.0], "Phenotype": ["CD21+"]})],
                   ignore_index=True)
    default = classify_cells(df, radius=150, b_labels="B-cell", t_labels="T-cell")
    custom = classify_cells(df, radius=150, b_labels="B-cell", t_labels="T-cell",
                            fdc_labels=["CD21+"])
    assert not (default == PRIMARY_FOLLICLE_LIKE_TLS).any()
    assert (custom == PRIMARY_FOLLICLE_LIKE_TLS).any()
