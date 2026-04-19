import numpy as np
import pandas as pd

from qpcr_analyzer.core.outliers import _flag_cluster, mark_outliers


def test_single_well_is_flagged():
    assert _flag_cluster(np.array([20.0]), 1.0).tolist() == [True]


def test_pair_close_kept():
    assert _flag_cluster(np.array([20.0, 20.5]), 1.0).tolist() == [False, False]


def test_pair_far_flagged():
    assert _flag_cluster(np.array([20.0, 22.0]), 1.0).tolist() == [True, True]


def test_triplet_one_outlier():
    flags = _flag_cluster(np.array([20.0, 20.2, 22.0]), 1.0).tolist()
    assert flags == [False, False, True]


def test_triplet_all_different_all_flagged():
    flags = _flag_cluster(np.array([20.0, 22.0, 24.0]), 1.0).tolist()
    assert flags == [True, True, True]


def test_triplet_prefers_tightest_cluster():
    # [20, 20.5, 21.5] has two length-2 runs; the tighter one wins.
    flags = _flag_cluster(np.array([20.0, 20.5, 21.5]), 1.0).tolist()
    assert flags == [False, False, True]


def test_five_replicates():
    flags = _flag_cluster(np.array([20.0, 20.1, 20.2, 22.0, 24.0]), 1.0).tolist()
    assert flags == [False, False, False, True, True]


def test_quadruplet_symmetric_is_ambiguous():
    # [20, 21, 22, 23] with tol=1 has runs of length 2 with equal range → flag all
    flags = _flag_cluster(np.array([20.0, 21.0, 22.0, 23.0]), 1.0).tolist()
    assert flags == [True, True, True, True]


def test_mark_outliers_with_nan():
    df = pd.DataFrame(
        {
            "Well": ["A1", "A2", "A3"],
            "Target": ["g", "g", "g"],
            "Sample": ["s1", "s1", "s1"],
            "Group": ["all", "all", "all"],
            "Cq": [20.0, 20.1, float("nan")],
        }
    )
    out = mark_outliers(df, tolerance=1.0)
    by_well = dict(zip(out["Well"], out["Outlier"]))
    assert by_well["A1"] == False  # noqa: E712
    assert by_well["A2"] == False  # noqa: E712
    assert by_well["A3"] == True  # noqa: E712
    assert out["Replicates"].iloc[0] == 2


def test_mark_outliers_variable_replicates():
    # Sample s1 has 4 replicates with one clear outlier, s2 has 2 matching
    df = pd.DataFrame(
        {
            "Well": ["A1", "A2", "A3", "A4", "B1", "B2"],
            "Target": ["g"] * 6,
            "Sample": ["s1", "s1", "s1", "s1", "s2", "s2"],
            "Group": ["all"] * 6,
            "Cq": [20.0, 20.1, 20.2, 23.0, 25.0, 25.3],
        }
    )
    out = mark_outliers(df, tolerance=1.0)
    by_well = dict(zip(out["Well"], out["Outlier"]))
    assert by_well["A4"] == True  # noqa: E712
    assert by_well["A1"] == False  # noqa: E712
    assert by_well["B1"] == False  # noqa: E712
    assert by_well["B2"] == False  # noqa: E712
