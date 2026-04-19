import pandas as pd
import pytest

from qpcr_analyzer.core.quant import compute_relative_expression


def _triplicate(values):
    return values * 3


def test_single_group_single_ref_sample():
    df = pd.DataFrame(
        {
            "Well": [f"W{i}" for i in range(12)],
            "Target": ["gapdh"] * 3 + ["geneX"] * 3 + ["gapdh"] * 3 + ["geneX"] * 3,
            "Sample": ["s1"] * 6 + ["s2"] * 6,
            "Group": ["all"] * 12,
            "Cq": [
                18.0, 18.1, 18.05,
                22.0, 22.1, 22.05,
                18.0, 18.1, 18.05,
                21.0, 21.1, 21.05,
            ],
            "Excluded": [False] * 12,
        }
    )
    res = compute_relative_expression(df, ["gapdh"], {"all": ["s1"]})
    r = res["gapdh"]
    s1 = r[r["Sample"] == "s1"]["Relative_Expr"].iloc[0]
    s2 = r[r["Sample"] == "s2"]["Relative_Expr"].iloc[0]
    assert abs(s1 - 1.0) < 1e-6
    assert abs(s2 - 2.0) < 0.05  # 2^~1


def test_per_group_normalization():
    # Each group has its own reference sample; relative expression is
    # normalized within the group.
    df = pd.DataFrame(
        {
            "Well": [f"W{i}" for i in range(12)],
            "Target": ["gapdh", "geneX"] * 6,
            "Sample": sum(([f"s{i}"] * 2 for i in [1, 2, 3, 4, 5, 6]), []),
            "Group": ["g1"] * 6 + ["g2"] * 6,
            "Cq": [18, 22, 18, 21, 18, 23, 20, 24, 20, 24, 20, 26],
            "Excluded": [False] * 12,
        }
    )
    res = compute_relative_expression(
        df, ["gapdh"], {"g1": ["s1"], "g2": ["s4"]}
    )
    r = res["gapdh"]
    by_sample = dict(zip(r["Sample"], r["Relative_Expr"]))
    assert abs(by_sample["s1"] - 1.0) < 1e-6
    assert abs(by_sample["s2"] - 2.0) < 1e-6
    assert abs(by_sample["s3"] - 0.5) < 1e-6
    assert abs(by_sample["s4"] - 1.0) < 1e-6
    assert abs(by_sample["s5"] - 1.0) < 1e-6
    assert abs(by_sample["s6"] - 0.25) < 1e-6


def test_missing_ref_raises():
    df = pd.DataFrame(
        {
            "Well": ["A1", "A2"],
            "Target": ["geneX", "geneX"],
            "Sample": ["s1", "s1"],
            "Group": ["all", "all"],
            "Cq": [22.0, 22.1],
            "Excluded": [False, False],
        }
    )
    with pytest.raises(ValueError):
        compute_relative_expression(df, ["gapdh"], {"all": ["s1"]})


def test_group_without_reference_raises():
    df = pd.DataFrame(
        {
            "Well": [f"W{i}" for i in range(8)],
            "Target": ["gapdh", "geneX"] * 4,
            "Sample": ["s1"] * 2 + ["s2"] * 2 + ["s3"] * 2 + ["s4"] * 2,
            "Group": ["g1"] * 4 + ["g2"] * 4,
            "Cq": [18, 22, 18, 21, 20, 24, 20, 23],
            "Excluded": [False] * 8,
        }
    )
    with pytest.raises(ValueError):
        compute_relative_expression(df, ["gapdh"], {"g1": ["s1"]})


def test_multiple_ref_genes():
    df = pd.DataFrame(
        {
            "Well": [f"W{i}" for i in range(9)],
            "Target": ["gapdh", "actb", "geneX"] * 3,
            "Sample": ["s1"] * 3 + ["s2"] * 3 + ["s3"] * 3,
            "Group": ["all"] * 9,
            "Cq": [18, 19, 22, 18, 19, 21, 18, 19, 23],
            "Excluded": [False] * 9,
        }
    )
    res = compute_relative_expression(
        df, ["gapdh", "actb"], {"all": ["s1"]}
    )
    assert set(res.keys()) == {"gapdh", "actb"}
    # s1 is reference, so its RelExpr = 1 for both refs
    for ref, r in res.items():
        s1_val = r[r["Sample"] == "s1"]["Relative_Expr"].iloc[0]
        assert abs(s1_val - 1.0) < 1e-6
