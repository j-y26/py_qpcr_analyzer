import numpy as np
import pandas as pd
import pytest

from qpcr_analyzer.core.quant import (
    compute_delta_ct,
    compute_delta_delta_ct,
    compute_mean_cq,
)


# ── Fixtures ──────────────────────────────────────────────────────────────────

def _make_df(targets, samples, groups, cq_values, excluded=None):
    n = len(cq_values)
    return pd.DataFrame(
        {
            "Well": [f"W{i}" for i in range(n)],
            "Target": targets,
            "Sample": samples,
            "Group": groups,
            "Cq": cq_values,
            "Excluded": excluded if excluded else [False] * n,
        }
    )


# ── compute_delta_ct ──────────────────────────────────────────────────────────

def test_dct_basic():
    """ΔCt = mean_Cq(target) − mean_Cq(HK) per sample."""
    df = _make_df(
        targets=["gapdh", "geneX", "gapdh", "geneX"],
        samples=["s1", "s1", "s2", "s2"],
        groups=["ctrl", "ctrl", "ctrl", "ctrl"],
        cq_values=[18.0, 22.0, 18.0, 21.0],
    )
    res = compute_delta_ct(df, ["gapdh"])
    r = res["gapdh"]
    by_sample = dict(zip(r["Sample"], r["dCt"]))
    assert abs(by_sample["s1"] - 4.0) < 1e-9
    assert abs(by_sample["s2"] - 3.0) < 1e-9


def test_dct_expr_vs_hk():
    """Expr_vs_HK = 2^(−dCt)."""
    df = _make_df(
        targets=["gapdh", "geneX"],
        samples=["s1", "s1"],
        groups=["ctrl", "ctrl"],
        cq_values=[18.0, 20.0],
    )
    res = compute_delta_ct(df, ["gapdh"])
    r = res["gapdh"]
    expected = 2.0 ** -2.0
    assert abs(r["Expr_vs_HK"].iloc[0] - expected) < 1e-9


def test_dct_missing_hk_raises():
    df = _make_df(
        targets=["geneX", "geneX"],
        samples=["s1", "s1"],
        groups=["ctrl", "ctrl"],
        cq_values=[22.0, 22.1],
    )
    with pytest.raises(ValueError, match="gapdh"):
        compute_delta_ct(df, ["gapdh"])


def test_dct_multiple_ref_genes():
    df = _make_df(
        targets=["gapdh", "actb", "geneX"],
        samples=["s1", "s1", "s1"],
        groups=["ctrl", "ctrl", "ctrl"],
        cq_values=[18.0, 19.0, 22.0],
    )
    res = compute_delta_ct(df, ["gapdh", "actb"])
    assert set(res.keys()) == {"gapdh", "actb"}
    # gapdh sheet: targets are actb (dCt=1) and geneX (dCt=4)
    gapdh_dct = dict(zip(res["gapdh"]["Target"], res["gapdh"]["dCt"]))
    assert abs(gapdh_dct["geneX"] - 4.0) < 1e-9
    assert abs(gapdh_dct["actb"] - 1.0) < 1e-9
    # actb sheet: targets are gapdh (dCt=-1) and geneX (dCt=3)
    actb_dct = dict(zip(res["actb"]["Target"], res["actb"]["dCt"]))
    assert abs(actb_dct["geneX"] - 3.0) < 1e-9


def test_dct_excludes_flagged_wells():
    df = _make_df(
        targets=["gapdh", "gapdh", "geneX"],
        samples=["s1", "s1", "s1"],
        groups=["ctrl", "ctrl", "ctrl"],
        cq_values=[18.0, 30.0, 22.0],  # second gapdh is an outlier
        excluded=[False, True, False],
    )
    res = compute_delta_ct(df, ["gapdh"])
    # Mean Cq(gapdh) should be 18.0 only (30.0 excluded)
    assert abs(res["gapdh"]["dCt"].iloc[0] - 4.0) < 1e-9


# ── compute_delta_delta_ct ────────────────────────────────────────────────────

def test_ddct_single_batch_reference_anchor_zero():
    """Reference group samples must have mean ΔΔCt = 0 (anchored).

    The mean of 2^(−ΔΔCt) is NOT guaranteed to equal 1 when ctrl samples
    have different dCt values (due to exponential non-linearity).  The true
    invariant is mean(ΔΔCt_ctrl) = 0, which is enforced by construction.
    """
    df = _make_df(
        targets=["gapdh", "geneX"] * 3,
        samples=["s1", "s1", "s2", "s2", "s3", "s3"],
        groups=["ctrl", "ctrl", "ctrl", "ctrl", "treat", "treat"],
        cq_values=[18, 22, 18, 21, 18, 23],
    )
    res = compute_delta_delta_ct(df, ["gapdh"], reference_group="ctrl")
    r = res["gapdh"]
    ctrl_ddct_mean = r[r["Group"] == "ctrl"]["ddCt"].mean()
    assert abs(ctrl_ddct_mean) < 1e-9


def test_ddct_reference_group_identical_samples_are_one():
    """When all ctrl samples have identical dCt, each Relative_Expr = 1."""
    df = _make_df(
        targets=["gapdh", "geneX"] * 4,
        samples=["s1", "s1", "s2", "s2", "s3", "s3", "s4", "s4"],
        groups=["ctrl"] * 4 + ["treat"] * 4,
        cq_values=[18, 22, 18, 22, 18, 23, 18, 24],  # s1 and s2 identical dCt=4
    )
    res = compute_delta_delta_ct(df, ["gapdh"], reference_group="ctrl")
    r = res["gapdh"]
    # s1 and s2 have dCt=4; ref_dCt = mean(4, 4) = 4 → ddCt = 0 → Expr = 1
    ctrl_expr = r[r["Group"] == "ctrl"]["Relative_Expr"]
    assert abs(ctrl_expr.mean() - 1.0) < 1e-9


def test_ddct_relative_values():
    """Verify specific ΔΔCt values against hand-calculated expectations."""
    # s1 (ctrl): dCt = 22 - 18 = 4  → ref_dCt = 4  → ddCt = 0 → Expr = 1
    # s2 (treat): dCt = 24 - 18 = 6 → ddCt = 6 - 4 = 2 → Expr = 2^(-2) = 0.25
    df = _make_df(
        targets=["gapdh", "geneX", "gapdh", "geneX"],
        samples=["s1", "s1", "s2", "s2"],
        groups=["ctrl", "ctrl", "treat", "treat"],
        cq_values=[18, 22, 18, 24],
    )
    res = compute_delta_delta_ct(df, ["gapdh"], reference_group="ctrl")
    r = res["gapdh"]
    by_sample = dict(zip(r["Sample"], r["Relative_Expr"]))
    assert abs(by_sample["s1"] - 1.0) < 1e-9
    assert abs(by_sample["s2"] - 0.25) < 1e-9


def test_ddct_batch_aware():
    """Each batch is normalised independently; both batches end up in output."""
    # Batch A: s1 (ctrl), s2 (treat)
    # Batch B: s3 (ctrl), s4 (treat)
    # Within each batch, ctrl anchors at 1.
    df = _make_df(
        targets=["gapdh", "geneX"] * 4,
        samples=["s1", "s1", "s2", "s2", "s3", "s3", "s4", "s4"],
        groups=["ctrl", "ctrl", "treat", "treat", "ctrl", "ctrl", "treat", "treat"],
        cq_values=[18, 22, 18, 25, 20, 24, 20, 27],
    )
    batches = {"s1": "A", "s2": "A", "s3": "B", "s4": "B"}
    res = compute_delta_delta_ct(
        df, ["gapdh"], reference_group="ctrl", sample_batches=batches
    )
    r = res["gapdh"]
    # All four samples must appear in the output
    assert set(r["Sample"].unique()) == {"s1", "s2", "s3", "s4"}
    # Ctrl samples in each batch → Relative_Expr = 1
    ctrl = r[r["Group"] == "ctrl"]
    assert (ctrl["Relative_Expr"] - 1.0).abs().max() < 1e-9
    # Treat samples: dCt difference = 3 → Expr = 2^(-3) = 0.125
    treat = r[r["Group"] == "treat"]
    assert (treat["Relative_Expr"] - 0.125).abs().max() < 1e-9


def test_ddct_batch_metadata_in_output():
    """Output DataFrame must carry Batch and Reference_Group columns."""
    df = _make_df(
        targets=["gapdh", "geneX"],
        samples=["s1", "s1"],
        groups=["ctrl", "ctrl"],
        cq_values=[18.0, 22.0],
    )
    res = compute_delta_delta_ct(
        df, ["gapdh"], reference_group="ctrl",
        sample_batches={"s1": "run1"}
    )
    r = res["gapdh"]
    assert "Batch" in r.columns
    assert "Reference_Group" in r.columns
    assert r["Batch"].iloc[0] == "run1"
    assert r["Reference_Group"].iloc[0] == "ctrl"


def test_ddct_invalid_reference_group_raises():
    df = _make_df(
        targets=["gapdh", "geneX"],
        samples=["s1", "s1"],
        groups=["ctrl", "ctrl"],
        cq_values=[18.0, 22.0],
    )
    with pytest.raises(ValueError, match="nonexistent"):
        compute_delta_delta_ct(df, ["gapdh"], reference_group="nonexistent")


def test_ddct_batch_missing_reference_raises():
    """A batch with no reference-group samples must raise ValueError."""
    df = _make_df(
        targets=["gapdh", "geneX"] * 2,
        samples=["s1", "s1", "s2", "s2"],
        groups=["ctrl", "ctrl", "treat", "treat"],
        cq_values=[18, 22, 18, 25],
    )
    # s1 (ctrl) in batch A, s2 (treat) in batch B — batch B has no ctrl
    batches = {"s1": "A", "s2": "B"}
    with pytest.raises(ValueError, match="batch"):
        compute_delta_delta_ct(
            df, ["gapdh"], reference_group="ctrl", sample_batches=batches
        )


def test_ddct_multiple_ref_genes():
    df = _make_df(
        targets=["gapdh", "actb", "geneX"] * 2,
        samples=["s1"] * 3 + ["s2"] * 3,
        groups=["ctrl"] * 3 + ["treat"] * 3,
        cq_values=[18, 19, 22, 18, 19, 24],
    )
    res = compute_delta_delta_ct(df, ["gapdh", "actb"], reference_group="ctrl")
    assert set(res.keys()) == {"gapdh", "actb"}
    for ref, r in res.items():
        ctrl_val = r[r["Sample"] == "s1"]["Relative_Expr"].iloc[0]
        assert abs(ctrl_val - 1.0) < 1e-9
