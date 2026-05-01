"""Tests for core/summary.py — sort, blocks, and summarisation helpers."""

import pandas as pd

from qpcr_analyzer.core.summary import (
    build_blocks,
    sort_wells,
    summarize_dataset,
    target_order,
    well_sort_key,
)


def _make_flagged(rows):
    df = pd.DataFrame(rows)
    if "Outlier" not in df.columns:
        df["Outlier"] = False
    if "Replicates" not in df.columns:
        df["Replicates"] = 1
    return df


def test_well_sort_key_natural_order():
    wells = ["A10", "A2", "B1", "A1"]
    sorted_wells = sorted(wells, key=well_sort_key)
    assert sorted_wells == ["A1", "A2", "A10", "B1"]


def test_well_sort_key_non_conforming_last():
    wells = ["A1", "X", "A2"]
    sorted_wells = sorted(wells, key=well_sort_key)
    assert sorted_wells[:2] == ["A1", "A2"]
    assert sorted_wells[-1] == "X"


def test_target_order_preserves_first_appearance():
    df = pd.DataFrame({"Target": ["geneB", "geneA", "geneB", "geneC"]})
    assert target_order(df) == ["geneB", "geneA", "geneC"]


def test_sort_wells_by_target_group_sample_well():
    df = pd.DataFrame(
        {
            "Well": ["A2", "A1", "B1", "A10"],
            "Target": ["t1", "t1", "t2", "t1"],
            "Group": ["g1", "g1", "g1", "g1"],
            "Sample": ["s1", "s1", "s1", "s1"],
            "Cq": [1.0, 1.0, 1.0, 1.0],
        }
    )
    out = sort_wells(df)
    # t1 first (file order), within t1: A1, A2, A10 (natural)
    assert out["Target"].tolist() == ["t1", "t1", "t1", "t2"]
    assert out["Well"].tolist()[:3] == ["A1", "A2", "A10"]


def test_build_blocks_floats_excluded_to_top():
    flagged = _make_flagged(
        [
            {"Well": "A1", "Target": "geneX", "Group": "ctrl", "Sample": "s1", "Cq": 20.0},
            {"Well": "A2", "Target": "geneX", "Group": "ctrl", "Sample": "s1", "Cq": 20.1},
            {"Well": "B1", "Target": "geneX", "Group": "ctrl", "Sample": "s2", "Cq": 21.0},
            {"Well": "B2", "Target": "geneX", "Group": "ctrl", "Sample": "s2", "Cq": 21.1},
        ]
    )
    excluded = {"A1"}
    blocks = build_blocks(flagged, excluded)
    assert len(blocks) == 2
    assert blocks[0]["sample"] == "s1"
    assert blocks[0]["has_exclusion"] is True
    assert blocks[1]["sample"] == "s2"
    assert blocks[1]["has_exclusion"] is False


def test_build_blocks_groups_replicates_together():
    flagged = _make_flagged(
        [
            {"Well": "A1", "Target": "geneX", "Group": "ctrl", "Sample": "s1", "Cq": 20.0},
            {"Well": "A2", "Target": "geneX", "Group": "ctrl", "Sample": "s1", "Cq": 20.1},
            {"Well": "A3", "Target": "geneX", "Group": "ctrl", "Sample": "s1", "Cq": 30.0},
        ]
    )
    excluded = {"A3"}
    blocks = build_blocks(flagged, excluded)
    assert len(blocks) == 1
    blk = blocks[0]
    assert len(blk["replicates"]) == 3
    assert blk["n_excluded"] == 1


def test_summarize_dataset_basic():
    df = pd.DataFrame(
        {
            "Well": ["A1", "A2", "A3", "A4"],
            "Target": ["geneA", "geneB", "geneA", "geneB"],
            "Sample": ["s1", "s1", "s2", "s2"],
            "Group": ["ctrl", "ctrl", "treat", "treat"],
            "Cq": [20.0, 22.0, float("nan"), 23.0],
        }
    )
    s = summarize_dataset(df, filename="test.csv")
    assert s["filename"] == "test.csv"
    assert s["n_wells"] == 4
    assert s["n_samples"] == 2
    assert s["samples"] == ["s1", "s2"]
    assert s["n_targets"] == 2
    assert s["targets"] == ["geneA", "geneB"]  # file order
    assert s["n_groups"] == 2
    assert s["n_na_wells"] == 1
    assert s["has_batch_column"] is False


def test_summarize_dataset_with_batch():
    df = pd.DataFrame(
        {
            "Well": ["A1", "A2"],
            "Target": ["g", "g"],
            "Sample": ["s1", "s2"],
            "Group": ["c", "c"],
            "Cq": [20.0, 20.5],
            "Batch": ["b1", "b2"],
        }
    )
    s = summarize_dataset(df)
    assert s["has_batch_column"] is True
    assert s["n_batches"] == 2
    assert set(s["batches"]) == {"b1", "b2"}
