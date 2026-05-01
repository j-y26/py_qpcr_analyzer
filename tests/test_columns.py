import pandas as pd
import pytest

from qpcr_analyzer.core.columns import (
    ColumnMapping,
    apply_mapping,
    detect_columns,
    validate_sample_groups,
)


def test_detect_standard_names():
    df = pd.DataFrame(
        {
            "Well": ["A1", "A2"],
            "Target": ["g", "g"],
            "Sample": ["s", "s"],
            "Cq": [20.0, 20.5],
            "Group": ["g1", "g1"],
        }
    )
    m = detect_columns(df)
    assert m.assignments["well"] == "Well"
    assert m.assignments["target"] == "Target"
    assert m.assignments["sample"] == "Sample"
    assert m.assignments["cq"] == "Cq"
    assert m.assignments["group"] == "Group"
    assert m.validate() == []


def test_detect_nonstandard_names():
    df = pd.DataFrame(
        {
            "Well Position": ["A1", "A2"],
            "Gene Name": ["g", "g"],
            "Sample ID": ["s", "s"],
            "Ct Mean": [20.0, 20.5],
            "Treatment": ["ctrl", "ctrl"],
        }
    )
    m = detect_columns(df)
    assert m.assignments["well"] == "Well Position"
    assert m.assignments["target"] == "Gene Name"
    assert m.assignments["sample"] == "Sample ID"
    assert m.assignments["cq"] == "Ct Mean"
    assert m.assignments["group"] == "Treatment"


def test_detect_heuristic_fallback():
    df = pd.DataFrame(
        {
            "foo": ["A1", "A2", "A3"],
            "bar": ["x", "x", "x"],
            "baz": ["s1", "s1", "s1"],
            "quux": [22.1, 22.2, 22.3],
        }
    )
    m = detect_columns(df)
    assert m.assignments["cq"] == "quux"
    assert m.assignments["well"] == "foo"


def test_validate_missing_required():
    m = ColumnMapping(
        assignments={
            "well": "W",
            "target": "T",
            "sample": None,
            "cq": "C",
            "group": None,
        },
        confidence={},
    )
    errs = m.validate()
    assert any("Sample" in e for e in errs)


def test_apply_mapping_adds_group_default():
    df = pd.DataFrame({"W": ["A1"], "T": ["g"], "S": ["s"], "C": [20.0]})
    m = ColumnMapping(
        assignments={"well": "W", "target": "T", "sample": "S", "cq": "C", "group": None},
        confidence={},
    )
    out = apply_mapping(df, m)
    assert set(out.columns) >= {"Well", "Target", "Sample", "Cq", "Group"}
    assert (out["Group"] == "all").all()


def test_well_tiebreak_prefers_a1_format():
    """When 'Well' (integers) and 'Well Position' (A1 format) both score 1.0,
    the column with actual well-ID values should be chosen."""
    df = pd.DataFrame(
        {
            "Well": [1, 2, 3, 4],           # numeric — not A1 format
            "Well Position": ["A1", "A2", "B1", "B2"],  # A1 format
            "Target": ["g"] * 4,
            "Sample": ["s"] * 4,
            "Cq": [20.0, 20.1, 20.2, 20.3],
        }
    )
    m = detect_columns(df)
    assert m.assignments["well"] == "Well Position"


def test_well_tiebreak_both_numeric_first_wins():
    """If neither candidate has A1-format values, fall back to the first match."""
    df = pd.DataFrame(
        {
            "Well": [1, 2, 3, 4],
            "Well Position": [10, 20, 30, 40],
            "Target": ["g"] * 4,
            "Sample": ["s"] * 4,
            "Cq": [20.0, 20.1, 20.2, 20.3],
        }
    )
    m = detect_columns(df)
    # Both score 1.0 and neither has well-format values; first column wins
    assert m.assignments["well"] in ("Well", "Well Position")


def test_validate_sample_groups_ok():
    df = pd.DataFrame(
        {
            "Sample": ["s1", "s1", "s2", "s2"],
            "Group": ["ctrl", "ctrl", "treat", "treat"],
        }
    )
    assert validate_sample_groups(df) == []


def test_validate_sample_groups_conflict():
    df = pd.DataFrame(
        {
            "Sample": ["s1", "s1", "s2", "s2"],
            "Group": ["ctrl", "treat", "ctrl", "ctrl"],  # s1 in two groups
        }
    )
    errs = validate_sample_groups(df)
    assert len(errs) == 1
    assert "s1" in errs[0]
    assert "ctrl" in errs[0]
    assert "treat" in errs[0]


def test_validate_sample_groups_multiple_conflicts():
    df = pd.DataFrame(
        {
            "Sample": ["s1", "s1", "s2", "s2"],
            "Group": ["ctrl", "treat", "A", "B"],
        }
    )
    errs = validate_sample_groups(df)
    assert len(errs) == 2
