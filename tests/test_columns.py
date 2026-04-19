import pandas as pd

from qpcr_analyzer.core.columns import (
    ColumnMapping,
    apply_mapping,
    detect_columns,
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
    df = pd.DataFrame(
        {"W": ["A1"], "T": ["g"], "S": ["s"], "C": [20.0]}
    )
    m = ColumnMapping(
        assignments={
            "well": "W",
            "target": "T",
            "sample": "S",
            "cq": "C",
            "group": None,
        },
        confidence={},
    )
    out = apply_mapping(df, m)
    assert set(out.columns) >= {"Well", "Target", "Sample", "Cq", "Group"}
    assert (out["Group"] == "all").all()
