"""Tests for the file loader.

Cover the messy real-world inputs that motivated the auto-detection in
``read_table``: prepended metadata blocks with various comment prefixes,
``Undetermined`` / ``N/A`` Cq markers, partially-annotated empty wells.
"""

from __future__ import annotations

import io as _io

import pandas as pd
import pytest

from qpcr_analyzer.core.columns import apply_mapping, detect_columns
from qpcr_analyzer.core.io import read_table


def _bytes(text: str) -> _io.BytesIO:
    return _io.BytesIO(text.encode("utf-8"))


def test_csv_plain_header_unchanged():
    csv = "Well,Target,Sample,Cq\nA1,g,s,23.5\nA2,g,s,23.7\n"
    df = read_table(_bytes(csv), "x.csv")
    assert list(df.columns) == ["Well", "Target", "Sample", "Cq"]
    assert len(df) == 2


def test_csv_skips_hash_comments():
    csv = (
        "# Run Name: Plate1\n"
        "# Date: 2026-01-01\n"
        "Well,Target,Sample,Cq\n"
        "A1,g,s,23.5\n"
    )
    df = read_table(_bytes(csv), "x.csv")
    assert list(df.columns) == ["Well", "Target", "Sample", "Cq"]
    assert len(df) == 1


def test_csv_skips_mixed_prefix_comments():
    """Metadata can use *, ;, //, ! or no prefix — header is still found."""
    csv = (
        "# Comment with hash\n"
        "* Block Type,384-Well Block\n"
        "; Other comment\n"
        "// C-style comment\n"
        "Run Name: Plate1\n"
        "\n"
        "Well,Target,Sample,Cq\n"
        "A1,g,s,23.5\n"
    )
    df = read_table(_bytes(csv), "x.csv")
    assert list(df.columns) == ["Well", "Target", "Sample", "Cq"]
    assert len(df) == 1


def test_csv_handles_undetermined_and_whitespace():
    csv = (
        "Well,Target,Sample,Cq\n"
        "A1,g,s, 23.5 \n"
        "A2,g,s,Undetermined\n"
        "A3,g,s,N/A\n"
        "A4,g,s,\n"
        "A5,g,s,23.8\n"
    )
    df = read_table(_bytes(csv), "x.csv")
    m = detect_columns(df)
    std = apply_mapping(df, m)
    assert std["Cq"].dtype.kind == "f"
    assert std["Cq"].notna().sum() == 2
    assert std.loc[std["Well"] == "A1", "Cq"].iloc[0] == pytest.approx(23.5)


def test_apply_mapping_only_drops_fully_empty_rows():
    """apply_mapping only silent-drops rows where Cq, Sample, AND Target are all blank.

    Rows with partial missingness (e.g. only Sample blank) are KEPT — they
    are surfaced for review via :func:`find_incomplete_rows` and the UI's
    data-completeness panel rather than dropped silently.
    """
    df = pd.DataFrame(
        {
            "Well": ["A1", "A2", "A3", "A4", "A5", "A6"],
            "Target": ["g", "g", "", "g", None, ""],
            "Sample": ["s", "", "s", "s", "s", ""],
            "Cq":    [23.5, 24.0, 25.0, None, 26.0, None],  # A6: all three blank
        }
    )
    m = detect_columns(df)
    std = apply_mapping(df, m)
    # A6 dropped silently (Cq+Target+Sample all empty). Every other row kept.
    assert list(std["Well"]) == ["A1", "A2", "A3", "A4", "A5"]


def test_find_incomplete_rows_flags_partial_missingness():
    from qpcr_analyzer.core import find_incomplete_rows
    df = pd.DataFrame(
        {
            "Well": ["A1", "A2", "A3", "A4"],
            "Target": ["g", "g", "", "g"],
            "Sample": ["s", "", "s", "s"],
            "Cq": [23.5, 24.0, 25.0, 26.0],
            "Group": ["a"] * 4,
        }
    )
    incomplete = find_incomplete_rows(df)
    assert list(incomplete["Well"]) == ["A2", "A3"]
    assert list(incomplete["Missing"]) == ["Sample", "Target"]


def test_find_incomplete_rows_ignores_missing_cq():
    """A row missing only Cq (Undetermined well) is NOT flagged as incomplete."""
    from qpcr_analyzer.core import find_incomplete_rows
    df = pd.DataFrame(
        {
            "Well": ["A1", "A2"],
            "Target": ["g", "g"],
            "Sample": ["s", "s"],
            "Cq": [23.5, None],
            "Group": ["a", "a"],
        }
    )
    assert find_incomplete_rows(df).empty


def test_annotation_coverage_reports_per_field():
    from qpcr_analyzer.core import annotation_coverage
    df = pd.DataFrame(
        {
            "Well": ["A1", "A2", ""],
            "Target": ["g", "g", "g"],
            "Sample": ["s", "", "s"],
            "Cq": [23.5, 24.0, None],
        }
    )
    cov = annotation_coverage(df)
    assert cov["Well"] == (2, 3)
    assert cov["Target"] == (3, 3)
    assert cov["Sample"] == (2, 3)
    assert cov["Cq"] == (2, 3)


def test_drop_incomplete_rows_round_trip():
    from qpcr_analyzer.core import drop_incomplete_rows
    df = pd.DataFrame(
        {
            "Well": ["A1", "A2", "A3"],
            "Target": ["g", "", "g"],
            "Sample": ["s", "s", "s"],
            "Cq": [23.5, 24.0, 25.0],
            "Group": ["a"] * 3,
        }
    )
    cleaned = drop_incomplete_rows(df)
    assert list(cleaned["Well"]) == ["A1", "A3"]


def test_apply_mapping_strips_whitespace_in_labels():
    df = pd.DataFrame(
        {
            "Well": [" A1 ", "A2"],
            "Target": ["gene1 ", " gene1"],
            "Sample": ["sx", "sx "],
            "Cq": [23.5, 23.7],
        }
    )
    m = detect_columns(df)
    std = apply_mapping(df, m)
    assert list(std["Well"]) == ["A1", "A2"]
    assert list(std["Target"]) == ["gene1", "gene1"]
    assert list(std["Sample"]) == ["sx", "sx"]


def test_tsv_with_metadata_header():
    tsv = (
        "# File Name\tplate.eds\n"
        "* Date\t2026-01-01\n"
        "\n"
        "Well\tTarget\tSample\tCq\n"
        "A1\tg\ts\t23.5\n"
    )
    df = read_table(_bytes(tsv), "x.tsv")
    assert list(df.columns) == ["Well", "Target", "Sample", "Cq"]
    assert len(df) == 1


def test_txt_sniffs_separator():
    txt = (
        "# Run Name: Plate1\n"
        "Well\tTarget\tSample\tCq\n"
        "A1\tg\ts\t23.5\n"
        "A2\tg\ts\t23.7\n"
    )
    df = read_table(_bytes(txt), "x.txt")
    assert list(df.columns) == ["Well", "Target", "Sample", "Cq"]
    assert len(df) == 2


def test_excel_with_metadata_rows(tmp_path):
    path = tmp_path / "plate.xlsx"
    rows = [
        ["# Run Name: Plate1", None, None, None],
        ["Date", "2026-01-01", None, None],
        [None, None, None, None],
        ["Well", "Target", "Sample", "Cq"],
        ["A1", "g", "s", 23.5],
        ["A2", "g", "s", "Undetermined"],
        ["A3", "g", None, 24.0],
    ]
    pd.DataFrame(rows).to_excel(path, index=False, header=False)

    df = read_table(str(path), str(path))
    assert list(df.columns) == ["Well", "Target", "Sample", "Cq"]
    m = detect_columns(df)
    std = apply_mapping(df, m)
    # apply_mapping keeps partially-incomplete rows so the UI can flag them.
    # A3 (missing Sample) is now retained until UI-level cleanup excludes it.
    assert list(std["Well"]) == ["A1", "A2", "A3"]
    assert std["Cq"].iloc[0] == pytest.approx(23.5)
    assert pd.isna(std["Cq"].iloc[1])
    from qpcr_analyzer.core import find_incomplete_rows
    incomplete = find_incomplete_rows(std)
    assert list(incomplete["Well"]) == ["A3"]
    assert list(incomplete["Missing"]) == ["Sample"]


def test_csv_bom_stripped():
    csv = "﻿Well,Target,Sample,Cq\nA1,g,s,23.5\n"
    df = read_table(_bytes(csv), "x.csv")
    assert list(df.columns) == ["Well", "Target", "Sample", "Cq"]


def test_unsupported_extension_rejected():
    with pytest.raises(ValueError, match="Unsupported file type"):
        read_table(_bytes("foo"), "x.json")
